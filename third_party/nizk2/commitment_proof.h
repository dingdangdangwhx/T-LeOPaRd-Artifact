/* Affine opening residuals retain the original joint L2 bound. */
#include "../../src/lazer/lnp-tbox.h"
#include <errno.h>
#include <sys/random.h>
#include <time.h>
#include "../../src/lazer/threading.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef NIZK2_KEY_NORM
#define NIZK2_KEY_NORM 1
#endif

enum {
  LAYOUT_OPEN = NIZK2_IMPLICIT_OPENING == 1 ? COM_RANK :
    (NIZK2_IMPLICIT_OPENING == 2 ? COM_RANK+COM_MESSAGE_ROWS : WIT_KEY_OFF),
  LAYOUT_ERROR = LAYOUT_OPEN,
  LAYOUT_DIGIT = LAYOUT_ERROR + M + N,
  LAYOUT_KEY = LAYOUT_DIGIT + 2*H,
  LAYOUT_SHORT = LAYOUT_KEY + N*NIZK2_KEY_DIGITS,
  LAYOUT_MESSAGE = NIZK2_KEY_DIGITS ? 0 : N,
  LAYOUT_TOTAL = LAYOUT_SHORT + LAYOUT_MESSAGE,
  LAYOUT_SKIP = NIZK2_IMPLICIT_OPENING == 1 ? COM_ROWS :
    (NIZK2_IMPLICIT_OPENING == 2 ? COM_RANK : 0),
  LAYOUT_NORMS = NIZK2_KEY_DIGITS && NIZK2_KEY_NORM ? 5 : 4
};

static void nizk2_private_seed(uint8_t seed[32]) {
  for(size_t off=0;off<32;) {
    ssize_t n=getrandom(seed+off,32-off,0);
    if(n<0 && errno==EINTR) continue;
    if(n<=0) { perror("NIZK2 private randomness"); abort(); }
    off+=(size_t)n;
  }
}

static double layout_seconds(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC,&t);
  return t.tv_sec+1e-9*t.tv_nsec;
}

static int layout_zero(poly_t p) {
  intvec_srcptr c=poly_get_coeffvec(p);
  for(unsigned int i=0;i<c->nelems;i++) {
    int_srcptr x=intvec_get_elem_src(c,i);
    for(unsigned int j=0;j<x->nlimbs;j++) if(x->limbs[j]) return 0;
  }
  return 1;
}

typedef struct {
  polymat_ptr Es[LAYOUT_NORMS], Em[LAYOUT_NORMS];
  polyvec_ptr v[LAYOUT_NORMS];
  spolymat_ptr *quadratic;
  spolyvec_ptr *linear;
  poly_ptr *constant;
  unsigned int count;
} layout_statement;

static polymat_ptr layout_matrix(polyring_srcptr ring, unsigned int rows,
                                 unsigned int cols) {
  polymat_ptr p=malloc(sizeof(polymat_t));
  if(!p) abort();
  polymat_alloc(p,ring,rows,cols);
  polymat_set_zero(p);
  return p;
}

static unsigned int layout_original_column(unsigned int col) {
  if(col<LAYOUT_OPEN)
    return NIZK2_IMPLICIT_OPENING == 2 && col>=COM_RANK ? col+COM_RANK : col;
  if(col<LAYOUT_KEY) return WIT_ES_OFF+col-LAYOUT_ERROR;
  return WIT_KEY_OFF+(col-LAYOUT_KEY)%N;
}

static poly_ptr layout_entry(polymat_t A, polymat_t transformed,
                             unsigned int row, unsigned int col) {
#if NIZK2_KEY_BITS > 0
  (void)A;
  return polymat_get_elem(transformed,row,col);
#else
  (void)transformed;
  return polymat_get_elem(A,row,layout_original_column(col));
#endif
}

static void layout_statement_init(layout_statement *st, polymat_t A, polyvec_t t) {
  const lnp_tbox_params_srcptr p=params->tbox_params;
  const unsigned int lambda=p->quad_eval->lambda;
  polyring_srcptr ring=p->tbox->ring;
  polymat_t transformed;
  memset(st,0,sizeof(*st));
  if(!NIZK2_DOUBLE_COMMITMENT || ring->d!=64 ||
     p->tbox->m1-p->Z!=LAYOUT_SHORT || p->tbox->l!=LAYOUT_MESSAGE ||
     p->Z!=LAYOUT_NORMS || p->nprime || p->nbin ||
     polymat_get_nrows(A)!=REL_ROWS || polymat_get_ncols(A)!=WIT_ROWS ||
     polyvec_get_nelems(t)!=REL_ROWS) abort();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int row=0;row<REL_ROWS;row++)
    for(unsigned int col=0;col<WIT_ROWS;col++)
      poly_fromcrt(polymat_get_elem(A,row,col));
#if NIZK2_KEY_BITS > 0
  polymat_alloc(transformed,ring,REL_ROWS,LAYOUT_TOTAL);
  polymat_set_zero(transformed);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int col=0;col<LAYOUT_TOTAL;col++) {
    INT_T(scale,ring->q->nlimbs);
    unsigned int src=layout_original_column(col);
    int_set_one(scale);
    if(NIZK2_KEY_DIGITS && col>=LAYOUT_KEY)
      int_lshift(scale,scale,((col-LAYOUT_KEY)/N)*NIZK2_KEY_BITS);
    for(unsigned int row=0;row<REL_ROWS;row++) {
      poly_ptr a=polymat_get_elem(A,row,src);
      poly_ptr b=polymat_get_elem(transformed,row,col);
      if(layout_zero(a)) continue;
      if(NIZK2_KEY_DIGITS && col>=LAYOUT_KEY) {
        poly_scale(b,scale,a);
        poly_mod(b,b);
      } else poly_set(b,a);
    }
  }
#endif
  const unsigned int rows[]={WIT_KEY_OFF,M+N,H,H,N*NIZK2_KEY_DIGITS};
  for(unsigned int g=0;g<LAYOUT_NORMS;g++) {
    if(p->n[g]!=rows[g]) abort();
    st->Es[g]=layout_matrix(ring,rows[g],LAYOUT_SHORT);
  }
  if(NIZK2_IMPLICIT_OPENING) {
    for(unsigned int i=0;i<COM_RANK;i++)
      poly_set_one(polymat_get_elem(st->Es[0],i,i));
    if(LAYOUT_MESSAGE && NIZK2_IMPLICIT_OPENING == 1)
      st->Em[0]=layout_matrix(ring,WIT_KEY_OFF,LAYOUT_MESSAGE);
    st->v[0]=malloc(sizeof(polyvec_t));
    if(!st->v[0]) abort();
    polyvec_alloc(st->v[0],ring,WIT_KEY_OFF);
    polyvec_set_zero(st->v[0]);
    /* Squaring -(A_rest*w+t) or +(A_rest*w+t) gives the same joint norm. */
    for(unsigned int row=0;row<LAYOUT_SKIP;row++) {
      for(unsigned int col=0;col<LAYOUT_SHORT;col++)
        poly_set(polymat_get_elem(st->Es[0],COM_RANK+row,col),
                 layout_entry(A,transformed,row,col));
      for(unsigned int col=0;st->Em[0] && col<LAYOUT_MESSAGE;col++)
        poly_set(polymat_get_elem(st->Em[0],COM_RANK+row,col),
                 layout_entry(A,transformed,row,LAYOUT_SHORT+col));
      poly_set(polyvec_get_elem(st->v[0],COM_RANK+row),polyvec_get_elem(t,row));
    }
    if(NIZK2_IMPLICIT_OPENING == 2)
      for(unsigned int j=0;j<COM_MESSAGE_ROWS;j++)
        poly_set_one(polymat_get_elem(st->Es[0],2*COM_RANK+j,COM_RANK+j));
  } else {
    for(unsigned int i=0;i<WIT_KEY_OFF;i++)
      poly_set_one(polymat_get_elem(st->Es[0],i,i));
  }
  const unsigned int offsets[]={0,LAYOUT_ERROR,LAYOUT_DIGIT,LAYOUT_DIGIT+H,LAYOUT_KEY};
  for(unsigned int g=1;g<LAYOUT_NORMS;g++)
    for(unsigned int j=0;j<rows[g];j++)
      poly_set_one(polymat_get_elem(st->Es[g],j,offsets[g]+j));
  st->count=REL_ROWS-LAYOUT_SKIP;
  const unsigned int sz=QUADEQ_INPUT_OFF+st->count;
  st->quadratic=calloc(sz,sizeof(*st->quadratic));
  st->linear=calloc(sz,sizeof(*st->linear));
  st->constant=calloc(sz,sizeof(*st->constant));
  if(!st->quadratic || !st->linear || !st->constant) abort();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int row=LAYOUT_SKIP;row<REL_ROWS;row++) {
    unsigned int idx=QUADEQ_INPUT_OFF+row-LAYOUT_SKIP;
    st->linear[idx]=malloc(sizeof(spolyvec_t));
    st->constant[idx]=malloc(sizeof(poly_t));
    if(!st->linear[idx] || !st->constant[idx]) abort();
    spolyvec_alloc(st->linear[idx],ring,2*LAYOUT_TOTAL,LAYOUT_TOTAL);
    for(unsigned int col=0;col<LAYOUT_TOTAL;col++) {
      poly_ptr a=layout_entry(A,transformed,row,col);
      if(!layout_zero(a)) poly_set(spolyvec_insert_elem(st->linear[idx],2*col),a);
    }
    spolyvec_sort(st->linear[idx]);
    poly_alloc(st->constant[idx],ring);
    poly_set(st->constant[idx],polyvec_get_elem(t,row));
  }
#if NIZK2_KEY_BITS > 0
  polymat_free(transformed);
#endif
}

static void layout_statement_clear(layout_statement *st) {
  const unsigned int lambda=params->tbox_params->quad_eval->lambda;
  for(unsigned int g=0;g<LAYOUT_NORMS;g++) {
    polymat_free(st->Es[g]); free(st->Es[g]);
    if(st->Em[g]) { polymat_free(st->Em[g]); free(st->Em[g]); }
    if(st->v[g]) { polyvec_free(st->v[g]); free(st->v[g]); }
  }
  for(unsigned int i=QUADEQ_INPUT_OFF;i<QUADEQ_INPUT_OFF+st->count;i++) {
    spolyvec_free(st->linear[i]); free(st->linear[i]);
    poly_free(st->constant[i]); free(st->constant[i]);
  }
  free(st->quadratic); free(st->linear); free(st->constant);
}

static void layout_witness(polyvec_t s, polyvec_t shortw, polyvec_t message) {
  for(unsigned int col=0;col<LAYOUT_KEY;col++)
    poly_set(polyvec_get_elem(shortw,col),polyvec_get_elem(s,layout_original_column(col)));
  if(!NIZK2_KEY_DIGITS) {
    for(unsigned int j=0;j<N;j++)
      poly_set(polyvec_get_elem(message,j),polyvec_get_elem(s,WIT_KEY_OFF+j));
    return;
  }
#if NIZK2_KEY_BITS > 0
  const __int128 base=(__int128)1<<NIZK2_KEY_BITS;
  for(unsigned int j=0;j<N;j++) {
    poly_ptr key=polyvec_get_elem(s,WIT_KEY_OFF+j);
    poly_fromcrt(key); poly_mod(key,key); poly_redc(key,key);
    for(unsigned int h=0;h<D;h++) {
      int_srcptr x=poly_get_coeff(key,h);
      if(x->nlimbs>2) abort();
      unsigned __int128 mag=x->limbs[0];
      if(x->nlimbs>1) mag|=(unsigned __int128)x->limbs[1]<<64;
      __int128 v=x->neg ? -(__int128)mag : (__int128)mag;
      for(unsigned int g=0;g<NIZK2_KEY_DIGITS;g++) {
        __int128 rem=v%base;
        if(rem>base/2) rem-=base;
        if(rem< -base/2) rem+=base;
        int_set_i64(poly_get_coeff(polyvec_get_elem(shortw,LAYOUT_KEY+g*N+j),h),
                    (int64_t)rem);
        v=(v-rem)/base;
      }
      if(v) abort();
    }
  }
#endif
}

static size_t mixed_prover_coins(uint8_t *proof, polyvec_t s, polymat_t A,
                                 polyvec_t t, const uint8_t pp[32],
                                 const uint8_t coins[32]) {
  layout_statement st;
  lnp_prover_state_t p;
  polyvec_t shortw,message;
  polyring_srcptr ring=params->tbox_params->tbox->ring;
  double start=layout_seconds();
  layout_statement_init(&st,A,t);
  double prepared=layout_seconds();
  lnp_prover_init(p,pp,params->tbox_params);
  lnp_prover_set_statement_l2(p,st.Es,st.Em,st.v);
  lnp_prover_set_statement_quadeqs(p,st.quadratic,st.linear,st.constant,st.count);
  polyvec_alloc(shortw,ring,LAYOUT_SHORT); polyvec_set_zero(shortw);
  if(LAYOUT_MESSAGE) polyvec_alloc(message,ring,LAYOUT_MESSAGE);
  layout_witness(s,shortw,message);
  if(getenv("NIZK2_LAYOUT_CHECK_NORMS")) {
    INT_T(norm,2*ring->q->nlimbs);
    for(unsigned int g=0;g<LAYOUT_NORMS;g++) {
      polyvec_t residual;
      polyvec_alloc(residual,ring,params->tbox_params->n[g]);
      if(st.v[g]) polyvec_set(residual,st.v[g]);
      else polyvec_set_zero(residual);
      polyvec_addmul(residual,st.Es[g],shortw,0);
      if(st.Em[g]) polyvec_addmul(residual,st.Em[g],message,0);
      polyvec_l2sqr(norm,residual);
      printf("norm group %u raw=",g); int_out_str(stdout,10,norm);
      polyvec_mod(residual,residual); polyvec_redc(residual,residual);
      polyvec_l2sqr(norm,residual);
      printf(" centered="); int_out_str(stdout,10,norm);
      printf(" limit="); int_out_str(stdout,10,params->tbox_params->l2Bsqr[g]);
      putchar('\n'); fflush(stdout);
      polyvec_free(residual);
    }
  }
  lnp_prover_set_witness(p,shortw,LAYOUT_MESSAGE ? message : NULL);
  double setup=layout_seconds();
  size_t len=0;
  lnp_prover_prove(p,proof,&len,coins);
  double core=layout_seconds();
  lnp_prover_clear(p);
  polyvec_free(shortw);
  if(LAYOUT_MESSAGE) polyvec_free(message);
  layout_statement_clear(&st);
  if(getenv("NIZK2_MIXED_PROFILE"))
    printf("layout prove: prepare=%.6f setup=%.6f core=%.6f total=%.6f\n",
           prepared-start,setup-prepared,core-setup,layout_seconds()-start);
  if(getenv("NIZK2_MIXED_PROFILE")) print_stopwatch_lnp_prover_prove(2);
  return len;
}

static size_t mixed_prover(uint8_t *proof, polyvec_t s, polymat_t A,
                           polyvec_t t, const uint8_t pp[32]) {
  uint8_t coins[32]; nizk2_private_seed(coins);
  return mixed_prover_coins(proof,s,A,t,pp,coins);
}

static int mixed_verifier(const uint8_t *proof, polymat_t A, polyvec_t t,
                          const uint8_t pp[32]) {
  layout_statement st;
  lnp_verifier_state_t v;
  double start=layout_seconds();
  layout_statement_init(&st,A,t);
  double prepared=layout_seconds();
  lnp_verifier_init(v,pp,params->tbox_params);
  lnp_verifier_set_statement_l2(v,st.Es,st.Em,st.v);
  lnp_verifier_set_statement_quadeqs(v,st.quadratic,st.linear,st.constant,st.count);
  double setup=layout_seconds();
  int ok=lnp_verifier_verify(v,proof,NULL);
  double core=layout_seconds();
  lnp_verifier_clear(v); layout_statement_clear(&st);
  if(getenv("NIZK2_MIXED_PROFILE"))
    printf("layout verify: prepare=%.6f setup=%.6f core=%.6f total=%.6f\n",
           prepared-start,setup-prepared,core-setup,layout_seconds()-start);
  if(getenv("NIZK2_MIXED_PROFILE")) print_stopwatch_lnp_verifier_verify(2);
  return ok;
}
