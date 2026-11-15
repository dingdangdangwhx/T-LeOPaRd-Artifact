/* Exact R_q equations, with optional small-modulus quotient range proofs.
 * Included by the NIZK2 frontend after its parameters and dimensions. */
#include "../../src/lazer/lnp-tbox.h"
#include <errno.h>
#include <sys/random.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef NIZK2_PUBLIC_PREP_FAST
#define NIZK2_PUBLIC_PREP_FAST 1
#endif
#ifndef NIZK2_PREHASH_L2
#define NIZK2_PREHASH_L2 1
#endif

/* Only for selectors constructed below from the same compiled parameters.
 * Never accept a statement hash from the peer. Call before setting a witness. */
static void mixed_fixed_l2(_lnp_state_t state, polymat_ptr norms[]) {
#if NIZK2_PREHASH_L2 && defined(NIZK2_FIXED_L2_HASH_AVAILABLE)
  if(state->params != params->tbox_params || state->statement_l2_set) abort();
  for(unsigned int i=0;i<state->params->Z;i++) {
    state->Es[i]=norms[i];
    state->Em[i]=NULL;
    state->v[i]=NULL;
  }
  memcpy(state->hash_l2,nizk2_fixed_l2_hash,32);
  state->statement_l2_set=1;
#else
  (void)state; (void)norms;
#endif
}

/* This early exit is only for public statement entries, not secret witnesses. */
static int mixed_public_zero(poly_t value, poly_t zero) {
#if NIZK2_PUBLIC_PREP_FAST
  intvec_srcptr c=poly_get_coeffvec(value);
  (void)zero;
  for(unsigned int i=0;i<c->nelems;i++) {
    int_srcptr x=intvec_get_elem_src(c,i);
    for(unsigned int j=0;j<x->nlimbs;j++) if(x->limbs[j]) return 0;
  }
  return 1;
#else
  return poly_eq(value,zero);
#endif
}

static double mixed_seconds(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec + 1e-9*(double)t.tv_nsec;
}

typedef struct {
  polymat_t quotient;
  polyvec_t constant;
  polymat_ptr *norms;
  spolymat_ptr *quadratic;
  spolyvec_ptr *linear;
  poly_ptr *rhs;
  unsigned int offset, count;
} nizk2_mixed_statement;

static void nizk2_private_seed(uint8_t seed[32]) {
  for(size_t off=0;off<32;) {
    ssize_t got=getrandom(seed+off,32-off,0);
    if(got<0 && errno==EINTR) continue;
    if(got<=0) { perror("NIZK2 private randomness"); abort(); }
    off+=(size_t)got;
  }
}

static void mixed_statement_init(nizk2_mixed_statement *st,
                                  polymat_t A, polyvec_t t) {
  lnp_tbox_params_srcptr p = params->tbox_params;
  polyring_srcptr ring = p->tbox->ring;
  const unsigned int nw = p->tbox->m1 - p->Z;
  const unsigned int nm = p->tbox->l;
  const unsigned int nt = nw + nm;
  const unsigned int nq = p->nprime;
  const unsigned int lambda = p->quad_eval->lambda;
  const unsigned int mask_slots=(256/D)*((p->nex>0)+(p->nprime>0));
  polymat_t top;
  polyvec_t top_t;
  poly_t zero;
  INT_T(inv, ring->q->nlimbs);
  INT_T(qc, ring->q->nlimbs);
  memset(st, 0, sizeof(*st));
  if(nt != WIT_ROWS || polymat_get_ncols(A) != nt
     || polymat_get_nrows(A) != REL_ROWS || polyvec_get_nelems(t) != REL_ROWS
     || params->same_modulus != 1 || params->k != 1
     || params->ns1_indices != nw || params->nm_indices != nm
     || params->Ps_nrows != 0 || ring->d != D
     || p->tbox->lext != mask_slots+1+lambda/2+1
     || p->quad_eval->quad_eval->l != nm+mask_slots+1
     || p->quad_eval->quad_many->l != nm+mask_slots+1+lambda/2
     || (nq && nm)
     || nq != (NIZK2_SMALL_COMMITMENT ? COM_ROWS : 0)) {
    fputs("NIZK2 relation/parameter mismatch; regenerate the parameter headers\n",stderr);
    abort();
  }
  unsigned int *column = calloc(nt,sizeof(*column));
  unsigned char *seen = calloc(nt,sizeof(*seen));
  if(!column || !seen) abort();
  for(unsigned int i=0;i<nt;i++) {
    unsigned int j = i<nw ? params->s1_indices[i] : params->m_indices[i-nw];
    if(j>=nt || seen[j] || (!nm && j!=i)) {
      fputs("NIZK2 invalid witness partition\n",stderr);
      abort();
    }
    seen[j]=1;
    column[i]=j;
  }
  free(seen);
  st->offset = QUADEQ_INPUT_OFF;
  st->count = polymat_get_nrows(A) - nq;
  st->norms = calloc(p->Z, sizeof(*st->norms));
  st->quadratic = calloc(st->offset + st->count, sizeof(*st->quadratic));
  st->linear = calloc(st->offset + st->count, sizeof(*st->linear));
  st->rhs = calloc(st->offset + st->count, sizeof(*st->rhs));
  if(!st->norms || !st->quadratic || !st->linear || !st->rhs) abort();
  poly_alloc(zero,ring);
  poly_set_zero(zero);
  if(nq) {
    polymat_alloc(st->quotient, ring, nq, nw);
    polyvec_alloc(st->constant, ring, nq);
    polymat_get_submat(top, A, 0, 0, nq, nw, 1, 1);
    polyvec_get_subvec(top_t, t, 0, nq, 1);
    int_set_i64(qc, NIZK2_COM_Q);
    int_invmod(inv, qc, ring->q);
    polymat_set_zero(st->quotient);
    /* Only the public nonzero blocks need quotient scaling and reduction. */
    for(unsigned int i=0;i<nq;i++)
      for(unsigned int j=0;j<nw;j++) {
        poly_ptr src=polymat_get_elem(top,i,j);
        if(mixed_public_zero(src,zero)) continue;
        poly_ptr dst=polymat_get_elem(st->quotient,i,j);
        poly_scale(dst,inv,src);
        poly_mod(dst,dst);
      }
    polyvec_scale(st->constant, inv, top_t);
    polyvec_mod(st->constant, st->constant);
  }
  for(unsigned int i=0;i<p->Z;i++) {
    st->norms[i] = malloc(sizeof(polymat_t));
    if(!st->norms[i]) abort();
    polymat_alloc(st->norms[i], ring, params->Es_nrows[i], nw);
    polymat_set_zero(st->norms[i]);
    for(unsigned int j=0;j<params->Es_nrows[i];j++)
      poly_set_one(polymat_get_elem(st->norms[i],j,params->Es[i][j]));
  }
#if defined(_OPENMP) && NIZK2_PUBLIC_PREP_FAST
#pragma omp parallel for schedule(static) if(!omp_in_parallel() && st->count>=32)
#endif
  for(unsigned int i=0;i<st->count;i++) {
    const unsigned int idx = st->offset+i;
    st->linear[idx] = malloc(sizeof(spolyvec_t));
    st->rhs[idx] = malloc(sizeof(poly_t));
    if(!st->linear[idx] || !st->rhs[idx]) abort();
    spolyvec_alloc(st->linear[idx], ring, 2*nt, nt);
    for(unsigned int j=0;j<nt;j++) {
      poly_ptr src = polymat_get_elem(A,nq+i,column[j]);
      if(mixed_public_zero(src,zero)) continue;
      /* The LNP toolbox orders (s_0,auto(s_0),s_1,auto(s_1),...). */
      poly_set(spolyvec_insert_elem(st->linear[idx],2*j),src);
    }
    spolyvec_sort(st->linear[idx]);
    poly_alloc(st->rhs[idx],ring);
    poly_set(st->rhs[idx],polyvec_get_elem(t,nq+i));
  }
  poly_free(zero);
  free(column);
}

static void mixed_statement_clear(nizk2_mixed_statement *st) {
  for(unsigned int i=0;i<params->tbox_params->Z;i++) {
    polymat_free(st->norms[i]); free(st->norms[i]);
  }
  for(unsigned int i=st->offset;i<st->offset+st->count;i++) {
    spolyvec_free(st->linear[i]); free(st->linear[i]);
    poly_free(st->rhs[i]); free(st->rhs[i]);
  }
  free(st->norms); free(st->quadratic); free(st->linear); free(st->rhs);
  if(params->tbox_params->nprime) {
    polymat_free(st->quotient); polyvec_free(st->constant);
  }
}

static size_t mixed_prover_coins(uint8_t *proof, polyvec_t s, polymat_t A,
                                 polyvec_t t, const uint8_t pp[32],
                                 const uint8_t coins[32]) {
  lnp_prover_state_t p;
  polyvec_t short_witness, message;
  nizk2_mixed_statement st;
  size_t len=0;
  const int profile=NIZK2_PROFILE_TIMERS || getenv("NIZK2_MIXED_PROFILE")!=NULL;
  double ts[7]={0};
  if(profile) ts[0]=mixed_seconds();
  mixed_statement_init(&st,A,t);
  if(profile) ts[1]=mixed_seconds();
  lnp_prover_init(p,pp,params->tbox_params);
  if(profile) ts[2]=mixed_seconds();
#if NIZK2_PREHASH_L2 && defined(NIZK2_FIXED_L2_HASH_AVAILABLE)
  if(p->witness_set) abort();
  mixed_fixed_l2(p->state,st.norms);
  if(getenv("NIZK2_CHECK_FIXED_L2_HASH")) {
    lnp_prover_set_statement_l2(p,st.norms,NULL,NULL);
    if(memcmp(p->state->hash_l2,nizk2_fixed_l2_hash,32)) abort();
  }
#else
  lnp_prover_set_statement_l2(p,st.norms,NULL,NULL);
#endif
  if(params->tbox_params->nprime)
    lnp_prover_set_statement_arp(p,st.quotient,NULL,st.constant);
  lnp_prover_set_statement_quadeqs(p,st.quadratic,st.linear,st.rhs,st.count);
  if(profile) ts[3]=mixed_seconds();
  const unsigned int ns=params->ns1_indices, nm=params->nm_indices;
  if(nm) {
    polyvec_alloc(short_witness,polyvec_get_ring(s),ns);
    polyvec_alloc(message,polyvec_get_ring(s),nm);
    for(unsigned int i=0;i<ns;i++)
      poly_set(polyvec_get_elem(short_witness,i),polyvec_get_elem(s,params->s1_indices[i]));
    for(unsigned int i=0;i<nm;i++)
      poly_set(polyvec_get_elem(message,i),polyvec_get_elem(s,params->m_indices[i]));
    lnp_prover_set_witness(p,short_witness,message);
    polyvec_free(short_witness);
    polyvec_free(message);
  } else {
    lnp_prover_set_witness(p,s,NULL);
  }
  if(profile) ts[4]=mixed_seconds();
  lnp_prover_prove(p,proof,&len,coins);
  if(profile) ts[5]=mixed_seconds();
  lnp_prover_clear(p); mixed_statement_clear(&st);
  if(profile) {
    ts[6]=mixed_seconds();
    printf("mixed prove: prepare=%.6f init=%.6f statement=%.6f witness=%.6f core=%.6f clear=%.6f total=%.6f\n",
           ts[1]-ts[0],ts[2]-ts[1],ts[3]-ts[2],ts[4]-ts[3],ts[5]-ts[4],ts[6]-ts[5],ts[6]-ts[0]);
    print_stopwatch_lnp_prover_prove(2);
  }
  return len;
}

static size_t mixed_prover(uint8_t *proof, polyvec_t s, polymat_t A,
                           polyvec_t t, const uint8_t pp[32]) {
  uint8_t coins[32];
  nizk2_private_seed(coins);
  return mixed_prover_coins(proof,s,A,t,pp,coins);
}

static int mixed_verifier(const uint8_t *proof, polymat_t A, polyvec_t t,
                           const uint8_t pp[32]) {
  lnp_verifier_state_t v;
  nizk2_mixed_statement st;
  const int profile=NIZK2_PROFILE_TIMERS || getenv("NIZK2_MIXED_PROFILE")!=NULL;
  double ts[6]={0};
  if(profile) ts[0]=mixed_seconds();
  mixed_statement_init(&st,A,t);
  if(profile) ts[1]=mixed_seconds();
  lnp_verifier_init(v,pp,params->tbox_params);
  if(profile) ts[2]=mixed_seconds();
#if NIZK2_PREHASH_L2 && defined(NIZK2_FIXED_L2_HASH_AVAILABLE)
  mixed_fixed_l2(v->state,st.norms);
  if(getenv("NIZK2_CHECK_FIXED_L2_HASH")) {
    lnp_verifier_set_statement_l2(v,st.norms,NULL,NULL);
    if(memcmp(v->state->hash_l2,nizk2_fixed_l2_hash,32)) abort();
  }
#else
  lnp_verifier_set_statement_l2(v,st.norms,NULL,NULL);
#endif
  if(params->tbox_params->nprime)
    lnp_verifier_set_statement_arp(v,st.quotient,NULL,st.constant);
  lnp_verifier_set_statement_quadeqs(v,st.quadratic,st.linear,st.rhs,st.count);
  if(profile) ts[3]=mixed_seconds();
  int ok=lnp_verifier_verify(v,proof,NULL);
  if(profile) ts[4]=mixed_seconds();
  lnp_verifier_clear(v); mixed_statement_clear(&st);
  if(profile) {
    ts[5]=mixed_seconds();
    printf("mixed verify: prepare=%.6f init=%.6f statement=%.6f core=%.6f clear=%.6f total=%.6f\n",
           ts[1]-ts[0],ts[2]-ts[1],ts[3]-ts[2],ts[4]-ts[3],ts[5]-ts[4],ts[5]-ts[0]);
    print_stopwatch_lnp_verifier_verify(2);
  }
  return ok;
}
