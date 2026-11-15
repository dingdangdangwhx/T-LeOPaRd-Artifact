#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "malloc.h"
#include "randombytes.h"
#include "fips202.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "chihuahua.h"
#include "pack.h"

#ifndef DEMO_L
#define DEMO_L 16
#endif
#ifndef SELECTOR_MAX_ATTEMPTS
#define SELECTOR_MAX_ATTEMPTS 64
#endif
#ifndef GADGET_BASE
#define GADGET_BASE 3
#endif
#if (GADGET_BASE < 3) || ((GADGET_BASE % 2) == 0)
#error "GADGET_BASE must be an odd integer >= 3"
#endif
#ifndef GADGET_LABEL
#if GADGET_BASE == 3
#define GADGET_LABEL "ternary"
#elif GADGET_BASE == 5
#define GADGET_LABEL "quinary"
#else
#define GADGET_LABEL "odd-base"
#endif
#endif
#ifndef DEMO_M
#define DEMO_M 24
#endif
#ifndef DEMO_KAPPA
#define DEMO_KAPPA 4
#endif
#ifndef SMALL_BLOCK
#define SMALL_BLOCK 4
#endif

enum {
  MAX_ATTEMPTS = SELECTOR_MAX_ATTEMPTS,
  GADGET_CENTER = (GADGET_BASE - 1)/2
};

static const uint8_t demo_fixed_seed[16] = {
  0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
  0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f
};

static size_t b_block(size_t t) {
  return t;
}

static size_t x_block(size_t t) {
  return DEMO_L + t;
}

static size_t z_block(size_t t) {
  return 2*DEMO_L + (t - 1);
}

static size_t y_block(size_t t) {
  return 3*DEMO_L - 1 + (t - 1);
}

static void absorb_polx(shake128incctx *shakectx, const polx *a) {
  __attribute__((aligned(16)))
  uint8_t hashbuf[N*QBYTES];
  polz t[1];

  polz_frompolx(t,a);
  polz_bitpack(hashbuf,t);
  shake128_inc_absorb(shakectx,hashbuf,sizeof(hashbuf));
}

static void set_const_coeff(polx *p, __int128 x) {
#if LOGQ > 64
  polz z = {0};
  polz_setcoeff_fromint128(&z,x,0);
  polz_topolx(p,&z);
#else
  polx_monomial(p,x,0);
#endif
}

static void poly_set_zero(poly *p) {
  polyvec_setzero(p,1);
}

static void poly_set_const(poly *p, int16_t x) {
  poly_set_zero(p);
  p->vec->c[0] = x;
}

static unsigned __int128 zz_to_u64(const zz *a) {
  unsigned __int128 r = 0;
  size_t i;

  for(i=0;i<L-1;i++)
    r |= (unsigned __int128)(uint16_t)a->limbs[i] << (14*i);
  r |= (unsigned __int128)(uint16_t)a->limbs[L-1] << (14*(L-1));
  return r;
}

static size_t gadget_len_base3(void) {
  unsigned __int128 q, bound, powb;
  size_t len = 0;

  q = ((unsigned __int128)1 << LOGQ) - QOFF;
  bound = 2*q - 2;
  powb = 1;
  while(powb <= bound) {
    powb *= GADGET_BASE;
    len++;
  }
  return len;
}

static __int128 pow3_u64(size_t e) {
  const unsigned __int128 q = ((unsigned __int128)1 << LOGQ) - QOFF;
  unsigned __int128 r = 1;

  while(e--)
    r = (r * GADGET_BASE) % q;
  if(r > q/2)
    return (__int128)r - (__int128)q;
  return (__int128)r;
}

static void gadget_decompose_balanced3(poly *out, const polx *in, size_t digits) {
  size_t i,j;
  polz z[1];
  unsigned __int128 corr = 0, p = 1;

  for(i=0;i<digits;i++) {
    corr += (unsigned __int128)GADGET_CENTER * p;
    p *= GADGET_BASE;
  }

  for(i=0;i<digits;i++)
    poly_set_zero(&out[i]);

  polz_frompolx(z,in);
  polz_caddq(z);
  for(j=0;j<N;j++) {
    zz coeff;
    unsigned __int128 w;

    polz_getcoeff(&coeff,z,j);
    w = zz_to_u64(&coeff);
    w += corr;
    for(i=0;i<digits && w>UINT64_MAX;i++) {
      out[i].vec->c[j] = (int16_t)(w % GADGET_BASE) - GADGET_CENTER;
      w /= GADGET_BASE;
    }
    uint64_t narrow=(uint64_t)w;
    for(;i<digits;i++) {
      out[i].vec->c[j] = (int16_t)(narrow % GADGET_BASE) - GADGET_CENTER;
      narrow /= GADGET_BASE;
    }
  }
}

static void main_witness_to_polx(polx *dst, const witness *wt) {
  size_t off = 0, t;

  for(t=0;t<DEMO_L;t++) {
    polxvec_frompolyvec(&dst[off],wt->s[b_block(t)],wt->n[b_block(t)]);
    off += wt->n[b_block(t)];
  }
  for(t=0;t<DEMO_L;t++) {
    polxvec_frompolyvec(&dst[off],wt->s[x_block(t)],wt->n[x_block(t)]);
    off += wt->n[x_block(t)];
  }
  for(t=1;t<DEMO_L;t++) {
    polxvec_frompolyvec(&dst[off],wt->s[z_block(t)],wt->n[z_block(t)]);
    off += wt->n[z_block(t)];
  }
}

static void derive_attempt_seed(uint8_t out[16], const char *label, size_t attempt) {
  uint8_t input[40];
  size_t i, labellen = strlen(label);
  const size_t header = sizeof(demo_fixed_seed);
  const size_t trailer = 2*sizeof(uint64_t);

  memcpy(input,demo_fixed_seed,sizeof(demo_fixed_seed));
  memset(input + sizeof(demo_fixed_seed),0,sizeof(input) - sizeof(demo_fixed_seed));
  if(labellen > sizeof(input) - header - trailer)
    labellen = sizeof(input) - header - trailer;
  memcpy(input + header,label,labellen);
  for(i=0;i<sizeof(uint64_t);i++)
    input[sizeof(input) - trailer + i] = (DEMO_L >> (8*i)) & 0xff;
  for(i=0;i<sizeof(uint64_t);i++)
    input[sizeof(input) - sizeof(uint64_t) + i] = (attempt >> (8*i)) & 0xff;
  shake128(out,16,input,sizeof(input));
}

static void prepare_selector_relation_seeded(prncplstmnt *st, witness *wt, const uint8_t seed_in[16]) {
  const size_t m = DEMO_M;
  const size_t gdigits = m*gadget_len_base3();
  const size_t r = 4*DEMO_L - 2;
  const size_t k = 3*DEMO_L - 1;
  const size_t deg = 1;
  const size_t ulen = (2*DEMO_L - 1)*gdigits + DEMO_L*SMALL_BLOCK;
  size_t i,j,row,t,cidx;
  size_t n[r];
  uint64_t betasq = 0;
  uint64_t nonce = 0;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  __attribute__((aligned(16)))
  uint8_t chal_seed[16];
  uint8_t *xbits = NULL;
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;
  polx *A0 = NULL, *A1 = NULL, *D = NULL;
  polx *kvec = NULL, *kprime = NULL, *rhovec = NULL;
  polx *sx = NULL, *mainx = NULL;
  polx *g0 = NULL, *g1 = NULL, *gdelta = NULL;
  poly *g1poly = NULL;
  polx *work[r];
  polx rhs[1], yx[1], c0[1], cdelta[1];
  polx ucomx[DEMO_KAPPA], minus_one[1];
  polz ucomz[DEMO_KAPPA];
  memcpy(seed,seed_in,sizeof(seed));

  for(t=0;t<DEMO_L;t++)
    n[b_block(t)] = gdigits;
  for(t=0;t<DEMO_L;t++)
    n[x_block(t)] = SMALL_BLOCK;
  for(t=1;t<DEMO_L;t++)
    n[z_block(t)] = gdigits;
  for(t=1;t<DEMO_L;t++)
    n[y_block(t)] = SMALL_BLOCK;

  init_witness_raw(wt,r,n);
  for(i=0;i<r;i++)
    polyvec_setzero(wt->s[i],wt->n[i]);
  *work = NULL;

  A0 = _aligned_alloc(64,m*gdigits*sizeof(polx));
  A1 = _aligned_alloc(64,m*gdigits*sizeof(polx));
  D  = _aligned_alloc(64,m*gdigits*sizeof(polx));
  xbits = _malloc(DEMO_L*sizeof(uint8_t));
  sx = _aligned_alloc(64,2*gdigits*sizeof(polx));
  mainx = _aligned_alloc(64,ulen*sizeof(polx));
  kvec = _aligned_alloc(64,(DEMO_L-1)*gdigits*sizeof(polx));
  kprime = _aligned_alloc(64,gdigits*sizeof(polx));
  rhovec = _aligned_alloc(64,(DEMO_L-1)*m*sizeof(polx));
  g0 = _aligned_alloc(64,gdigits*sizeof(polx));
  g1 = _aligned_alloc(64,gdigits*sizeof(polx));
  gdelta = _aligned_alloc(64,gdigits*sizeof(polx));
  g1poly = _aligned_alloc(64,gdigits*sizeof(poly));
  if(!A0 || !A1 || !D || !xbits || !sx || !mainx || !kvec || !kprime || !rhovec || !g0 || !g1 || !gdelta || !g1poly) {
    fprintf(stderr,"ERROR: allocation failed in prepare_selector_relation\n");
    exit(1);
  }

  shake128(xbits,DEMO_L,seed_in,16);

  polxvec_almostuniform(A0,m*gdigits,demo_fixed_seed,nonce);
  nonce += (m*gdigits + 31)/32;
  polxvec_almostuniform(A1,m*gdigits,demo_fixed_seed,nonce);
  for(i=0;i<m*gdigits;i++)
    polx_sub(&D[i],&A1[i],&A0[i]);

  /* public decompositions of the first columns */
  for(row=0;row<m;row++) {
    gadget_decompose_balanced3(&wt->s[b_block(DEMO_L-1)][row*(gdigits/m)],&A0[row*gdigits],gdigits/m);
    gadget_decompose_balanced3(&g1poly[row*(gdigits/m)],&A1[row*gdigits],gdigits/m);
  }
  polxvec_frompolyvec(g0,wt->s[b_block(DEMO_L-1)],gdigits);
  polxvec_frompolyvec(g1,g1poly,gdigits);
  for(i=0;i<gdigits;i++)
    polx_sub(&gdelta[i],&g1[i],&g0[i]);

  /* sample x_{L-1} and select b_{L-1} = g0 or g1 */
  poly_set_const(&wt->s[x_block(DEMO_L-1)][0],xbits[DEMO_L-1] & 1);
  if(wt->s[x_block(DEMO_L-1)][0].vec->c[0])
    polyvec_copy(wt->s[b_block(DEMO_L-1)],g1poly,gdigits);

  /* for t = L-2,...,0, define z_{t+1} = x_t * b_{t+1} and compute b_t */
  for(t=DEMO_L-1;t-- > 0;) {
    poly_set_const(&wt->s[x_block(t)][0],xbits[t] & 1);
    for(i=0;i<gdigits;i++) {
      poly_set_zero(&wt->s[z_block(t+1)][i]);
      poly_set_zero(&wt->s[b_block(t)][i]);
      if(wt->s[x_block(t)][0].vec->c[0])
        wt->s[z_block(t+1)][i] = wt->s[b_block(t+1)][i];
    }

    polxvec_frompolyvec(&sx[0],wt->s[b_block(t+1)],gdigits);
    polxvec_frompolyvec(&sx[gdigits],wt->s[z_block(t+1)],gdigits);
    for(row=0;row<m;row++) {
      polxvec_sprod(rhs,&A0[row*gdigits],&sx[0],gdigits);
      polxvec_sprod_add(rhs,&D[row*gdigits],&sx[gdigits],gdigits);
      gadget_decompose_balanced3(&wt->s[b_block(t)][row*(gdigits/m)],rhs,gdigits/m);
    }
  }

  /* u_com = Com(b_i || b_{i+1} || x_i || x_{i+1} || z_{i+1}) */
  main_witness_to_polx(mainx,wt);
  memset(ucomx,0,sizeof(ucomx));
  init_comkey(ulen*DEMO_KAPPA);
  {
    const comparams cpp = { .fu = 1, .kappa = DEMO_KAPPA };
    (void)commit_raw(ucomx,NULL,1,ulen,(const polx (*)[ulen])mainx,0,&cpp,1);
  }
  polzvec_frompolxvec(ucomz,ucomx,DEMO_KAPPA);
  {
    __attribute__((aligned(16)))
    uint8_t hashbuf[DEMO_KAPPA*N*QBYTES];
    polzvec_bitpack(hashbuf,ucomz,DEMO_KAPPA);
    shake128(chal_seed,sizeof(chal_seed),hashbuf,sizeof(hashbuf));
  }

  for(t=1;t<DEMO_L;t++)
    polxvec_challenge(&kvec[(t-1)*gdigits],gdigits,chal_seed,t-1);
  for(t=0;t<DEMO_L-1;t++)
    polxvec_challenge(&rhovec[t*m],m,chal_seed,(DEMO_L-1)+t);
  polxvec_challenge(kprime,gdigits,chal_seed,2*(DEMO_L-1));

  /* y_t = <k_t, b_t> for t=1,...,L-1 */
  for(t=1;t<DEMO_L;t++) {
    polxvec_frompolyvec(&sx[0],wt->s[b_block(t)],gdigits);
    polxvec_sprod(yx,&kvec[(t-1)*gdigits],&sx[0],gdigits);
    {
      polz yz[1];
      polz_frompolx(yz,yx);
      polz_topoly(&wt->s[y_block(t)][0],yz);
    }
  }

  /* public constants for the b_{i+1}-selection relation */
  polxvec_sprod(c0,kprime,g0,gdigits);
  polxvec_sprod(cdelta,kprime,gdelta,gdigits);

  for(i=0;i<r;i++) {
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    betasq += wt->normsq[i];
  }
  betasq += 16*N;

  init_prncplstmnt_raw(st,r,n,betasq,k,1);
  shake128_inc_init(&shakectx);
  set_const_coeff(minus_one,-1);

  cidx = 0;
  for(t=1;t<DEMO_L;t++) {
    cnst = &st->cnst[cidx++];
    buf = init_sparsecnst_half(cnst,r,2,gdigits + 1,deg,0,0);
    cnst->idx[0] = b_block(t); cnst->off[0] = 0; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
    cnst->idx[1] = y_block(t); cnst->off[1] = 0; cnst->len[1] = 1; cnst->mult[1] = 1; cnst->phi[1] = buf;
    for(j=0;j<gdigits;j++) {
      polx_neg(&cnst->phi[0][j],&kvec[(t-1)*gdigits + j]);
      absorb_polx(&shakectx,&cnst->phi[0][j]);
    }
    set_const_coeff(&cnst->phi[1][0],1);
    absorb_polx(&shakectx,&cnst->phi[1][0]);
    sparsecnst_eval(cnst->b,cnst,work,wt);
    absorb_polx(&shakectx,cnst->b);
  }

  for(t=1;t<DEMO_L;t++) {
    cnst = &st->cnst[cidx++];
    buf = init_sparsecnst_half(cnst,r,1,gdigits,deg,1,0);
    cnst->idx[0] = z_block(t); cnst->off[0] = 0; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf;
    for(j=0;j<gdigits;j++) {
      polxvec_copy(&cnst->phi[0][j],&kvec[(t-1)*gdigits + j],1);
      absorb_polx(&shakectx,&cnst->phi[0][j]);
    }
    cnst->a->len = 1;
    cnst->a->rows[0] = x_block(t-1);
    cnst->a->cols[0] = y_block(t);
    polxvec_copy(&cnst->a->coeffs[0],minus_one,1);
    absorb_polx(&shakectx,&cnst->a->coeffs[0]);
    sparsecnst_eval(cnst->b,cnst,work,wt);
    absorb_polx(&shakectx,cnst->b);
  }

  for(t=0;t<DEMO_L-1;t++) {
    cnst = &st->cnst[cidx++];
    buf = init_sparsecnst_half(cnst,r,3,3*gdigits,deg,0,0);
    cnst->idx[0] = b_block(t);   cnst->off[0] = 0; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
    cnst->idx[1] = b_block(t+1); cnst->off[1] = 0; cnst->len[1] = gdigits; cnst->mult[1] = 1; cnst->phi[1] = buf; buf += gdigits;
    cnst->idx[2] = z_block(t+1); cnst->off[2] = 0; cnst->len[2] = gdigits; cnst->mult[2] = 1; cnst->phi[2] = buf;
    for(j=0;j<gdigits;j++) {
      polxvec_setzero(&cnst->phi[0][j],1);
      polxvec_setzero(&cnst->phi[1][j],1);
      polxvec_setzero(&cnst->phi[2][j],1);
    }
    for(row=0;row<m;row++) {
      for(j=0;j<gdigits/m;j++)
        polx_scale_add(&cnst->phi[0][row*(gdigits/m) + j],&rhovec[t*m + row],pow3_u64(j));
      for(j=0;j<gdigits;j++) {
        polx tmp[1];
        polx_mul(tmp,&rhovec[t*m + row],&A0[row*gdigits + j]);
        polx_sub(&cnst->phi[1][j],&cnst->phi[1][j],tmp);
        polx_mul(tmp,&rhovec[t*m + row],&D[row*gdigits + j]);
        polx_sub(&cnst->phi[2][j],&cnst->phi[2][j],tmp);
      }
    }
    for(j=0;j<gdigits;j++) {
      absorb_polx(&shakectx,&cnst->phi[0][j]);
      absorb_polx(&shakectx,&cnst->phi[1][j]);
      absorb_polx(&shakectx,&cnst->phi[2][j]);
    }
    sparsecnst_eval(cnst->b,cnst,work,wt);
    absorb_polx(&shakectx,cnst->b);
  }

  /* <k', b_{L-1}> - cdelta*x_{L-1} - c0 = 0 */
  cnst = &st->cnst[cidx++];
  buf = init_sparsecnst_half(cnst,r,2,gdigits + 1,deg,0,0);
  cnst->idx[0] = b_block(DEMO_L-1); cnst->off[0] = 0; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
  cnst->idx[1] = x_block(DEMO_L-1); cnst->off[1] = 0; cnst->len[1] = 1; cnst->mult[1] = 1; cnst->phi[1] = buf;
  for(j=0;j<gdigits;j++) {
    polxvec_copy(&cnst->phi[0][j],&kprime[j],1);
    absorb_polx(&shakectx,&cnst->phi[0][j]);
  }
  polx_neg(&cnst->phi[1][0],cdelta);
  absorb_polx(&shakectx,&cnst->phi[1][0]);
  polxvec_copy(cnst->b,c0,1);
  absorb_polx(&shakectx,cnst->b);

  /* u_com = A*(b_i || b_{i+1} || x_i || x_{i+1} || z_{i+1}) */
  cnst = &st->cnst[cidx++];
  buf = init_sparsecnst_half(cnst,r,3*DEMO_L-1,ulen,DEMO_KAPPA,0,0);
  i = 0;
  for(t=0;t<DEMO_L;t++) {
    cnst->idx[t] = b_block(t);
    cnst->off[t] = 0;
    cnst->len[t] = gdigits;
    cnst->mult[t] = 1;
    cnst->phi[t] = buf;
    for(j=0;j<gdigits;j++) {
      polxvec_copy(&cnst->phi[t][j],&comkey[i++],1);
      absorb_polx(&shakectx,&cnst->phi[t][j]);
    }
    buf += gdigits;
  }
  for(t=0;t<DEMO_L;t++) {
    const size_t idx = DEMO_L + t;
    cnst->idx[idx] = x_block(t);
    cnst->off[idx] = 0;
    cnst->len[idx] = SMALL_BLOCK;
    cnst->mult[idx] = 1;
    cnst->phi[idx] = buf;
    for(j=0;j<SMALL_BLOCK;j++) {
      polxvec_copy(&cnst->phi[idx][j],&comkey[i++],1);
      absorb_polx(&shakectx,&cnst->phi[idx][j]);
    }
    buf += SMALL_BLOCK;
  }
  for(t=1;t<DEMO_L;t++) {
    const size_t idx = 2*DEMO_L + (t-1);
    cnst->idx[idx] = z_block(t);
    cnst->off[idx] = 0;
    cnst->len[idx] = gdigits;
    cnst->mult[idx] = 1;
    cnst->phi[idx] = buf;
    for(j=0;j<gdigits;j++) {
      polxvec_copy(&cnst->phi[idx][j],&comkey[i++],1);
      absorb_polx(&shakectx,&cnst->phi[idx][j]);
    }
    buf += gdigits;
  }
  polxvec_copy(cnst->b,ucomx,DEMO_KAPPA);
  for(j=0;j<DEMO_KAPPA;j++)
    absorb_polx(&shakectx,&cnst->b[j]);

  if(cidx != k) {
    fprintf(stderr,"ERROR: internal constraint count mismatch (%zu != %zu)\n",cidx,k);
    exit(1);
  }

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);

  free(A0);
  free(A1);
  free(D);
  free(xbits);
  free(kvec);
  free(kprime);
  free(rhovec);
  free(g0);
  free(g1);
  free(gdelta);
  free(g1poly);
  free(*work);
  free(sx);
  free(mainx);
}

static void prepare_selector_relation(prncplstmnt *st, witness *wt) {
  prepare_selector_relation_seeded(st,wt,demo_fixed_seed);
}

static int test_twolayer(void) {
  int ret = 1;
  prncplstmnt st0 = {};
  statement st1 = {};
  proof pi0 = {};
  witness wt0 = {}, wt1 = {};
  double size = 0;
  size_t attempt;
  uint8_t seed[16];

  printf("Testing Chihuahua hashed selector %s frontend\n\n",GADGET_LABEL);
  for(attempt=0;attempt<MAX_ATTEMPTS;attempt++) {
    derive_attempt_seed(seed,"two-layer",attempt);
    prepare_selector_relation_seeded(&st0,&wt0,seed);
    ret = principle_verify(&st0,&wt0);
    if(!ret)
      break;
    free_prncplstmnt(&st0);
    free_witness(&wt0);
  }
  if(ret) {
    fprintf(stderr,"ERROR: Verification of selector %s relation instance failed after %d attempts: %d\n",GADGET_LABEL,MAX_ATTEMPTS,ret);
    goto end;
  }
  printf("[two-layer] attempt %zu succeeded\n\n",attempt);

  print_prncplstmnt_pp(&st0);
  ret = principle_prove(&st1,&wt1,&pi0,&st0,&wt0,0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt0);
  size += print_proof_pp(&pi0);
  print_statement_pp(&st1);
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st1);
  ret = principle_reduce(&st1,&pi0,&st0);
  free_prncplstmnt(&st0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of Chihuahua reduction failed: %d\n",ret);
    goto end;
  }

  size += print_witness_pp(&wt1);
  printf("Reduced proof size: %.2f KB\n\n",size);

end:
  free_prncplstmnt(&st0);
  free_statement(&st1);
  free_proof(&pi0);
  free_witness(&wt0);
  free_witness(&wt1);
  return ret;
}

static int test_pack(void) {
  int ret = 1;
  prncplstmnt st = {};
  witness wt = {};
  composite p = {};
  size_t attempt;
  uint8_t seed[16];

  printf("Testing Chihuahua hashed selector %s composite frontend\n\n",GADGET_LABEL);
  for(attempt=0;attempt<MAX_ATTEMPTS;attempt++) {
    derive_attempt_seed(seed,"pack",attempt);
    prepare_selector_relation_seeded(&st,&wt,seed);
    ret = principle_verify(&st,&wt);
    if(ret) {
      free_prncplstmnt(&st);
      free_witness(&wt);
      continue;
    }
    print_prncplstmnt_pp(&st);

    ret = composite_prove_principle(&p,&st,&wt);
    if(!ret)
      ret = composite_verify_principle(&p,&st);
    if(!ret)
      break;

    free_prncplstmnt(&st);
    free_composite(&p);
    free_witness(&wt);
  }
  if(ret) {
    fprintf(stderr,"ERROR: Hashed selector composite proof/verification failed after %d attempts: %d\n",MAX_ATTEMPTS,ret);
    goto end;
  }
  printf("[composite] attempt %zu succeeded\n\n",attempt);

end:
  free_prncplstmnt(&st);
  free_composite(&p);
  free_witness(&wt);
  return ret;
}

#ifndef SELECTOR_NO_MAIN
int main(void) {
  int ret;

  ret = test_twolayer();
  if(ret) goto end;
  ret = test_pack();
  free_comkey();
end:
  return ret;
}
#endif
