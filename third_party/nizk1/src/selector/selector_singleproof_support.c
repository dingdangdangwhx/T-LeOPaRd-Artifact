#ifndef DEMO_L
#define DEMO_L 4
#endif
#ifndef SINGLEPROOF_MAX_ATTEMPTS
#define SINGLEPROOF_MAX_ATTEMPTS 64
#endif
#ifndef SELECTOR_MAX_ATTEMPTS
#define SELECTOR_MAX_ATTEMPTS SINGLEPROOF_MAX_ATTEMPTS
#endif
#ifndef SINGLEPROOF_ENABLE_FALLBACK
#define SINGLEPROOF_ENABLE_FALLBACK 0
#endif
#ifndef SINGLEPROOF_SEARCH_ALL_ATTEMPTS
#define SINGLEPROOF_SEARCH_ALL_ATTEMPTS 0
#endif
#ifndef SINGLEPROOF_GOOD_ENOUGH_KB
#define SINGLEPROOF_GOOD_ENOUGH_KB 60.0
#endif
#ifndef SINGLEPROOF_EXTRA_ATTEMPTS_AFTER_SUCCESS
#define SINGLEPROOF_EXTRA_ATTEMPTS_AFTER_SUCCESS 8
#endif
#ifndef SINGLEPROOF_PROGRESS_LOG
#define SINGLEPROOF_PROGRESS_LOG 0
#endif
#ifndef SINGLEPROOF_FIXED_PUBLIC_MATRICES
#define SINGLEPROOF_FIXED_PUBLIC_MATRICES 1
#endif
#ifndef SINGLEPROOF_STOP_ON_FALLBACK
#define SINGLEPROOF_STOP_ON_FALLBACK 0
#endif
#ifndef SINGLEPROOF_COMPOSITE_TIMEOUT_SEC
#define SINGLEPROOF_COMPOSITE_TIMEOUT_SEC 0
#endif
#ifndef Y_SPLIT_GAMMA
#define Y_SPLIT_GAMMA 200
#endif
#include <signal.h>
#include <setjmp.h>
#include <unistd.h>
#define SELECTOR_NO_MAIN
#include "selector_base.c"

enum {
  B_BLOCK = 0,
  X_BLOCK = 1,
  Z_BLOCK = 2,
  Y0_BLOCK = 3,
  Y1_BLOCK = 4,
  INNER_COMPOSITE_ATTEMPTS = 1
};

#if SINGLEPROOF_COMPOSITE_TIMEOUT_SEC > 0
static jmp_buf singleproof_timeout_jmp;
static volatile sig_atomic_t singleproof_timeout_active = 0;

static void singleproof_timeout_handler(int sig) {
  (void)sig;
  if(singleproof_timeout_active)
    longjmp(singleproof_timeout_jmp,1);
}

static int singleproof_composite_call_with_timeout(
    composite *p, const prncplstmnt *st, const witness *wt,
    int (*fn)(composite *, const prncplstmnt *, const witness *)) {
  void (*old_handler)(int);
  int ret;

  old_handler = signal(SIGALRM,singleproof_timeout_handler);
  singleproof_timeout_active = 1;
  if(setjmp(singleproof_timeout_jmp)) {
    alarm(0);
    singleproof_timeout_active = 0;
    signal(SIGALRM,old_handler);
    free_composite(p);
    memset(p,0,sizeof(*p));
    return 90;
  }

  alarm(SINGLEPROOF_COMPOSITE_TIMEOUT_SEC);
  ret = fn(p,st,wt);
  alarm(0);
  singleproof_timeout_active = 0;
  signal(SIGALRM,old_handler);
  return ret;
}
#else
static int singleproof_composite_call_with_timeout(
    composite *p, const prncplstmnt *st, const witness *wt,
    int (*fn)(composite *, const prncplstmnt *, const witness *)) {
  return fn(p,st,wt);
}
#endif

static void singleproof_main_witness_to_polx(polx *dst, const witness *wt) {
  size_t off = 0;

  polxvec_frompolyvec(&dst[off],wt->s[B_BLOCK],wt->n[B_BLOCK]);
  off += wt->n[B_BLOCK];
  polxvec_frompolyvec(&dst[off],wt->s[X_BLOCK],wt->n[X_BLOCK]);
  off += wt->n[X_BLOCK];
  polxvec_frompolyvec(&dst[off],wt->s[Z_BLOCK],wt->n[Z_BLOCK]);
}

static void prepare_selector_singleproof_seeded(prncplstmnt *st, witness *wt, const uint8_t seed_in[16]) {
  const size_t m = DEMO_M;
  const size_t gdigits = m*gadget_len_base3();
  const size_t Lchain = DEMO_L;
  const size_t r = 5;
  const size_t k = 2*Lchain + 1;
  const size_t deg = 1;
  const size_t b_len = Lchain*gdigits;
  const size_t x_len = Lchain*SMALL_BLOCK;
  const size_t z_len = (Lchain-1)*gdigits;
  const size_t y_len = Lchain*SMALL_BLOCK;
  const size_t ulen = b_len + x_len + z_len;
  size_t i,j,row,t,cidx;
  size_t n[r];
  uint64_t betasq = 0;
  uint64_t nonce = 0;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  __attribute__((aligned(16)))
  uint8_t chal_seed[16];
  __attribute__((aligned(16)))
  uint8_t matrix_seed[16];
  uint8_t *xbits = NULL;
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;
  polx *A0 = NULL, *A1 = NULL, *D = NULL;
  polx *kvec = NULL, *kprime = NULL, *rhovec = NULL, *alphavec = NULL;
  polx *sx = NULL, *mainx = NULL;
  polx *g0 = NULL, *g1 = NULL, *gdelta = NULL;
  poly *g1poly = NULL;
  polx *work[r];
  polx rhs[1], yx[1], c0[1], cdelta[1], scaled[1];
  polx ucomx[DEMO_KAPPA], minus_one[1], gamma_coeff[1], minus_gamma[1];
  polz ucomz[DEMO_KAPPA];

  memcpy(seed,seed_in,sizeof(seed));
  if(SINGLEPROOF_FIXED_PUBLIC_MATRICES)
    memcpy(matrix_seed,demo_fixed_seed,sizeof(matrix_seed));
  else
    memcpy(matrix_seed,seed_in,sizeof(matrix_seed));

  n[B_BLOCK] = b_len;
  n[X_BLOCK] = x_len;
  n[Z_BLOCK] = z_len;
  n[Y0_BLOCK] = y_len;
  n[Y1_BLOCK] = y_len;

  init_witness_raw(wt,r,n);
  for(i=0;i<r;i++)
    polyvec_setzero(wt->s[i],wt->n[i]);
  *work = NULL;

  A0 = _aligned_alloc(64,m*gdigits*sizeof(polx));
  A1 = _aligned_alloc(64,m*gdigits*sizeof(polx));
  D  = _aligned_alloc(64,m*gdigits*sizeof(polx));
  xbits = _malloc(Lchain*sizeof(uint8_t));
  sx = _aligned_alloc(64,2*gdigits*sizeof(polx));
  mainx = _aligned_alloc(64,ulen*sizeof(polx));
  kvec = _aligned_alloc(64,(Lchain-1)*gdigits*sizeof(polx));
  kprime = _aligned_alloc(64,gdigits*sizeof(polx));
  rhovec = _aligned_alloc(64,(Lchain-1)*m*sizeof(polx));
  alphavec = _aligned_alloc(64,(Lchain-1)*sizeof(polx));
  g0 = _aligned_alloc(64,gdigits*sizeof(polx));
  g1 = _aligned_alloc(64,gdigits*sizeof(polx));
  gdelta = _aligned_alloc(64,gdigits*sizeof(polx));
  g1poly = _aligned_alloc(64,gdigits*sizeof(poly));
  if(!A0 || !A1 || !D || !xbits || !sx || !mainx || !kvec || !kprime || !rhovec || !alphavec || !g0 || !g1 || !gdelta || !g1poly) {
    fprintf(stderr,"ERROR: allocation failed in prepare_selector_singleproof\n");
    exit(1);
  }

  shake128(xbits,Lchain,seed_in,16);

  polxvec_almostuniform(A0,m*gdigits,matrix_seed,nonce);
  nonce += (m*gdigits + 31)/32;
  polxvec_almostuniform(A1,m*gdigits,matrix_seed,nonce);
  for(i=0;i<m*gdigits;i++)
    polx_sub(&D[i],&A1[i],&A0[i]);

  /* public decompositions of the first columns */
  for(row=0;row<m;row++) {
    gadget_decompose_balanced3(&wt->s[B_BLOCK][(Lchain-1)*gdigits + row*(gdigits/m)],&A0[row*gdigits],gdigits/m);
    gadget_decompose_balanced3(&g1poly[row*(gdigits/m)],&A1[row*gdigits],gdigits/m);
  }
  polxvec_frompolyvec(g0,&wt->s[B_BLOCK][(Lchain-1)*gdigits],gdigits);
  polxvec_frompolyvec(g1,g1poly,gdigits);
  for(i=0;i<gdigits;i++)
    polx_sub(&gdelta[i],&g1[i],&g0[i]);

  poly_set_const(&wt->s[X_BLOCK][(Lchain-1)*SMALL_BLOCK],xbits[Lchain-1] & 1);
  if(wt->s[X_BLOCK][(Lchain-1)*SMALL_BLOCK].vec->c[0])
    polyvec_copy(&wt->s[B_BLOCK][(Lchain-1)*gdigits],g1poly,gdigits);

  for(t=Lchain-1;t-- > 0;) {
    poly_set_const(&wt->s[X_BLOCK][t*SMALL_BLOCK],xbits[t] & 1);
    if(wt->s[X_BLOCK][t*SMALL_BLOCK].vec->c[0])
      polyvec_copy(&wt->s[Z_BLOCK][t*gdigits],&wt->s[B_BLOCK][(t+1)*gdigits],gdigits);
    else
      polyvec_setzero(&wt->s[Z_BLOCK][t*gdigits],gdigits);

    polxvec_frompolyvec(&sx[0],&wt->s[B_BLOCK][(t+1)*gdigits],gdigits);
    polxvec_frompolyvec(&sx[gdigits],&wt->s[Z_BLOCK][t*gdigits],gdigits);
    for(row=0;row<m;row++) {
      polxvec_sprod(rhs,&A0[row*gdigits],&sx[0],gdigits);
      polxvec_sprod_add(rhs,&D[row*gdigits],&sx[gdigits],gdigits);
      gadget_decompose_balanced3(&wt->s[B_BLOCK][t*gdigits + row*(gdigits/m)],rhs,gdigits/m);
    }
  }

  singleproof_main_witness_to_polx(mainx,wt);
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

  for(t=0;t<Lchain-1;t++) {
    polxvec_challenge(&kvec[t*gdigits],gdigits,chal_seed,(uint64_t)t);
    polxvec_challenge(&rhovec[t*m],m,chal_seed,32 + (uint64_t)t);
  }
  polxvec_challenge(kprime,gdigits,chal_seed,64);
  polxvec_challenge(alphavec,Lchain-1,chal_seed,65);

  /* Y0 + Gamma*Y1 stores alpha_t * <k_t, b_{t+1}> in segment t */
  for(t=0;t<Lchain-1;t++) {
    poly tmpy[1];
    polz yz[1];

    polyvec_setzero(tmpy,1);
    polxvec_frompolyvec(&sx[0],&wt->s[B_BLOCK][(t+1)*gdigits],gdigits);
    polxvec_sprod(yx,&kvec[t*gdigits],&sx[0],gdigits);
    polx_mul(scaled,&alphavec[t],yx);
    polz_frompolx(yz,scaled);
    polz_topoly(tmpy,yz);
    for(j=0;j<N;j++) {
      int64_t c = tmpy->vec->c[j];
      int64_t q;
      int64_t rem;

      if(c >= 0)
        q = (c + Y_SPLIT_GAMMA/2) / Y_SPLIT_GAMMA;
      else
        q = -(((-c) + Y_SPLIT_GAMMA/2) / Y_SPLIT_GAMMA);
      rem = c - Y_SPLIT_GAMMA*q;
      wt->s[Y0_BLOCK][t*SMALL_BLOCK].vec->c[j] = (int16_t)rem;
      wt->s[Y1_BLOCK][t*SMALL_BLOCK].vec->c[j] = (int16_t)q;
    }
  }

  polxvec_frompolyvec(&sx[0],&wt->s[B_BLOCK][(Lchain-1)*gdigits],gdigits);
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
  set_const_coeff(gamma_coeff,Y_SPLIT_GAMMA);
  set_const_coeff(minus_gamma,-Y_SPLIT_GAMMA);
  cidx = 0;

  /* top selector: <k', b_{L-1}> - cdelta*x_{L-1} - c0 = 0 */
  cnst = &st->cnst[cidx++];
  buf = init_sparsecnst_half(cnst,r,2,gdigits + 1,deg,0,0);
  cnst->idx[0] = B_BLOCK; cnst->off[0] = (Lchain-1)*gdigits; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
  cnst->idx[1] = X_BLOCK; cnst->off[1] = (Lchain-1)*SMALL_BLOCK; cnst->len[1] = 1; cnst->mult[1] = 1; cnst->phi[1] = buf;
  for(j=0;j<gdigits;j++) {
    polxvec_copy(&cnst->phi[0][j],&kprime[j],1);
    absorb_polx(&shakectx,&cnst->phi[0][j]);
  }
  polx_neg(&cnst->phi[1][0],cdelta);
  absorb_polx(&shakectx,&cnst->phi[1][0]);
  polxvec_copy(cnst->b,c0,1);
  absorb_polx(&shakectx,cnst->b);

  /* for t=0..L-2: Y0_t + Gamma*Y1_t - alpha_t <k_t, b_{t+1}> = 0 */
  for(t=0;t<Lchain-1;t++) {
    cnst = &st->cnst[cidx++];
    buf = init_sparsecnst_half(cnst,r,3,gdigits + 2,deg,0,0);
    cnst->idx[0] = B_BLOCK; cnst->off[0] = (t+1)*gdigits; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
    cnst->idx[1] = Y0_BLOCK; cnst->off[1] = t*SMALL_BLOCK; cnst->len[1] = 1; cnst->mult[1] = 1; cnst->phi[1] = buf; buf += 1;
    cnst->idx[2] = Y1_BLOCK; cnst->off[2] = t*SMALL_BLOCK; cnst->len[2] = 1; cnst->mult[2] = 1; cnst->phi[2] = buf;
    for(j=0;j<gdigits;j++) {
      polx_mul(&cnst->phi[0][j],&alphavec[t],&kvec[t*gdigits + j]);
      polx_neg(&cnst->phi[0][j],&cnst->phi[0][j]);
      absorb_polx(&shakectx,&cnst->phi[0][j]);
    }
    set_const_coeff(&cnst->phi[1][0],1);
    absorb_polx(&shakectx,&cnst->phi[1][0]);
    polxvec_copy(&cnst->phi[2][0],gamma_coeff,1);
    absorb_polx(&shakectx,&cnst->phi[2][0]);
    sparsecnst_eval(cnst->b,cnst,work,wt);
    absorb_polx(&shakectx,cnst->b);
  }

  /* packed multiplication: sum alpha_t <k_t,z_{t+1}> - <X,Y0> - Gamma*<X,Y1> = 0 */
  cnst = &st->cnst[cidx++];
  buf = init_sparsecnst_half(cnst,r,1,z_len,deg,1,0);
  cnst->idx[0] = Z_BLOCK; cnst->off[0] = 0; cnst->len[0] = z_len; cnst->mult[0] = 1; cnst->phi[0] = buf;
  for(t=0;t<Lchain-1;t++) {
    for(j=0;j<gdigits;j++) {
      polx_mul(&cnst->phi[0][t*gdigits + j],&alphavec[t],&kvec[t*gdigits + j]);
      absorb_polx(&shakectx,&cnst->phi[0][t*gdigits + j]);
    }
  }
  cnst->a->len = 2;
  cnst->a->rows[0] = X_BLOCK;
  cnst->a->cols[0] = Y0_BLOCK;
  polxvec_copy(&cnst->a->coeffs[0],minus_one,1);
  absorb_polx(&shakectx,&cnst->a->coeffs[0]);
  cnst->a->rows[1] = X_BLOCK;
  cnst->a->cols[1] = Y1_BLOCK;
  polxvec_copy(&cnst->a->coeffs[1],minus_gamma,1);
  absorb_polx(&shakectx,&cnst->a->coeffs[1]);
  sparsecnst_eval(cnst->b,cnst,work,wt);
  absorb_polx(&shakectx,cnst->b);

  /* for t=0..L-2: <rho_t, G3 b_t - A0 b_{t+1} - D z_{t+1}> = 0 */
  for(t=0;t<Lchain-1;t++) {
    cnst = &st->cnst[cidx++];
    buf = init_sparsecnst_half(cnst,r,3,3*gdigits,deg,0,0);
    cnst->idx[0] = B_BLOCK; cnst->off[0] = t*gdigits; cnst->len[0] = gdigits; cnst->mult[0] = 1; cnst->phi[0] = buf; buf += gdigits;
    cnst->idx[1] = B_BLOCK; cnst->off[1] = (t+1)*gdigits; cnst->len[1] = gdigits; cnst->mult[1] = 1; cnst->phi[1] = buf; buf += gdigits;
    cnst->idx[2] = Z_BLOCK; cnst->off[2] = t*gdigits; cnst->len[2] = gdigits; cnst->mult[2] = 1; cnst->phi[2] = buf;
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

  /* u_com well-formedness for B||X||Z */
  cnst = &st->cnst[cidx++];
  buf = init_sparsecnst_half(cnst,r,3,ulen,DEMO_KAPPA,0,0);
  cnst->idx[0] = B_BLOCK; cnst->off[0] = 0; cnst->len[0] = b_len; cnst->mult[0] = 1; cnst->phi[0] = buf;
  for(j=0;j<b_len;j++) {
    polxvec_copy(&cnst->phi[0][j],&comkey[j],1);
    absorb_polx(&shakectx,&cnst->phi[0][j]);
  }
  buf += b_len;
  cnst->idx[1] = X_BLOCK; cnst->off[1] = 0; cnst->len[1] = x_len; cnst->mult[1] = 1; cnst->phi[1] = buf;
  for(j=0;j<x_len;j++) {
    polxvec_copy(&cnst->phi[1][j],&comkey[b_len + j],1);
    absorb_polx(&shakectx,&cnst->phi[1][j]);
  }
  buf += x_len;
  cnst->idx[2] = Z_BLOCK; cnst->off[2] = 0; cnst->len[2] = z_len; cnst->mult[2] = 1; cnst->phi[2] = buf;
  for(j=0;j<z_len;j++) {
    polxvec_copy(&cnst->phi[2][j],&comkey[b_len + x_len + j],1);
    absorb_polx(&shakectx,&cnst->phi[2][j]);
  }
  polxvec_copy(cnst->b,ucomx,DEMO_KAPPA);
  for(j=0;j<DEMO_KAPPA;j++)
    absorb_polx(&shakectx,&cnst->b[j]);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);

  free(A0);
  free(A1);
  free(D);
  free(xbits);
  free(sx);
  free(mainx);
  free(kvec);
  free(kprime);
  free(rhovec);
  free(alphavec);
  free(g0);
  free(g1);
  free(gdelta);
  free(g1poly);
  free(*work);
}

static void prepare_selector_singleproof(prncplstmnt *st, witness *wt) {
  prepare_selector_singleproof_seeded(st,wt,demo_fixed_seed);
}

static int test_singleproof_twolayer(void) {
  int ret = 1;
  size_t attempt;
  uint8_t seed[16];

  printf("Testing Chihuahua hashed selector %s singleproof frontend\n\n",GADGET_LABEL);
  for(attempt=0;attempt<MAX_ATTEMPTS;attempt++) {
    prncplstmnt st0 = {};
    statement st1 = {};
    proof pi0 = {};
    witness wt0 = {}, wt1 = {};
    double size = 0;

    derive_attempt_seed(seed,"singleproof-two-layer",attempt);
    prepare_selector_singleproof_seeded(&st0,&wt0,seed);
    ret = principle_verify(&st0,&wt0);
    if(ret) {
      free_prncplstmnt(&st0);
      free_witness(&wt0);
      continue;
    }

    ret = principle_prove(&st1,&wt1,&pi0,&st0,&wt0,0);
    if(ret) {
      free_prncplstmnt(&st0);
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt0);
      free_witness(&wt1);
      continue;
    }
    ret = verify(&st1,&wt1);
    if(ret) {
      free_prncplstmnt(&st0);
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt0);
      free_witness(&wt1);
      continue;
    }

    printf("[singleproof two-layer] attempt %zu succeeded\n\n",attempt);
    print_prncplstmnt_pp(&st0);
    free_witness(&wt0);
    size += print_proof_pp(&pi0);
    print_statement_pp(&st1);

    free_statement(&st1);
    ret = principle_reduce(&st1,&pi0,&st0);
    free_prncplstmnt(&st0);
    if(ret) {
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt1);
      continue;
    }
    ret = verify(&st1,&wt1);
    if(ret) {
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt1);
      continue;
    }

    size += print_witness_pp(&wt1);
    printf("Reduced singleproof size: %.2f KB\n\n",size);
    free_statement(&st1);
    free_proof(&pi0);
    free_witness(&wt1);
    break;
  }
  if(ret) {
    fprintf(stderr,"ERROR: Verification of selector %s singleproof instance failed after %d attempts: %d\n",GADGET_LABEL,MAX_ATTEMPTS,ret);
    goto end;
  }
  ret = 0;

end:
  return ret;
}

static int composite_prove_principle_singleproof_stable(composite *p, const prncplstmnt *st, const witness *wt) {
  int ret;
  size_t i = 0;
  statement tst[2] = {};
  witness twt[2] = {};
  statement chk = {};
  double twtsize[2], pisize;

  p->l = 0;
  p->size = 0.0;
  memset(&p->owt,0,sizeof(witness));
  memset(p->pi,0,sizeof(p->pi));

  p->pi[p->l] = _malloc(sizeof(proof));
  ret = principle_prove(tst,twt,p->pi[p->l],st,wt,0);
  if(ret)
    goto err;
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l += 1;

  while(p->l < 16) {
    p->pi[p->l] = _malloc(sizeof(proof));
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],0);
    if(ret) {
      free(p->pi[p->l]);
      p->pi[p->l] = NULL;
      break;
    }

    ret = reduce(&chk,p->pi[p->l],&tst[i]);
    if(ret || verify(&chk,&twt[i^1])) {
      free_statement(&chk);
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
      p->pi[p->l] = NULL;
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }
    free_statement(&chk);

    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
      p->pi[p->l] = NULL;
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }

    free_statement(&tst[i]);
    free_witness(&twt[i]);
    p->size += pisize;
    p->l += 1;
    i ^= 1;
  }

  if(p->l < 16) {
    p->pi[p->l] = _malloc(sizeof(proof));
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],1);
    if(!ret) {
      ret = reduce(&chk,p->pi[p->l],&tst[i]);
      if(!ret)
        ret = verify(&chk,&twt[i^1]);
      free_statement(&chk);
    }
    if(ret) {
      free(p->pi[p->l]);
      p->pi[p->l] = NULL;
    }
    else {
      pisize = print_proof_pp(p->pi[p->l]);
      print_statement_pp(&tst[i^1]);
      twtsize[i^1] = print_witness_pp(&twt[i^1]);
      if(pisize + twtsize[i^1] >= twtsize[i]) {
        free_proof(p->pi[p->l]);
        free(p->pi[p->l]);
        p->pi[p->l] = NULL;
      }
      else {
        p->size += pisize;
        p->l += 1;
        free_statement(&tst[i]);
        free_witness(&twt[i]);
        i ^= 1;
      }
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
    }
  }

  free_statement(&tst[i]);
  p->owt = twt[i];
  p->size += twtsize[i];
  return 0;

err:
  free_statement(&chk);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free_composite(p);
  return ret;
}

static int composite_verify_principle_prefix(size_t l, proof *const pi[], const witness *owt,
                                             const prncplstmnt *st) {
  int ret;
  size_t i = 0, j;
  statement tst[2] = {};

  if(!l)
    return 1;

  ret = principle_reduce(tst,pi[0],st);
  if(ret)
    goto err;

  for(j=1;j<l;j++) {
    ret = reduce(&tst[i^1],pi[j],&tst[i]);
    free_statement(&tst[i]);
    if(ret)
      goto err;
    i ^= 1;
  }

  ret = verify(&tst[i],owt);

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

static int composite_prove_principle_singleproof_backoff(composite *p, const prncplstmnt *st, const witness *wt) {
  int ret = 1;
  size_t i, l, count = 0;
  statement tst[16] = {};
  witness twt[16] = {};
  proof *pi[16] = {};
  double pisize[16] = {0}, twtsize[16] = {0};
  double prefsize = 0.0, bestsize = 1e300;
  size_t bestl = 0;
  clock_t t;

  memset(&p->owt,0,sizeof(witness));
  p->l = 0;
  p->size = 0.0;

  t = clock();
  pi[0] = _malloc(sizeof(proof));
  ret = principle_prove(tst,twt,pi[0],st,wt,0);
  if(ret)
    goto err;
  pisize[0] = print_proof_pp(pi[0]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  count = 1;

  while(count < 16) {
    pi[count] = _malloc(sizeof(proof));
    ret = prove(&tst[count],&twt[count],pi[count],&tst[count-1],&twt[count-1],0);
    if(ret) {
      free(pi[count]);
      pi[count] = NULL;
      ret = 0;
      break;
    }

    pisize[count] = print_proof_pp(pi[count]);
    print_statement_pp(&tst[count]);
    twtsize[count] = print_witness_pp(&twt[count]);
    if(pisize[count] + twtsize[count] >= twtsize[count-1]) {
      free_proof(pi[count]);
      free(pi[count]);
      pi[count] = NULL;
      free_statement(&tst[count]);
      free_witness(&twt[count]);
      break;
    }
    count += 1;
  }

  if(count < 16) {
    pi[count] = _malloc(sizeof(proof));
    ret = prove(&tst[count],&twt[count],pi[count],&tst[count-1],&twt[count-1],1);
    if(ret) {
      free(pi[count]);
      pi[count] = NULL;
      ret = 0;
    }
    else {
      pisize[count] = print_proof_pp(pi[count]);
      print_statement_pp(&tst[count]);
      twtsize[count] = print_witness_pp(&twt[count]);
      if(pisize[count] + twtsize[count] >= twtsize[count-1]) {
        free_proof(pi[count]);
        free(pi[count]);
        pi[count] = NULL;
        free_statement(&tst[count]);
        free_witness(&twt[count]);
      }
      else
        count += 1;
    }
  }

  for(l=1;l<=count;l++) {
    prefsize += pisize[l-1];
    if(prefsize + twtsize[l-1] >= bestsize)
      continue;
    if(!composite_verify_principle_prefix(l,pi,&twt[l-1],st)) {
      bestsize = prefsize + twtsize[l-1];
      bestl = l;
    }
  }

  if(!bestl) {
    ret = 1;
    goto err;
  }

  p->l = bestl;
  p->size = bestsize;
  for(i=0;i<bestl;i++) {
    p->pi[i] = pi[i];
    pi[i] = NULL;
  }
  p->owt = twt[bestl-1];
  memset(&twt[bestl-1],0,sizeof(witness));
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Chihuahua Pack members: %zu\n",p->l);
  printf("Chihuahua Pack size: %.2f KB\n",p->size);
  printf("Chihuahua Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);

  ret = 0;

err:
  for(i=0;i<16;i++) {
    free_statement(&tst[i]);
    free_witness(&twt[i]);
    if(pi[i]) {
      free_proof(pi[i]);
      free(pi[i]);
    }
  }
  if(ret)
    free_composite(p);
  return ret;
}

static int test_singleproof_pack(void) {
  int ret = 1;
  prncplstmnt st = {};
  witness wt = {};
  composite p = {};
  composite best = {};
  size_t best_attempt = 0;
  int have_best = 0;
  size_t stop_after_attempt = (size_t)-1;
#if SINGLEPROOF_ENABLE_FALLBACK
  composite fallback = {};
  size_t fallback_attempt = 0;
  int have_fallback = 0;
#endif
  size_t attempt;
  uint8_t seed[16];

  printf("Testing Chihuahua hashed selector %s singleproof composite frontend\n\n",GADGET_LABEL);
  for(attempt=0;attempt<MAX_ATTEMPTS;attempt++) {
    statement st1 = {};
    proof pi0 = {};
    witness wt1 = {};

    derive_attempt_seed(seed,"singleproof-pack",attempt);
    if(SINGLEPROOF_PROGRESS_LOG)
      fprintf(stderr,"[singleproof composite] attempt %zu/%zu\n",attempt,(size_t)(MAX_ATTEMPTS-1));
    prepare_selector_singleproof_seeded(&st,&wt,seed);
    ret = principle_verify(&st,&wt);
    if(ret) {
      free_prncplstmnt(&st);
      free_witness(&wt);
      continue;
    }

    ret = principle_prove(&st1,&wt1,&pi0,&st,&wt,0);
    if(ret) {
      free_prncplstmnt(&st);
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt);
      free_witness(&wt1);
      continue;
    }
    ret = verify(&st1,&wt1);
    if(ret) {
      free_prncplstmnt(&st);
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt);
      free_witness(&wt1);
      continue;
    }
    free_statement(&st1);
    ret = principle_reduce(&st1,&pi0,&st);
    if(ret) {
      free_prncplstmnt(&st);
      free_statement(&st1);
      free_proof(&pi0);
      free_witness(&wt);
      free_witness(&wt1);
      continue;
    }
    ret = verify(&st1,&wt1);
    free_statement(&st1);
    free_proof(&pi0);
    free_witness(&wt1);
    if(ret) {
      free_prncplstmnt(&st);
      free_witness(&wt);
      continue;
    }

    print_prncplstmnt_pp(&st);
    {
      size_t inner;
      int found = 0;
      int last_ret = ret;
      composite attempt_best = {};

      for(inner=0;inner<INNER_COMPOSITE_ATTEMPTS;inner++) {
        composite cand = {};
        ret = singleproof_composite_call_with_timeout(&cand,&st,&wt,
            composite_prove_principle);
        if(!ret)
          ret = composite_verify_principle(&cand,&st);
        if(!ret) {
          if(!found || cand.size < attempt_best.size) {
            if(found)
              free_composite(&attempt_best);
            attempt_best = cand;
            found = 1;
            print_prncplstmnt_pp(&st);
          }
          else
            free_composite(&cand);
        }
        else {
          last_ret = ret;
          free_composite(&cand);
        }
      }

      if(found) {
        ret = 0;
        if(!have_best || attempt_best.size < best.size) {
          if(have_best)
            free_composite(&best);
          best = attempt_best;
          memset(&attempt_best,0,sizeof(attempt_best));
          best_attempt = attempt;
          have_best = 1;
          if(best.size <= SINGLEPROOF_GOOD_ENOUGH_KB)
            stop_after_attempt = attempt;
          else if(stop_after_attempt == (size_t)-1)
            stop_after_attempt = attempt + SINGLEPROOF_EXTRA_ATTEMPTS_AFTER_SUCCESS;
        }
        else
          free_composite(&attempt_best);
      }
      else {
        ret = last_ret;
#if SINGLEPROOF_ENABLE_FALLBACK
        {
          composite cand = {};
          int fb_ret;

          fb_ret = singleproof_composite_call_with_timeout(&cand,&st,&wt,
              composite_prove_principle_singleproof_backoff);
          if(!fb_ret)
            fb_ret = composite_verify_principle(&cand,&st);
          if(fb_ret) {
            free_composite(&cand);
            memset(&cand,0,sizeof(cand));
            fb_ret = singleproof_composite_call_with_timeout(&cand,&st,&wt,
                composite_prove_principle_singleproof_stable);
            if(!fb_ret)
              fb_ret = composite_verify_principle(&cand,&st);
          }

          if(!fb_ret) {
            if(!have_fallback || cand.size < fallback.size) {
              if(have_fallback)
                free_composite(&fallback);
              fallback = cand;
              fallback_attempt = attempt;
              have_fallback = 1;
            }
            else
              free_composite(&cand);
#if SINGLEPROOF_STOP_ON_FALLBACK
            stop_after_attempt = attempt;
#endif
          }
          else
            free_composite(&cand);
        }
#endif
      }
    }
    free_prncplstmnt(&st);
    free_witness(&wt);
#if !SINGLEPROOF_SEARCH_ALL_ATTEMPTS
#if SINGLEPROOF_ENABLE_FALLBACK
    if(have_best || (SINGLEPROOF_STOP_ON_FALLBACK && have_fallback))
      break;
#else
    if(have_best)
      break;
#endif
#endif
#if SINGLEPROOF_SEARCH_ALL_ATTEMPTS
    if(have_best && attempt >= stop_after_attempt)
      break;
#endif
  }
#if !SINGLEPROOF_ENABLE_FALLBACK
  if(!have_best)
    goto no_success;
  p = best;
  memset(&best,0,sizeof(best));
  have_best = 0;
  ret = 0;
  attempt = best_attempt;
#else
  if(have_best && (!have_fallback || best.size <= fallback.size)) {
    p = best;
    memset(&best,0,sizeof(best));
    have_best = 0;
    ret = 0;
    attempt = best_attempt;
  }
  else if(ret && have_fallback) {
    p = fallback;
    memset(&fallback,0,sizeof(fallback));
    have_fallback = 0;
    ret = 0;
    attempt = fallback_attempt;
  }
  else if(have_best) {
    p = best;
    memset(&best,0,sizeof(best));
    have_best = 0;
    ret = 0;
    attempt = best_attempt;
  }
  else if(have_fallback) {
    p = fallback;
    memset(&fallback,0,sizeof(fallback));
    have_fallback = 0;
    ret = 0;
    attempt = fallback_attempt;
  }
#endif

no_success:
  if(ret) {
    fprintf(stderr,"ERROR: Hashed selector singleproof composite proof/verification failed after %d attempts: %d\n",MAX_ATTEMPTS,ret);
    goto end;
  }
  printf("Chihuahua Pack members: %zu\n",p.l);
  printf("Chihuahua Pack size: %.2f KB\n",p.size);
  printf("[singleproof composite] attempt %zu succeeded\n\n",attempt);

end:
  free_prncplstmnt(&st);
  if(have_best)
    free_composite(&best);
#if SINGLEPROOF_ENABLE_FALLBACK
  if(have_fallback)
    free_composite(&fallback);
#endif
  free_composite(&p);
  free_witness(&wt);
  return ret;
}

#ifndef SELECTOR_SINGLEPROOF_NO_MAIN
int main(void) {
  int ret;

  ret = test_singleproof_twolayer();
  if(ret) goto end;
  ret = test_singleproof_pack();
  free_comkey();
end:
  return ret;
}
#endif
