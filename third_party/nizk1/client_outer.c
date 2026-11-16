#define _POSIX_C_SOURCE 200809L
#ifndef LOGQ
#define LOGQ 48
#endif
#ifndef DEMO_L
#define DEMO_L 80
#endif
#ifndef GADGET_BASE
#if LOGQ >= 64
#define GADGET_BASE 319
#else
#define GADGET_BASE 255
#endif
#endif
#ifndef LNP22_BATCH_ROUNDS
#define LNP22_BATCH_ROUNDS 1
#endif
#ifndef OUTER_DEBUG_PRECHECK
#define OUTER_DEBUG_PRECHECK 0
#endif
#ifndef OUTER_PREP_PROFILE
#define OUTER_PREP_PROFILE 0
#endif
#ifndef OUTER_RUN_TWOLAYER
#define OUTER_RUN_TWOLAYER 0
#endif
#ifndef OUTER_VERIFY_REQUEST_BENCH
#define OUTER_VERIFY_REQUEST_BENCH 1
#endif
#ifndef DIRECT_OUTER_VALIDATE
#define DIRECT_OUTER_VALIDATE 0
#endif
#ifndef OUTER_BZ_DIGITS
#define OUTER_BZ_DIGITS 2
#endif
#if OUTER_BZ_DIGITS != 2 && OUTER_BZ_DIGITS != 3
#error "OUTER_BZ_DIGITS must be 2 or 3"
#endif
#ifndef SELECTOR_LNP22_INNER_NO_MAIN
#define SELECTOR_LNP22_INNER_NO_MAIN
#endif
#include "client_inner.c"

#ifndef OUTER_RBIT_PRODUCT_DIGITS
#define OUTER_RBIT_PRODUCT_DIGITS 0
#endif
#if OUTER_RBIT_PRODUCT_DIGITS != 0 && (OUTER_RBIT_PRODUCT_DIGITS < 3 || OUTER_RBIT_PRODUCT_DIGITS > 5)
#error "OUTER_RBIT_PRODUCT_DIGITS must be between 3 and 5"
#endif
#ifndef OUTER_RBIT_PRODUCT_BASE
#define OUTER_RBIT_PRODUCT_BASE 32768
#endif
#ifndef OUTER_RBIT_PRODUCT_CHECK_ROUNDS
#define OUTER_RBIT_PRODUCT_CHECK_ROUNDS 1
#endif

enum {
  OUT_B0 = 0,
  OUT_B1,
#if LNP22_DIRECT_CHAIN
  OUT_X0,
  OUT_X1,
  OUT_P0,
  OUT_P1,
  OUT_P2,
  OUT_P3,
  OUT_P4,
  OUT_T0,
  OUT_T1,
#if LNP22_USE_NATIVE_FOLD_BINARY
  OUT_RAUT0,
  OUT_RAUT1,
#else
  OUT_RBIT0,
  OUT_RBIT1,
#if OUTER_RBIT_PRODUCT_DIGITS > 0
  OUT_PBIT0,
  OUT_PBIT1,
  OUT_PBIT2,
#endif
#if OUTER_RBIT_PRODUCT_DIGITS >= 4
  OUT_PBIT3,
#endif
#if OUTER_RBIT_PRODUCT_DIGITS == 5
  OUT_PBIT4,
#endif
#endif
  OUT_Z20,
  OUT_Z21,
#else
#if OUTER_BZ_DIGITS == 3
  OUT_B2,
#endif
  OUT_X0,
  OUT_X1,
  OUT_Z0,
  OUT_Z1,
#if OUTER_BZ_DIGITS == 3
  OUT_Z2,
#endif
  OUT_Y00,
  OUT_Y01,
  OUT_Y10,
  OUT_Y11,
  OUT_T0,
  OUT_T1,
  OUT_RBIT0,
  OUT_RBIT1,
  OUT_Z20,
  OUT_Z21,
#endif
  OUT_R
};

#ifndef OUTER_BZ_DECOMP_BASE
#if OUTER_BZ_DIGITS == 2
#if LOGQ == 72
#define OUTER_BZ_DECOMP_BASE 40960
#elif LOGQ >= 64
#define OUTER_BZ_DECOMP_BASE 49152
#elif LOGQ >= 56
#define OUTER_BZ_DECOMP_BASE 32768
#else
#define OUTER_BZ_DECOMP_BASE 24576
#endif
#else
#define OUTER_BZ_DECOMP_BASE 256
#endif
#endif

#ifndef OUTER_XY_DECOMP_BASE
#if LOGQ == 72
#define OUTER_XY_DECOMP_BASE 40960
#elif LOGQ >= 64
#define OUTER_XY_DECOMP_BASE 49152
#elif LOGQ >= 56
#define OUTER_XY_DECOMP_BASE 32768
#else
#define OUTER_XY_DECOMP_BASE 24576
#endif
#endif

#ifndef OUTER_Z2_DECOMP_BASE
#define OUTER_Z2_DECOMP_BASE 2048
#endif
#ifndef OUTER_RAUT_ALPHA_ROUNDS
#define OUTER_RAUT_ALPHA_ROUNDS 3
#endif
#ifndef OUTER_DIRECT_P_DIGITS
#define OUTER_DIRECT_P_DIGITS 5
#endif
#if OUTER_DIRECT_P_DIGITS != 5
#error "The direct-chain frontend currently requires five P digits"
#endif
#ifndef OUTER_DIRECT_P_BASE
#if LOGQ > 64
#define OUTER_DIRECT_P_BASE 32768
#else
#define OUTER_DIRECT_P_BASE 8192
#endif
#endif

/* Public admission bounds, not a substitute for the composed MSIS estimate.
 * For z = z0 + base*z1, use |z0| <= base/2 and ||z||^2 <= 2*n*sigma^2. */
static uint64_t outer_pair_budget(size_t coeffs, uint64_t base, int logsigma) {
  const unsigned __int128 denom = (unsigned __int128)base*base;
  const unsigned __int128 high = ((unsigned __int128)16*coeffs << (2*logsigma));
  return (uint64_t)((unsigned __int128)coeffs*denom/4 +
                    (high+denom-1)/denom + coeffs);
}

static uint64_t outer_admissible_betasq(const selector_lnp22_ctx *ctx) {
  const size_t tail = ctx->rhat_len + ctx->umc_eua_len + ctx->umc_eub_len +
                      ctx->umc_r_len + ctx->umc_va_len + ctx->umc_vb_len;
  uint64_t bound = outer_pair_budget(ctx->b_len*N,OUTER_BZ_DECOMP_BASE,LNP22_LOG2STDEV1);
  bound += outer_pair_budget((ctx->x_len+tail+ctx->rhat_bits_len)*N,
                             OUTER_XY_DECOMP_BASE,LNP22_LOG2STDEV1);
  bound += outer_pair_budget(ctx->m2*N,OUTER_Z2_DECOMP_BASE,LNP22_LOG2STDEV2);
  const unsigned __int128 q = ((unsigned __int128)1<<LOGQ)-QOFF;
  unsigned __int128 power = 1;
  for(size_t i=1;i<OUTER_DIRECT_P_DIGITS;i++) power *= OUTER_DIRECT_P_BASE;
  const uint64_t top = (uint64_t)((q/2)/power)+2;
  bound += ctx->x_len*N*((uint64_t)(OUTER_DIRECT_P_DIGITS-1)*
                         OUTER_DIRECT_P_BASE*OUTER_DIRECT_P_BASE/4 + top*top);
#if LNP22_DIRECT_CHAIN && LNP22_USE_RBIT_BINARY && OUTER_RBIT_PRODUCT_DIGITS > 0
  power = 1;
  for(size_t i=1;i<OUTER_RBIT_PRODUCT_DIGITS;i++) power *= OUTER_RBIT_PRODUCT_BASE;
  const unsigned __int128 denom = power*power;
  /* ||beta*z|| <= ||beta||_1 ||z|| <= 2*N ||z||. */
  const unsigned __int128 high = (unsigned __int128)16*(2*N)*(2*N)*
                                 ctx->rhat_bits_len*N*((unsigned __int128)1<<(2*LNP22_LOG2STDEV1));
  bound += ctx->rhat_bits_len*N*((uint64_t)(OUTER_RBIT_PRODUCT_DIGITS-1)*
             OUTER_RBIT_PRODUCT_BASE*OUTER_RBIT_PRODUCT_BASE/4 + 2) +
             (uint64_t)((high+denom-1)/denom);
#endif
  return bound + 16*N;
}

#ifndef OUTER_HASH_EXPANDED_STATEMENT
#define OUTER_HASH_EXPANDED_STATEMENT 0
#endif
#ifndef OUTER_COMPACT_STATEMENT
#define OUTER_COMPACT_STATEMENT 1
#endif
#define OUTER_USE_COMPACT_STATEMENT (OUTER_COMPACT_STATEMENT && !OUTER_HASH_EXPANDED_STATEMENT)
#if OUTER_HASH_EXPANDED_STATEMENT
#define OUTER_STMT_ABSORB_POLX(st, p) absorb_polx((st),(p))
#else
#define OUTER_STMT_ABSORB_POLX(st, p) do { (void)(st); (void)(p); } while(0)
#endif
typedef struct {
  selector_lnp22_ctx ctx;
  polx *s2;
  polx *tA;
  polx *a3;
  polx *w;
  polx *z1;
  polx *z2;
#if !LNP22_COEFF_RESPONSE
  polx *cs1;
  polx *cs2;
#endif
  polx *y1;
  polx *y2;
  int64_t *z1coeffs;
  int64_t *z2coeffs;
  int64_t *cs1coeffs;
  int64_t *cs2coeffs;
  polx c[1];
  polx vlin[1];
  polx tq[1];
  polx vq[1];
  polx ufold_public[1];
  polx uquad_public[1];
  polx a3s2[1];
  polx g1[1];
  polx g0[1];
  size_t rejection_rounds;
} outer_instance;

typedef struct {
  size_t count;
  polx *owned[12];
} outer_frontend_storage;

static void free_outer_frontend_storage(void *ptr) {
  outer_frontend_storage *storage = ptr;
  size_t i;

  if(!storage)
    return;
  for(i=0;i<storage->count;i++)
    free(storage->owned[i]);
  free(storage);
}

static void attach_outer_frontend_storage(prncplstmnt *st, polx *owned[12]) {
  outer_frontend_storage *storage = _malloc(sizeof(*storage));
  size_t i;

  storage->count = 12;
  for(i=0;i<storage->count;i++)
    storage->owned[i] = owned[i];
  st->frontend_data = storage;
  st->free_frontend_data = free_outer_frontend_storage;
}

static int64_t balanced_round_div_i64(int64_t x, int64_t base) {
  if(x >= 0)
    return (x + base/2) / base;
  return -(((-x) + base/2) / base);
}

static uint64_t outer_abs_i64(int64_t x) {
  return (x < 0) ? (uint64_t)(-x) : (uint64_t)x;
}

static int decompose_segment_two_digits(poly *dst0, poly *dst1,
                                        const int64_t *src, size_t len,
                                        int64_t base, const char *label) {
  const size_t count = len*N;
  size_t i;
  int invalid = 0;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(count >= LNP22_PARALLEL_MIN_LEN) reduction(|:invalid)
#endif
  for(i=0;i<count;i++) {
    const int64_t q = balanced_round_div_i64(src[i],base);
    const int64_t r = src[i] - base*q;

    invalid |= outer_abs_i64(r) > 32767 || outer_abs_i64(q) > 32767;
    dst0[i/N].vec[0].c[i%N] = (int16_t)r;
    dst1[i/N].vec[0].c[i%N] = (int16_t)q;
  }
  if(invalid)
    fprintf(stderr,"ERROR: %s decomposition digit out of int16 range\n",label);
  return invalid;
}

static int decompose_segment_three_digits(poly *dst0, poly *dst1, poly *dst2,
                                          const int64_t *src, size_t len,
                                          int64_t base, const char *label) {
  const size_t count = len*N;
  size_t i;
  int invalid = 0;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(count >= LNP22_PARALLEL_MIN_LEN) reduction(|:invalid)
#endif
  for(i=0;i<count;i++) {
    const int64_t q1 = balanced_round_div_i64(src[i],base);
    const int64_t r0 = src[i] - base*q1;
    const int64_t q2 = balanced_round_div_i64(q1,base);
    const int64_t r1 = q1 - base*q2;

    invalid |= outer_abs_i64(r0) > 32767 || outer_abs_i64(r1) > 32767 ||
               outer_abs_i64(q2) > 32767;
    dst0[i/N].vec[0].c[i%N] = (int16_t)r0;
    dst1[i/N].vec[0].c[i%N] = (int16_t)r1;
    dst2[i/N].vec[0].c[i%N] = (int16_t)q2;
  }
  if(invalid)
    fprintf(stderr,"ERROR: %s decomposition digit out of int16 range\n",label);
  return invalid;
}

static int decompose_segment_five_digits(poly *dst[OUTER_DIRECT_P_DIGITS],
                                         const __int128 *src, size_t len,
                                         int64_t base, const char *label) {
  const size_t count = len*N;
  size_t i;
  int invalid = 0;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(count >= LNP22_PARALLEL_MIN_LEN) reduction(|:invalid)
#endif
  for(i=0;i<count;i++) {
    __int128 q = src[i];
    __int128 digit[OUTER_DIRECT_P_DIGITS];
    size_t j;

    for(j=0;j<OUTER_DIRECT_P_DIGITS-1;j++) {
      const __int128 next = q >= 0 ? (q+base/2)/base : -((-q+base/2)/base);
      digit[j] = q - base*next;
      q = next;
    }
    digit[OUTER_DIRECT_P_DIGITS-1] = q;
    for(j=0;j<OUTER_DIRECT_P_DIGITS;j++) {
      invalid |= digit[j] < -32767 || digit[j] > 32767;
      dst[j][i/N].vec[0].c[i%N] = (int16_t)digit[j];
    }
  }
  if(invalid)
    fprintf(stderr,"ERROR: %s decomposition digit out of int16 range\n",label);
  return invalid;
}

static int decompose_rbit_product_digits(
    poly *dst[OUTER_RBIT_PRODUCT_DIGITS], const int64_t *src,
    size_t len, const char *label) {
#if OUTER_RBIT_PRODUCT_DIGITS == 0
  (void)dst; (void)src; (void)len; (void)label;
  return 0;
#else
  const size_t count = len*N;
  size_t i;
  int invalid = 0;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(count >= LNP22_PARALLEL_MIN_LEN) reduction(|:invalid)
#endif
  for(i=0;i<count;i++) {
    int64_t q = src[i];
    int64_t digit[OUTER_RBIT_PRODUCT_DIGITS];
    size_t j;

    for(j=0;j<OUTER_RBIT_PRODUCT_DIGITS-1;j++) {
      const int64_t next =
          balanced_round_div_i64(q,OUTER_RBIT_PRODUCT_BASE);
      digit[j] = q - OUTER_RBIT_PRODUCT_BASE*next;
      q = next;
    }
    digit[OUTER_RBIT_PRODUCT_DIGITS-1] = q;
    for(j=0;j<OUTER_RBIT_PRODUCT_DIGITS;j++) {
      invalid |= outer_abs_i64(digit[j]) > 32767;
      dst[j][i/N].vec[0].c[i%N] = (int16_t)digit[j];
    }
  }
  if(invalid)
    fprintf(stderr,"ERROR: %s decomposition digit out of int16 range\n",label);
  return invalid;
#endif
}

static void set_scaled_polxvec(polx *dst, const polx *src, size_t len, int64_t scale) {
  size_t i;

  if(scale == 1) {
    polxvec_copy(dst,src,len);
    return;
  }

#if defined(_OPENMP) && LNP22_ROW_PARALLEL && LNP22_PARALLEL_MIN_LEN > 0
#pragma omp parallel for schedule(static) if(len >= LNP22_PARALLEL_MIN_LEN)
#elif defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(i=0;i<len;i++)
    polx_scale(&dst[i],&src[i],scale);
}

static polx *append_cnst_block(sparsecnst *cnst, size_t *pos, polx *buf,
                               size_t idx, size_t len) {
  const size_t p = *pos;

  cnst->idx[p] = idx;
  cnst->off[p] = 0;
  cnst->len[p] = len;
  cnst->mult[p] = 1;
  cnst->phi[p] = buf;
  *pos = p + 1;
  return buf + len;
}

static void append_cnst_view(sparsecnst *cnst, size_t *pos, size_t idx,
                             size_t len, polx *src, int64_t scale) {
  const size_t p = *pos;

  cnst->idx[p] = idx;
  cnst->off[p] = 0;
  cnst->len[p] = len;
  cnst->mult[p] = 1;
  cnst->phi[p] = src;
  cnst->scale[p] = scale;
  *pos = p + 1;
}

static void append_cnst_offset_view(sparsecnst *cnst, size_t *pos, size_t idx,
                                    size_t off, size_t len, polx *src,
                                    int64_t scale) {
  const size_t p = *pos;

  cnst->idx[p] = idx;
  cnst->off[p] = off;
  cnst->len[p] = len;
  cnst->mult[p] = 1;
  cnst->phi[p] = src;
  cnst->scale[p] = scale;
  *pos = p + 1;
}

static int polxvec_equal_slow(const polx *a, const polx *b, size_t len) {
  size_t i;
  polx tmp[1];

  for(i=0;i<len;i++) {
    polx_sub(tmp,&a[i],&b[i]);
    if(!polx_iszero(tmp))
      return 0;
  }
  return 1;
}

static int64_t outer_pow_i64(int64_t base, int exp) {
  int64_t s = 1;
  int i;

  for(i=0;i<exp;i++)
    s *= base;
  return s;
}

static void derive_outer_openfold(polx *openfold, const selector_lnp22_ctx *ctx,
                                  const polx *tA, const polx *w, const polx *vlin,
                                  const polx *tq, const polx *vq, const polx *c) {
  shake128incctx st;
  uint8_t seed[16];

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-outer-openfold",23);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,ctx->kmsis);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  absorb_polxvec(&st,ctx->Cout,ctx->m);
  absorb_rbit_const_transcript(&st,ctx);
  absorb_polxvec(&st,tA,ctx->kmsis);
  absorb_polxvec(&st,w,ctx->kmsis);
  absorb_polx(&st,vlin);
  absorb_polx(&st,tq);
  absorb_polx(&st,vq);
  absorb_polx(&st,c);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(seed,sizeof(seed),&st);
  polxvec_challenge(openfold,ctx->kmsis,seed,0);
}


static void absorb_outer_compact_statement(shake128incctx *st, const selector_lnp22_ctx *ctx,
                                           const polx *tA, const polx *w, const polx *vlin,
                                           const polx *tq, const polx *vq, const polx *c) {
  absorb_bytes(st,"selector-outer-statement-compact",sizeof("selector-outer-statement-compact")-1);
  absorb_u64(st,LOGQ);
  absorb_u64(st,ctx->Lchain);
  absorb_u64(st,GADGET_BASE);
  absorb_u64(st,ctx->kmsis);
  absorb_u64(st,ctx->m2);
  absorb_u64(st,OUTER_BZ_DIGITS);
  absorb_u64(st,OUTER_BZ_DECOMP_BASE);
  absorb_u64(st,OUTER_XY_DECOMP_BASE);
  absorb_u64(st,OUTER_Z2_DECOMP_BASE);
  absorb_u64(st,LNP22_DIRECT_CHAIN);
  absorb_u64(st,OUTER_DIRECT_P_DIGITS);
  absorb_u64(st,OUTER_DIRECT_P_BASE);
  absorb_u64(st,LNP22_RHAT_BINARY_BACKEND);
  absorb_polxvec(st,ctx->CX,ctx->cx_rows);
  absorb_polxvec(st,ctx->Cout,ctx->m);
  absorb_rbit_const_transcript(st,ctx);
  absorb_polxvec(st,tA,ctx->kmsis);
  absorb_polxvec(st,w,ctx->kmsis);
  absorb_polx(st,vlin);
  absorb_polx(st,tq);
  absorb_polx(st,vq);
  absorb_polx(st,c);
}

static void compute_openfold_coeffs(polx *open1, polx *open2, polx *rhs,
                                    const selector_lnp22_ctx *ctx,
                                    const polx *tA, const polx *w, const polx *c,
                                    const polx *openfold) {
  size_t row;
  polx rhsrow[1], tmp[1];

  polxvec_setzero(open1,ctx->m1);
  polxvec_setzero(open2,ctx->m2);
  polxvec_setzero(rhs,1);

  for(row=0;row<ctx->kmsis;row++) {
    polx_mul(rhsrow,c,&tA[row]);
    polx_add(rhsrow,rhsrow,&w[row]);
    polx_mul(tmp,&openfold[row],rhsrow);
    polx_add(rhs,rhs,tmp);

    polxvec_polx_mul_add(open2,&openfold[row],&ctx->A2open[row*ctx->m2prime],ctx->m2prime);
    polx_add(&open2[ctx->m2prime + row],&open2[ctx->m2prime + row],&openfold[row]);
  }

  {
    const size_t chunks = (ctx->m1 + 31)/32;
    size_t blk;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
    for(blk=0;blk<chunks;blk++) {
      const size_t off = 32*blk;
      const size_t len = MIN((size_t)32,ctx->m1 - off);
      size_t rrow;
      polx rowchunk[32];

      for(rrow=0;rrow<ctx->kmsis;rrow++) {
        if(ctx->A1open != NULL) {
          polxvec_polx_mul_add(&open1[off],&openfold[rrow],
                               &ctx->A1open[rrow*ctx->m1 + off],len);
        } else {
          const size_t row_base = rrow*chunks;
          polxvec_almostuniform(rowchunk,len,ctx->lnp_seed,NONCE_LNP_A1_BASE + row_base + blk);
          polxvec_polx_mul_add(&open1[off],&openfold[rrow],rowchunk,len);
        }
      }
    }
  }
}

static void compute_quad_linear_coeffs(polx *zcoeff, polx *xcoeff,
                                       polx *rbitcoeff, polx *z2coeff,
                                       const selector_lnp22_ctx *ctx,
                                       const polx *a3, const polx *c) {
  size_t t, i;
  polx tmp[1];

  polxvec_setzero(zcoeff,ctx->z_len);
  for(t=0;t<ctx->Lchain-1;t++) {
    for(i=0;i<ctx->gdigits;i++) {
      polx_mul(tmp,&ctx->alphavec[t],&ctx->kvec[t*ctx->gdigits + i]);
      polx_mul(tmp,tmp,c);
      polx_add(&zcoeff[t*ctx->gdigits + i],&zcoeff[t*ctx->gdigits + i],tmp);
    }
  }

  polxvec_setzero(xcoeff,ctx->x_len);
  if(LNP22_PROVE_X_BITS) {
    for(i=0;i<ctx->x_len;i++)
      polx_neg(&xcoeff[i],c);
  }

  polxvec_setzero(rbitcoeff,ctx->rhat_bits_len);
  for(i=0;i<ctx->rhat_bits_len;i++)
    polx_neg(&rbitcoeff[i],c);

  polxvec_copy(z2coeff,a3,ctx->m2);
}

static int generate_inner_instance(outer_instance *inst, uint64_t seed_offset) {
  uint8_t seeded_demo_seed[16];
  int ok = 1;
  lnp22_round_ws *batch_ws = NULL;
  int64_t *s1coeffs = NULL;
  int64_t *s2coeffs = NULL;
  int64_t ccoeffs[N];
#if defined(_OPENMP)
  const size_t batch_rounds = LNP22_ENABLE_REJECTION ? LNP22_BATCH_ROUNDS : 1;
#else
  const size_t batch_rounds = 1;
#endif
  size_t batch_workers = 1;
#if OUTER_PREP_PROFILE
  struct timespec gen_prof_last, gen_prof_now;
#define GEN_PROFILE_START() clock_gettime(CLOCK_MONOTONIC,&gen_prof_last)
#define GEN_PROFILE_MARK(label) do {                                        \
    clock_gettime(CLOCK_MONOTONIC,&gen_prof_now);                           \
    printf("    inner gen profile: %-22s %.6fs\n",label,                    \
           elapsed_sec(&gen_prof_last,&gen_prof_now));                      \
    gen_prof_last = gen_prof_now;                                            \
  } while(0)
#else
#define GEN_PROFILE_START() do { } while(0)
#define GEN_PROFILE_MARK(label) do { } while(0)
#endif

  GEN_PROFILE_START();
  memset(inst,0,sizeof(*inst));
  add_seed_offset(seeded_demo_seed,demo_fixed_seed,seed_offset);

  selector_ctx_init(&inst->ctx,seeded_demo_seed);
  GEN_PROFILE_MARK("ctx init");
  prepare_selector_lnp22(&inst->ctx);
  GEN_PROFILE_MARK("prepare stmt");
#if LNP22_VALIDATE_BASE_WITNESS
  if(first_linear_relation_failure(&inst->ctx,inst->ctx.s1)) {
    fprintf(stderr,"ERROR: inner base witness does not satisfy intended relations\n");
    return 1;
  }
#endif
  GEN_PROFILE_MARK("basecheck");

  inst->s2 = alloc_polxvec(inst->ctx.m2);
  inst->tA = alloc_polxvec(inst->ctx.kmsis);
  inst->a3 = alloc_polxvec(inst->ctx.m2);
  if(batch_rounds <= 1) {
    inst->w = alloc_polxvec(inst->ctx.kmsis);
    inst->z1 = alloc_polxvec(inst->ctx.m1);
    inst->z2 = alloc_polxvec(inst->ctx.m2);
#if !LNP22_COEFF_RESPONSE
    inst->cs1 = alloc_polxvec(inst->ctx.m1);
    inst->cs2 = alloc_polxvec(inst->ctx.m2);
#endif
    inst->y1 = alloc_polxvec(inst->ctx.m1);
    inst->y2 = alloc_polxvec(inst->ctx.m2);
    inst->z1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->z1coeffs));
    inst->z2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->z2coeffs));
    inst->cs1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->cs1coeffs));
    inst->cs2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->cs2coeffs));
  }
#if LNP22_COEFF_RESPONSE
  s1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*s1coeffs));
  s2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*s2coeffs));
  witness_coeffs_to_i64(s1coeffs,&inst->ctx);
#endif

  if(batch_rounds > 1) {
    size_t bi;

#if defined(_OPENMP)
    batch_workers = (size_t)omp_get_max_threads();
#endif
    if(batch_workers > batch_rounds)
      batch_workers = batch_rounds;
    batch_ws = _malloc(batch_rounds * sizeof(*batch_ws));
    for(bi=0;bi<batch_rounds;bi++)
      round_ws_init(&batch_ws[bi],&inst->ctx);
  }
  GEN_PROFILE_MARK("alloc");

  sample_ternary_polxvec(inst->s2,inst->ctx.m2,inst->ctx.lnp_seed,NONCE_S2_BASE);
#if LNP22_COEFF_RESPONSE
  polxvec_centered_to_i64(s2coeffs,inst->s2,inst->ctx.m2);
#endif
  compute_a3(inst->a3,&inst->ctx);
  inner_product_a3(inst->a3s2,inst->a3,inst->s2,inst->ctx.m2);
  GEN_PROFILE_MARK("s2+a3");
  if(LNP22_DIRECT_CHAIN) {
    compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);
    bind_prechallenge_commitment(&inst->ctx,inst->tA);
  }
  if(prepare_linear_fold_for_proof(&inst->ctx,inst->ufold_public,inst->s2)) {
    fprintf(stderr,"ERROR: failed to prepare inner folded linear relation\n");
    return 1;
  }
  compute_public_quad_rhs(inst->uquad_public,&inst->ctx);
  GEN_PROFILE_MARK("linear setup");
  if(batch_rounds > 1 && !LNP22_DIRECT_CHAIN)
    trim_ctx_for_batch(&inst->ctx);
#if LNP22_COMBINE_TA_W
  if(!LNP22_DIRECT_CHAIN && (batch_rounds > 1 || batch_ws != NULL))
    compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);
#else
  if(!LNP22_DIRECT_CHAIN)
    compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);
#endif
  if(LNP22_DIRECT_CHAIN) {
    polx direct_check[1];

    derive_direct_chain_fold(&inst->ctx,inst->tA);
    compute_quad_g1_coeffs(&inst->ctx);
    quad_value(direct_check,&inst->ctx,inst->ctx.s1);
    if(LNP22_USE_NATIVE_FOLD_BINARY) {
      polx native_term[1];

      polxvec_sprod(native_term,inst->ctx.native_binary_lin2,
                    inst->s2,inst->ctx.m2);
      polx_add(direct_check,direct_check,native_term);
      polx_sub(direct_check,direct_check,inst->uquad_public);
    }
    if(!polx_iszero(direct_check)) {
      fprintf(stderr,"ERROR: direct-chain base quadratic relation failed\n");
      return 1;
    }
    if(batch_rounds > 1)
      trim_ctx_for_batch(&inst->ctx);
  }
  GEN_PROFILE_MARK("compute tA");

  if(batch_rounds <= 1 || batch_ws == NULL) {
    for(inst->rejection_rounds = 0;;inst->rejection_rounds++) {
      size_t i;
      int rej1, rej2;

      sample_gaussian_polxvec(inst->y1,inst->z1coeffs,inst->ctx.m1,LNP22_LOG2STDEV1,
                              inst->ctx.lnp_seed,NONCE_Y1_BASE + inst->rejection_rounds);
      sample_gaussian_polxvec(inst->y2,inst->z2coeffs,inst->ctx.m2,LNP22_LOG2STDEV2,
                              inst->ctx.lnp_seed,NONCE_Y2_BASE + inst->rejection_rounds);
#if LNP22_COMBINE_TA_W
      if(LNP22_DIRECT_CHAIN)
        compute_w_rows(inst->w,&inst->ctx,inst->y1,inst->y2);
      else
        compute_tA_w_rows(inst->tA,inst->w,&inst->ctx,inst->ctx.s1,inst->s2,
                          inst->y1,inst->y2);
#else
      compute_w_rows(inst->w,&inst->ctx,inst->y1,inst->y2);
#endif
      folded_linear_eval(inst->vlin,&inst->ctx,inst->y1,inst->y2);
      quad_mask_terms(inst->g1,inst->g0,&inst->ctx,inst->ctx.s1,inst->y1);
      add_native_binary_s2_mask(inst->g1,&inst->ctx,inst->y2);
      polx_add(inst->tq,inst->a3s2,inst->g1);
      inner_product_a3(inst->vq,inst->a3,inst->y2,inst->ctx.m2);
      polx_add(inst->vq,inst->vq,inst->g0);
      fs_challenge(inst->c,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
#if LNP22_COEFF_RESPONSE
      polxvec_centered_to_i64(ccoeffs,inst->c,1);
      challenge_mul_coeffs(inst->cs1coeffs,ccoeffs,s1coeffs,inst->ctx.m1);
      challenge_mul_coeffs(inst->cs2coeffs,ccoeffs,s2coeffs,inst->ctx.m2);
#else
      polxvec_polx_mul(inst->cs1,inst->c,inst->ctx.s1,inst->ctx.m1);
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_polx_mul(inst->cs2,inst->c,inst->s2,inst->ctx.m2);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);

      polxvec_centered_to_i64(inst->cs1coeffs,inst->cs1,inst->ctx.m1);
      polxvec_centered_to_i64(inst->cs2coeffs,inst->cs2,inst->ctx.m2);
#endif
      for(i=0;i<inst->ctx.m1*N;i++)
        inst->z1coeffs[i] += inst->cs1coeffs[i];
      for(i=0;i<inst->ctx.m2*N;i++)
        inst->z2coeffs[i] += inst->cs2coeffs[i];

      if(LNP22_ENABLE_REJECTION) {
        rej1 = reject_standard_polxvec(inst->z1coeffs,inst->cs1coeffs,inst->ctx.m1*N,
                                       sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                       LNP22_GAMMA1,
                                       inst->ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + inst->rejection_rounds);
        rej2 = reject_z2_polxvec(inst->z2coeffs,inst->cs2coeffs,inst->ctx.m2*N,
                                   sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                   LNP22_GAMMA2,
                                   inst->ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + inst->rejection_rounds);
      } else {
        rej1 = 0;
        rej2 = 0;
      }
      if(!rej1 && !rej2) {
#if LNP22_COEFF_RESPONSE
        polxvec_fromint64vec(inst->z1,inst->ctx.m1,1,inst->z1coeffs);
        polxvec_fromint64vec(inst->z2,inst->ctx.m2,1,inst->z2coeffs);
#endif
        break;
      }
    }
  } else {
    uint8_t *accept = _malloc(batch_rounds);
    size_t batch_base = 0;
    size_t found_bi = 0;
    int found = 0;

    while(!found) {
      size_t bi;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) num_threads(batch_workers)
#endif
      for(bi=0;bi<batch_rounds;bi++) {
        lnp22_round_ws *ws = &batch_ws[bi];
        size_t i;

        sample_gaussian_polxvec(ws->y1,ws->z1coeffs,inst->ctx.m1,LNP22_LOG2STDEV1,
                                inst->ctx.lnp_seed,NONCE_Y1_BASE + batch_base + bi);
        sample_gaussian_polxvec(ws->y2,ws->z2coeffs,inst->ctx.m2,LNP22_LOG2STDEV2,
                                inst->ctx.lnp_seed,NONCE_Y2_BASE + batch_base + bi);
        compute_w_rows(ws->w,&inst->ctx,ws->y1,ws->y2);
        folded_linear_eval(ws->vlin,&inst->ctx,ws->y1,ws->y2);
        quad_mask_terms(ws->g1,ws->g0,&inst->ctx,inst->ctx.s1,ws->y1);
        add_native_binary_s2_mask(ws->g1,&inst->ctx,ws->y2);
        polx_add(ws->tq,inst->a3s2,ws->g1);
        inner_product_a3(ws->vq,inst->a3,ws->y2,inst->ctx.m2);
        polx_add(ws->vq,ws->vq,ws->g0);
        fs_challenge(ws->c,&inst->ctx,inst->tA,ws->w,ws->vlin,ws->tq,ws->vq);
#if LNP22_COEFF_RESPONSE
        polxvec_centered_to_i64(ws->ccoeffs,ws->c,1);
        challenge_mul_coeffs(ws->cs1coeffs,ws->ccoeffs,s1coeffs,inst->ctx.m1);
        challenge_mul_coeffs(ws->cs2coeffs,ws->ccoeffs,s2coeffs,inst->ctx.m2);
#else
        polxvec_polx_mul(ws->cs1,ws->c,inst->ctx.s1,inst->ctx.m1);
        polxvec_polx_mul(ws->cs2,ws->c,inst->s2,inst->ctx.m2);
        polxvec_centered_to_i64(ws->cs1coeffs,ws->cs1,inst->ctx.m1);
        polxvec_centered_to_i64(ws->cs2coeffs,ws->cs2,inst->ctx.m2);
#endif
        for(i=0;i<inst->ctx.m1*N;i++)
          ws->z1coeffs[i] += ws->cs1coeffs[i];
        for(i=0;i<inst->ctx.m2*N;i++)
          ws->z2coeffs[i] += ws->cs2coeffs[i];
        ws->rej1 = reject_standard_polxvec(ws->z1coeffs,ws->cs1coeffs,inst->ctx.m1*N,
                                           sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                           LNP22_GAMMA1,
                                           inst->ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + batch_base + bi);
        ws->rej2 = reject_z2_polxvec(ws->z2coeffs,ws->cs2coeffs,inst->ctx.m2*N,
                                       sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                       LNP22_GAMMA2,
                                       inst->ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + batch_base + bi);
        accept[bi] = (uint8_t)(!ws->rej1 && !ws->rej2);
      }

      for(bi=0;bi<batch_rounds;bi++) {
        if(accept[bi]) {
          inst->rejection_rounds = batch_base + bi;
          found_bi = bi;
          found = 1;
          break;
        }
      }
      batch_base += batch_rounds;
    }

    {
      const lnp22_round_ws *winner = &batch_ws[found_bi];
      size_t i;

      inst->w = alloc_polxvec(inst->ctx.kmsis);
      inst->z1 = alloc_polxvec(inst->ctx.m1);
      inst->z2 = alloc_polxvec(inst->ctx.m2);
#if !LNP22_COEFF_RESPONSE
      inst->cs1 = alloc_polxvec(inst->ctx.m1);
      inst->cs2 = alloc_polxvec(inst->ctx.m2);
#endif
      inst->y1 = alloc_polxvec(inst->ctx.m1);
      inst->y2 = alloc_polxvec(inst->ctx.m2);
      inst->z1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->z1coeffs));
      inst->z2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->z2coeffs));
      inst->cs1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->cs1coeffs));
      inst->cs2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->cs2coeffs));
      polxvec_copy(inst->y1,winner->y1,inst->ctx.m1);
      polxvec_copy(inst->y2,winner->y2,inst->ctx.m2);
      polxvec_copy(inst->w,winner->w,inst->ctx.kmsis);
#if !LNP22_COEFF_RESPONSE
      polxvec_copy(inst->cs1,winner->cs1,inst->ctx.m1);
      polxvec_copy(inst->cs2,winner->cs2,inst->ctx.m2);
#endif
      polxvec_copy(inst->c,winner->c,1);
      polxvec_copy(inst->vlin,winner->vlin,1);
      polxvec_copy(inst->tq,winner->tq,1);
      polxvec_copy(inst->vq,winner->vq,1);
      polxvec_copy(inst->g1,winner->g1,1);
      polxvec_copy(inst->g0,winner->g0,1);
      for(i=0;i<inst->ctx.m1*N;i++) {
        inst->z1coeffs[i] = winner->z1coeffs[i];
        inst->cs1coeffs[i] = winner->cs1coeffs[i];
      }
      for(i=0;i<inst->ctx.m2*N;i++) {
        inst->z2coeffs[i] = winner->z2coeffs[i];
        inst->cs2coeffs[i] = winner->cs2coeffs[i];
      }
#if LNP22_COEFF_RESPONSE
      polxvec_fromint64vec(inst->z1,inst->ctx.m1,1,inst->z1coeffs);
      polxvec_fromint64vec(inst->z2,inst->ctx.m2,1,inst->z2coeffs);
#else
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);
#endif
    }
    free(accept);
  }
  GEN_PROFILE_MARK("round search");

#if OUTER_DEBUG_PRECHECK
  {
    polx c_check[1], lhs[1], rhs[1], tmp[1];
    size_t i;
    fs_challenge(c_check,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
    if(memcmp(inst->c,c_check,sizeof(polx)) != 0) {
      fprintf(stderr,"ERROR: inner transcript challenge mismatch\n");
      ok = 0;
    }

    for(i=0;i<inst->ctx.kmsis && ok;i++) {
      polx row_eval[1];
      if(inst->ctx.cache_a1open)
        polxvec_sprod(row_eval,&inst->ctx.A1open[i*inst->ctx.m1],inst->z1,inst->ctx.m1);
      else
        stream_rowdot_almostuniform(row_eval,inst->ctx.lnp_seed,NONCE_LNP_A1_BASE,
                                    i,inst->ctx.m1,inst->z1);
      polxvec_sprod(tmp,&inst->ctx.A2open[i*inst->ctx.m2prime],inst->z2,inst->ctx.m2prime);
      polx_add(row_eval,row_eval,tmp);
      polx_add(row_eval,row_eval,&inst->z2[inst->ctx.m2prime + i]);
      polx_mul(rhs,inst->c,&inst->tA[i]);
      polx_sub(row_eval,row_eval,rhs);
      polx_sub(row_eval,row_eval,&inst->w[i]);
      if(!polx_iszero(row_eval)) {
        fprintf(stderr,"ERROR: inner opening check failed at row %zu\n",i);
        ok = 0;
      }
    }

    if(ok) {
      folded_linear_eval(lhs,&inst->ctx,inst->z1,inst->z2);
      polx_mul(rhs,inst->c,inst->ufold_public);
      polx_sub(lhs,lhs,rhs);
      polx_sub(lhs,lhs,inst->vlin);
      if(!polx_iszero(lhs)) {
        fprintf(stderr,"ERROR: inner folded linear check failed\n");
        ok = 0;
      }
    }

    if(ok) {
      polx lin_term[1], bil_term[1], bits_term[1], mask_term[1], a3z2[1], c2[1], c2u[1];
      linZ_value(lin_term,&inst->ctx,&inst->z1[inst->ctx.off[LNP22_Z_BLOCK]]);
      polx_mul(lin_term,inst->c,lin_term);
      bilinear_value(bil_term,&inst->ctx,&inst->z1[inst->ctx.off[LNP22_X_BLOCK]],
                     &inst->z1[inst->ctx.off[LNP22_Y0_BLOCK]],
                     &inst->z1[inst->ctx.off[LNP22_Y1_BLOCK]]);
      polx_sub(lhs,lin_term,bil_term);
      if(LNP22_PROVE_X_BITS) {
        x_bits_verifier_value(bits_term,&inst->ctx,
                              &inst->z1[inst->ctx.off[LNP22_X_BLOCK]],
                              inst->c);
        polx_add(lhs,lhs,bits_term);
      }
      rhat_bits_verifier_value(bits_term,&inst->ctx,
                               &inst->z1[inst->ctx.off[LNP22_RHAT_BITS_BLOCK]],
                               inst->c);
      polx_add(lhs,lhs,bits_term);
      polx_mul(c2,inst->c,inst->c);
      polx_mul(c2u,c2,inst->uquad_public);
      polx_sub(lhs,lhs,c2u);
      inner_product_a3(a3z2,inst->a3,inst->z2,inst->ctx.m2);
      polx_mul(mask_term,inst->c,inst->tq);
      polx_sub(mask_term,mask_term,a3z2);
      polx_sub(lhs,lhs,mask_term);
      polx_sub(lhs,lhs,inst->vq);
      if(!polx_iszero(lhs)) {
        fprintf(stderr,"ERROR: inner quadratic check failed\n");
        ok = 0;
      }
    }
  }
#endif

  if(batch_ws != NULL) {
    size_t bi;

    for(bi=0;bi<batch_rounds;bi++)
      round_ws_clear(&batch_ws[bi]);
    free(batch_ws);
  }
  GEN_PROFILE_MARK("batch cleanup");
  free(s1coeffs);
  free(s2coeffs);

#undef GEN_PROFILE_START
#undef GEN_PROFILE_MARK
  return ok ? 0 : 1;
}

static void release_inner_secret_buffers(outer_instance *inst) {
#define RELEASE(p) do { free(p); (p)=NULL; } while(0)
  RELEASE(inst->z1);
  RELEASE(inst->z2);
#if !LNP22_COEFF_RESPONSE
  RELEASE(inst->cs1);
  RELEASE(inst->cs2);
#endif
  RELEASE(inst->y1);
  RELEASE(inst->y2);
  RELEASE(inst->z1coeffs);
  RELEASE(inst->z2coeffs);
  RELEASE(inst->cs1coeffs);
  RELEASE(inst->cs2coeffs);
  RELEASE(inst->ctx.s1);
  RELEASE(inst->ctx.g1coeff);
  for(size_t i=0;i<LNP22_NUM_BLOCKS;i++) RELEASE(inst->ctx.witness[i]);
#undef RELEASE
}

static void free_outer_instance(outer_instance *inst) {
  free(inst->s2);
  free(inst->tA);
  free(inst->a3);
  free(inst->w);
  free(inst->z1);
  free(inst->z2);
#if !LNP22_COEFF_RESPONSE
  free(inst->cs1);
  free(inst->cs2);
#endif
  free(inst->y1);
  free(inst->y2);
  free(inst->z1coeffs);
  free(inst->z2coeffs);
  free(inst->cs1coeffs);
  free(inst->cs2coeffs);
  selector_ctx_clear(&inst->ctx);
  memset(inst,0,sizeof(*inst));
}


static void prepare_selector_lnp22_verifier_static(selector_lnp22_ctx *ctx) {
  poly *g0poly = alloc_polyvec(ctx->gdigits);
  poly *g1poly = alloc_polyvec(ctx->gdigits);
  uint64_t nonce = 0;
  size_t row, j;

  if(!g0poly || !g1poly) {
    fprintf(stderr,"ERROR: verifier public allocation failure\n");
    exit(1);
  }

  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&g0poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,g0poly,ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);

  sample_umc_matrices(ctx);
  fill_almostuniform_matrix(ctx->A2open,ctx->kmsis,ctx->m2prime,
                            ctx->lnp_seed,NONCE_LNP_A2_BASE);

  free(g0poly);
  free(g1poly);
}

static int prepare_selector_lnp22_verifier_request(selector_lnp22_ctx *ctx,
                                                   const selector_lnp22_ctx *pub,
                                                   polx *ufold_public, const polx *tA) {
  uint8_t ar_seed[16];
  size_t t;

  ctx->selected_j = pub->selected_j;
  polxvec_copy(ctx->CX,pub->CX,ctx->cx_rows);
  polxvec_copy(ctx->Cout,pub->Cout,ctx->m);
  derive_ar_seed_from_commitment(ar_seed,ctx,ctx->CX);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,ar_seed,NONCE_AR_BASE);

  derive_selector_chal_seed_from_public(ctx);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(t=0;t<ctx->Lchain-1;t++) {
    if(!LNP22_DIRECT_CHAIN) {
      polxvec_challenge(&ctx->kvec[t*ctx->gdigits],ctx->gdigits,
                        ctx->selector_chal_seed,(uint64_t)t);
      polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,
                        ctx->selector_chal_seed,32 + (uint64_t)t);
    }
  }
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->selector_chal_seed,64);
  polxvec_challenge(ctx->alphavec,ctx->Lchain-1,ctx->selector_chal_seed,65);

  bind_prechallenge_commitment(ctx,tA);

  if(ctx->rbit_const_lambda != 0) {
    polxvec_almostuniform(ctx->rbit_const_Bg,ctx->rbit_const_lambda * ctx->m2,
                          ctx->selector_chal_seed,NONCE_RBIT_CONST_BG_BASE);
    polxvec_copy(ctx->rbit_const_tg,pub->rbit_const_tg,ctx->rbit_const_lambda);
    polxvec_copy(ctx->rbit_const_h,pub->rbit_const_h,ctx->rbit_const_lambda);
    if(validate_rbit_const_transcript(ctx))
      return 1;
    derive_rbit_const_seed(ctx);
    if(LNP22_USE_NATIVE_FOLD_BINARY) {
      for(t=0;t<ctx->rbit_const_lambda;t++)
        sample_full_scalar_challenge(&ctx->native_binary_gamma[t],
                                     ctx->rbit_const_seed,
                                     NONCE_NATIVE_BINARY_GAMMA_BASE + t);
    }
  }

  compute_linear_constants(ctx);
  compute_linear_fold_constants(ctx);
  compute_native_binary_fold(ctx);
  compute_folded_linear_coeffs(ctx);
  compute_public_fold_rhs(ufold_public,ctx);
  return 0;
}

#if LNP22_DIRECT_CHAIN
#include "src/selector/selector_direct_outer.inc"
#else
static void build_selector_lnp22_outer_statement_public(prncplstmnt *st,
                                                        const selector_lnp22_ctx *ctx,
                                                        const polx *tA, const polx *w,
                                                        const polx *vlin, const polx *tq,
                                                        const polx *vq, const polx *c,
                                                        const polx *a3,
                                                        const polx *ufold_public,
                                                        uint64_t betasq) {
  const size_t r = OUT_R;
  size_t n[OUT_R];
  shake128incctx shakectx;
  sparsecnst *cnst;
#if !OUTER_USE_COMPACT_STATEMENT
  polx *buf;
#endif
  polx *openfold = NULL, *open1 = NULL, *open2 = NULL, *openrhs = NULL;
  polx *quadz = NULL, *quadx = NULL, *quadbits = NULL, *quadz2 = NULL;
  polx *linear1 = NULL, *linear2 = NULL;
  polx linrhs[1], quadrhs[1];
  size_t tail_len, rbit_len, i;
#if OUTER_BZ_DIGITS == 3
  const int64_t bzbase2 = OUTER_BZ_DECOMP_BASE * OUTER_BZ_DECOMP_BASE;
#endif

  tail_len = ctx->rhat_len + ctx->umc_eua_len + ctx->umc_eub_len + ctx->umc_r_len +
             ctx->umc_va_len + ctx->umc_vb_len;
  rbit_len = ctx->rhat_bits_len;

  n[OUT_B0] = ctx->b_len;  n[OUT_B1] = ctx->b_len;
#if OUTER_BZ_DIGITS == 3
  n[OUT_B2] = ctx->b_len;
#endif
  n[OUT_X0] = ctx->x_len;  n[OUT_X1] = ctx->x_len;
  n[OUT_Z0] = ctx->z_len;  n[OUT_Z1] = ctx->z_len;
#if OUTER_BZ_DIGITS == 3
  n[OUT_Z2] = ctx->z_len;
#endif
  n[OUT_Y00] = ctx->y_len; n[OUT_Y01] = ctx->y_len;
  n[OUT_Y10] = ctx->y_len; n[OUT_Y11] = ctx->y_len;
  n[OUT_T0] = tail_len;    n[OUT_T1] = tail_len;
  n[OUT_RBIT0] = rbit_len; n[OUT_RBIT1] = rbit_len;
  n[OUT_Z20] = ctx->m2;    n[OUT_Z21] = ctx->m2;

  init_prncplstmnt_raw(st,r,n,betasq,3,1);
  shake128_inc_init(&shakectx);
#if !OUTER_HASH_EXPANDED_STATEMENT
  absorb_outer_compact_statement(&shakectx,ctx,tA,w,vlin,tq,vq,c);
#endif

  openfold = alloc_polxvec(ctx->kmsis);
  open1 = alloc_polxvec(ctx->m1);
  open2 = alloc_polxvec(ctx->m2);
  openrhs = alloc_polxvec(1);
  quadz = alloc_polxvec(ctx->z_len);
  quadx = alloc_polxvec(ctx->x_len);
  quadbits = alloc_polxvec(rbit_len);
  quadz2 = alloc_polxvec(ctx->m2);
  if(!openfold || !open1 || !open2 || !openrhs || !quadz || !quadx || !quadbits || !quadz2) {
    fprintf(stderr,"ERROR: verifier outer statement allocation failure\n");
    exit(1);
  }

  derive_outer_openfold(openfold,ctx,tA,w,vlin,tq,vq,c);
  compute_openfold_coeffs(open1,open2,openrhs,ctx,tA,w,c,openfold);
  polx_mul(linrhs,c,ufold_public);
  polx_add(linrhs,linrhs,vlin);
  compute_quad_linear_coeffs(quadz,quadx,quadbits,quadz2,ctx,a3,c);
  polx_mul(quadrhs,c,tq);
  polx_add(quadrhs,quadrhs,vq);
#if OUTER_USE_COMPACT_STATEMENT
  /* The verifier context outlives statement verification, so these large
   * transcript-dependent coefficient vectors can be borrowed read-only. */
  linear1 = ctx->lincoeff;
  linear2 = ctx->lincoeff2;
#endif

  cnst = &st->cnst[0];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUT_R,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,ctx->b_len,&open1[ctx->off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,ctx->b_len,&open1[ctx->off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_B2,ctx->b_len,&open1[ctx->off[LNP22_B_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,&open1[ctx->off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,&open1[ctx->off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,ctx->z_len,&open1[ctx->off[LNP22_Z_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Z1,ctx->z_len,&open1[ctx->off[LNP22_Z_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,ctx->z_len,&open1[ctx->off[LNP22_Z_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_Y00,ctx->y_len,&open1[ctx->off[LNP22_Y0_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y01,ctx->y_len,&open1[ctx->off[LNP22_Y0_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Y10,ctx->y_len,&open1[ctx->off[LNP22_Y1_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y11,ctx->y_len,&open1[ctx->off[LNP22_Y1_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,&open1[ctx->off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,&open1[ctx->off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,&open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,&open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,open2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,open2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUT_R,
                             OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
                             OUTER_BZ_DIGITS*ctx->z_len + 2*ctx->y_len + 2*ctx->y_len +
                             2*tail_len + 2*rbit_len + 2*ctx->m2,1,0,0);
  i = 0;
  buf = append_cnst_block(cnst,&i,buf,OUT_B0,ctx->b_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_B1,ctx->b_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_B2,ctx->b_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_X0,ctx->x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_X1,ctx->x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z0,ctx->z_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z1,ctx->z_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_Z2,ctx->z_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_Y00,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y01,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y10,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y11,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T0,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T1,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT0,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT1,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z20,ctx->m2);
  (void)append_cnst_block(cnst,&i,buf,OUT_Z21,ctx->m2);
  set_scaled_polxvec(cnst->phi[OUT_B0],&open1[ctx->off[LNP22_B_BLOCK]],ctx->b_len,1);
  set_scaled_polxvec(cnst->phi[OUT_B1],&open1[ctx->off[LNP22_B_BLOCK]],ctx->b_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_B2],&open1[ctx->off[LNP22_B_BLOCK]],ctx->b_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_X0],&open1[ctx->off[LNP22_X_BLOCK]],ctx->x_len,1);
  set_scaled_polxvec(cnst->phi[OUT_X1],&open1[ctx->off[LNP22_X_BLOCK]],ctx->x_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z0],&open1[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Z1],&open1[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_Z2],&open1[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_Y00],&open1[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y01],&open1[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Y10],&open1[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y11],&open1[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_T0],&open1[ctx->off[LNP22_RHAT_BLOCK]],tail_len,1);
  set_scaled_polxvec(cnst->phi[OUT_T1],&open1[ctx->off[LNP22_RHAT_BLOCK]],tail_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_RBIT0],&open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,1);
  set_scaled_polxvec(cnst->phi[OUT_RBIT1],&open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z20],open2,ctx->m2,1);
  set_scaled_polxvec(cnst->phi[OUT_Z21],open2,ctx->m2,OUTER_Z2_DECOMP_BASE);
#endif
  polxvec_copy(cnst->b,openrhs,1);

  cnst = &st->cnst[1];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUT_R,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,ctx->b_len,&linear1[ctx->off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,ctx->b_len,&linear1[ctx->off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_B2,ctx->b_len,&linear1[ctx->off[LNP22_B_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,&linear1[ctx->off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,&linear1[ctx->off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,ctx->z_len,&linear1[ctx->off[LNP22_Z_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Z1,ctx->z_len,&linear1[ctx->off[LNP22_Z_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,ctx->z_len,&linear1[ctx->off[LNP22_Z_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_Y00,ctx->y_len,&linear1[ctx->off[LNP22_Y0_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y01,ctx->y_len,&linear1[ctx->off[LNP22_Y0_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Y10,ctx->y_len,&linear1[ctx->off[LNP22_Y1_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y11,ctx->y_len,&linear1[ctx->off[LNP22_Y1_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,&linear1[ctx->off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,&linear1[ctx->off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,&linear1[ctx->off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,&linear1[ctx->off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,linear2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,linear2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUT_R,
                             OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
                             OUTER_BZ_DIGITS*ctx->z_len + 2*ctx->y_len + 2*ctx->y_len +
                             2*tail_len + 2*rbit_len + 2*ctx->m2,1,0,0);
  i = 0;
  buf = append_cnst_block(cnst,&i,buf,OUT_B0,ctx->b_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_B1,ctx->b_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_B2,ctx->b_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_X0,ctx->x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_X1,ctx->x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z0,ctx->z_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z1,ctx->z_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_Z2,ctx->z_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_Y00,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y01,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y10,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y11,ctx->y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T0,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T1,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT0,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT1,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z20,ctx->m2);
  (void)append_cnst_block(cnst,&i,buf,OUT_Z21,ctx->m2);
  set_scaled_polxvec(cnst->phi[OUT_B0],&ctx->lincoeff[ctx->off[LNP22_B_BLOCK]],ctx->b_len,1);
  set_scaled_polxvec(cnst->phi[OUT_B1],&ctx->lincoeff[ctx->off[LNP22_B_BLOCK]],ctx->b_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_B2],&ctx->lincoeff[ctx->off[LNP22_B_BLOCK]],ctx->b_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_X0],&ctx->lincoeff[ctx->off[LNP22_X_BLOCK]],ctx->x_len,1);
  set_scaled_polxvec(cnst->phi[OUT_X1],&ctx->lincoeff[ctx->off[LNP22_X_BLOCK]],ctx->x_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z0],&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Z1],&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_Z2],&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_Y00],&ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y01],&ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Y10],&ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y11],&ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_T0],&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK]],tail_len,1);
  set_scaled_polxvec(cnst->phi[OUT_T1],&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK]],tail_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_RBIT0],&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,1);
  set_scaled_polxvec(cnst->phi[OUT_RBIT1],&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z20],ctx->lincoeff2,ctx->m2,1);
  set_scaled_polxvec(cnst->phi[OUT_Z21],ctx->lincoeff2,ctx->m2,OUTER_Z2_DECOMP_BASE);
#endif
  polxvec_copy(cnst->b,linrhs,1);

  cnst = &st->cnst[2];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUTER_BZ_DIGITS + 6,1,1,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,quadx,1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,quadx,OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,ctx->z_len,quadz,1);
  append_cnst_view(cnst,&i,OUT_Z1,ctx->z_len,quadz,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,ctx->z_len,quadz,bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,quadbits,1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,quadbits,OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,quadz2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,quadz2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUTER_BZ_DIGITS + 6,
                             2*ctx->x_len + OUTER_BZ_DIGITS*ctx->z_len + 2*rbit_len +
                             2*ctx->m2,1,1,0);
  {
    size_t pidx = 0;
    const size_t q_x0 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_X0,ctx->x_len);
    const size_t q_x1 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_X1,ctx->x_len);
    const size_t q_z0 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_Z0,ctx->z_len);
    const size_t q_z1 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_Z1,ctx->z_len);
#if OUTER_BZ_DIGITS == 3
    const size_t q_z2 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_Z2,ctx->z_len);
#endif
    const size_t q_rbit0 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_RBIT0,rbit_len);
    const size_t q_rbit1 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_RBIT1,rbit_len);
    const size_t q_z20 = pidx; buf = append_cnst_block(cnst,&pidx,buf,OUT_Z20,ctx->m2);
    const size_t q_z21 = pidx; (void)append_cnst_block(cnst,&pidx,buf,OUT_Z21,ctx->m2);
    set_scaled_polxvec(cnst->phi[q_x0],quadx,ctx->x_len,1);
    set_scaled_polxvec(cnst->phi[q_x1],quadx,ctx->x_len,OUTER_XY_DECOMP_BASE);
    set_scaled_polxvec(cnst->phi[q_z0],quadz,ctx->z_len,1);
    set_scaled_polxvec(cnst->phi[q_z1],quadz,ctx->z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
    set_scaled_polxvec(cnst->phi[q_z2],quadz,ctx->z_len,bzbase2);
#endif
    set_scaled_polxvec(cnst->phi[q_rbit0],quadbits,rbit_len,1);
    set_scaled_polxvec(cnst->phi[q_rbit1],quadbits,rbit_len,OUTER_XY_DECOMP_BASE);
    set_scaled_polxvec(cnst->phi[q_z20],quadz2,ctx->m2,1);
    set_scaled_polxvec(cnst->phi[q_z21],quadz2,ctx->m2,OUTER_Z2_DECOMP_BASE);
  }
#endif
  cnst->a->len = 14;
  {
    size_t idx = 0;
    int xa, yb;
    for(xa=0;xa<2;xa++) {
      for(yb=0;yb<2;yb++) {
        const int64_t scale = -outer_pow_i64(OUTER_XY_DECOMP_BASE,xa + yb);
        cnst->a->rows[idx] = (xa == 0) ? OUT_X0 : OUT_X1;
        cnst->a->cols[idx] = (yb == 0) ? OUT_Y00 : OUT_Y01;
        set_const_coeff(&cnst->a->coeffs[idx],scale);
        idx++;
      }
    }
    for(xa=0;xa<2;xa++) {
      for(yb=0;yb<2;yb++) {
        const int64_t scale = -(int64_t)Y_SPLIT_GAMMA * outer_pow_i64(OUTER_XY_DECOMP_BASE,xa + yb);
        cnst->a->rows[idx] = (xa == 0) ? OUT_X0 : OUT_X1;
        cnst->a->cols[idx] = (yb == 0) ? OUT_Y10 : OUT_Y11;
        set_const_coeff(&cnst->a->coeffs[idx],scale);
        idx++;
      }
    }
    cnst->a->rows[idx] = OUT_X0; cnst->a->cols[idx] = OUT_X0; set_const_coeff(&cnst->a->coeffs[idx],1); idx++;
    cnst->a->rows[idx] = OUT_X0; cnst->a->cols[idx] = OUT_X1; set_const_coeff(&cnst->a->coeffs[idx],2*OUTER_XY_DECOMP_BASE); idx++;
    cnst->a->rows[idx] = OUT_X1; cnst->a->cols[idx] = OUT_X1; set_const_coeff(&cnst->a->coeffs[idx],outer_pow_i64(OUTER_XY_DECOMP_BASE,2)); idx++;
    cnst->a->rows[idx] = OUT_RBIT0; cnst->a->cols[idx] = OUT_RBIT0; set_const_coeff(&cnst->a->coeffs[idx],1); idx++;
    cnst->a->rows[idx] = OUT_RBIT0; cnst->a->cols[idx] = OUT_RBIT1; set_const_coeff(&cnst->a->coeffs[idx],2*OUTER_XY_DECOMP_BASE); idx++;
    cnst->a->rows[idx] = OUT_RBIT1; cnst->a->cols[idx] = OUT_RBIT1; set_const_coeff(&cnst->a->coeffs[idx],outer_pow_i64(OUTER_XY_DECOMP_BASE,2));
  }
  polxvec_copy(cnst->b,quadrhs,1);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);

  free(openfold);
#if OUTER_USE_COMPACT_STATEMENT
  {
    polx *owned[12] = {
      open1,open2,NULL,NULL,quadz,quadx,quadbits,quadz2,
      NULL,NULL,NULL,NULL
    };
    attach_outer_frontend_storage(st,owned);
  }
  open1 = open2 = linear1 = linear2 = NULL;
  quadz = quadx = quadbits = quadz2 = NULL;
#endif
  free(open1);
  free(open2);
  free(openrhs);
  free(quadz);
  free(quadx);
  free(quadbits);
  free(quadz2);
  free(linear1);
  free(linear2);
}

static void prepare_selector_lnp22_outer(prncplstmnt *st, witness *wt, uint64_t seed_offset) {
  outer_instance inst = {};
  struct timespec prepare_start, inner_done, prepare_done;
  const size_t r = OUT_R;
  size_t n[OUT_R];
  uint64_t betasq = 0;
  shake128incctx shakectx;
  sparsecnst *cnst;
#if !OUTER_USE_COMPACT_STATEMENT
  polx *buf;
#endif
  polx *openfold = NULL, *open1 = NULL, *open2 = NULL, *openrhs = NULL;
  polx *quadz = NULL, *quadx = NULL, *quadbits = NULL, *quadz2 = NULL;
  polx *linear1 = NULL, *linear2 = NULL;
  size_t i;
  polx linrhs[1], quadrhs[1];
#if OUTER_BZ_DIGITS == 3
  const int64_t bzbase2 = OUTER_BZ_DECOMP_BASE * OUTER_BZ_DECOMP_BASE;
#endif
  size_t tail_len;
  size_t rbit_len;
#if OUTER_PREP_PROFILE
  struct timespec outer_prof_last, outer_prof_now;
#define OUTER_PROFILE_START() clock_gettime(CLOCK_MONOTONIC,&outer_prof_last)
#define OUTER_PROFILE_MARK(label) do {                                      \
    clock_gettime(CLOCK_MONOTONIC,&outer_prof_now);                         \
    printf("  outer prepare profile: %-24s %.6fs\n",label,                  \
           elapsed_sec(&outer_prof_last,&outer_prof_now));                  \
    outer_prof_last = outer_prof_now;                                        \
  } while(0)
#else
#define OUTER_PROFILE_START() do { } while(0)
#define OUTER_PROFILE_MARK(label) do { } while(0)
#endif

  OUTER_PROFILE_START();
  clock_gettime(CLOCK_MONOTONIC,&prepare_start);
  if(generate_inner_instance(&inst,seed_offset)) {
    fprintf(stderr,"ERROR: failed to generate valid inner LNP22 instance for outer proof\n");
    exit(1);
  }
  clock_gettime(CLOCK_MONOTONIC,&inner_done);
  OUTER_PROFILE_MARK("inner generation");
  printf("Selector outer inner generation: batch=%d, rejection rounds=%zu\n",
         LNP22_BATCH_ROUNDS,inst.rejection_rounds + 1);
  tail_len =
    inst.ctx.rhat_len + inst.ctx.umc_eua_len + inst.ctx.umc_eub_len + inst.ctx.umc_r_len +
    inst.ctx.umc_va_len + inst.ctx.umc_vb_len;
  rbit_len = inst.ctx.rhat_bits_len;

  n[OUT_B0] = inst.ctx.b_len;  n[OUT_B1] = inst.ctx.b_len;
#if OUTER_BZ_DIGITS == 3
  n[OUT_B2] = inst.ctx.b_len;
#endif
  n[OUT_X0] = inst.ctx.x_len;  n[OUT_X1] = inst.ctx.x_len;
  n[OUT_Z0] = inst.ctx.z_len;  n[OUT_Z1] = inst.ctx.z_len;
#if OUTER_BZ_DIGITS == 3
  n[OUT_Z2] = inst.ctx.z_len;
#endif
  n[OUT_Y00] = inst.ctx.y_len; n[OUT_Y01] = inst.ctx.y_len;
  n[OUT_Y10] = inst.ctx.y_len; n[OUT_Y11] = inst.ctx.y_len;
  n[OUT_T0] = tail_len;        n[OUT_T1] = tail_len;
  n[OUT_RBIT0] = rbit_len;     n[OUT_RBIT1] = rbit_len;
  n[OUT_Z20] = inst.ctx.m2;    n[OUT_Z21] = inst.ctx.m2;

  init_witness_raw(wt,r,n);
  OUTER_PROFILE_MARK("dims+alloc");

  if(
#if OUTER_BZ_DIGITS == 2
     decompose_segment_two_digits(wt->s[OUT_B0],wt->s[OUT_B1],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_B_BLOCK]*N],
                                  inst.ctx.b_len,OUTER_BZ_DECOMP_BASE,"B") ||
#else
     decompose_segment_three_digits(wt->s[OUT_B0],wt->s[OUT_B1],wt->s[OUT_B2],
                                    &inst.z1coeffs[inst.ctx.off[LNP22_B_BLOCK]*N],
                                    inst.ctx.b_len,OUTER_BZ_DECOMP_BASE,"B") ||
#endif
     decompose_segment_two_digits(wt->s[OUT_X0],wt->s[OUT_X1],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_X_BLOCK]*N],
                                  inst.ctx.x_len,OUTER_XY_DECOMP_BASE,"X") ||
#if OUTER_BZ_DIGITS == 2
     decompose_segment_two_digits(wt->s[OUT_Z0],wt->s[OUT_Z1],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_Z_BLOCK]*N],
                                  inst.ctx.z_len,OUTER_BZ_DECOMP_BASE,"Z") ||
#else
     decompose_segment_three_digits(wt->s[OUT_Z0],wt->s[OUT_Z1],wt->s[OUT_Z2],
                                    &inst.z1coeffs[inst.ctx.off[LNP22_Z_BLOCK]*N],
                                    inst.ctx.z_len,OUTER_BZ_DECOMP_BASE,"Z") ||
#endif
     decompose_segment_two_digits(wt->s[OUT_Y00],wt->s[OUT_Y01],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_Y0_BLOCK]*N],
                                  inst.ctx.y_len,OUTER_XY_DECOMP_BASE,"Y0") ||
     decompose_segment_two_digits(wt->s[OUT_Y10],wt->s[OUT_Y11],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_Y1_BLOCK]*N],
                                  inst.ctx.y_len,OUTER_XY_DECOMP_BASE,"Y1") ||
     decompose_segment_two_digits(wt->s[OUT_T0],wt->s[OUT_T1],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_RHAT_BLOCK]*N],
                                  tail_len,OUTER_XY_DECOMP_BASE,"tail") ||
     decompose_segment_two_digits(wt->s[OUT_RBIT0],wt->s[OUT_RBIT1],
                                  &inst.z1coeffs[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]*N],
                                  rbit_len,OUTER_XY_DECOMP_BASE,"RBIT") ||
     decompose_segment_two_digits(wt->s[OUT_Z20],wt->s[OUT_Z21],inst.z2coeffs,
                                  inst.ctx.m2,OUTER_Z2_DECOMP_BASE,"z2")) {
    fprintf(stderr,"ERROR: failed to decompose outer witness vectors\n");
    exit(1);
  }
  OUTER_PROFILE_MARK("decompose to witness");

  for(i=0;i<r;i++) {
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    betasq += wt->normsq[i];
  }
  betasq += 16*N;

  init_prncplstmnt_raw(st,r,n,betasq,3,1);
  shake128_inc_init(&shakectx);
#if !OUTER_HASH_EXPANDED_STATEMENT
  absorb_outer_compact_statement(&shakectx,&inst.ctx,inst.tA,inst.w,inst.vlin,inst.tq,inst.vq,inst.c);
#endif

  openfold = alloc_polxvec(inst.ctx.kmsis);
  open1 = alloc_polxvec(inst.ctx.m1);
  open2 = alloc_polxvec(inst.ctx.m2);
  openrhs = alloc_polxvec(1);
  quadz = alloc_polxvec(inst.ctx.z_len);
  quadx = alloc_polxvec(inst.ctx.x_len);
  quadbits = alloc_polxvec(rbit_len);
  quadz2 = alloc_polxvec(inst.ctx.m2);
  if(!openfold || !open1 || !open2 || !openrhs || !quadz || !quadx || !quadbits || !quadz2) {
    fprintf(stderr,"ERROR: outer allocation failure\n");
    exit(1);
  }
  OUTER_PROFILE_MARK("set witness+alloc");

  derive_outer_openfold(openfold,&inst.ctx,inst.tA,inst.w,inst.vlin,inst.tq,inst.vq,inst.c);
  compute_openfold_coeffs(open1,open2,openrhs,&inst.ctx,inst.tA,inst.w,inst.c,openfold);
  polx_mul(linrhs,inst.c,inst.ufold_public);
  polx_add(linrhs,linrhs,inst.vlin);
  compute_quad_linear_coeffs(quadz,quadx,quadbits,quadz2,&inst.ctx,inst.a3,inst.c);
  polx_mul(quadrhs,inst.c,inst.tq);
  polx_add(quadrhs,quadrhs,inst.vq);
#if OUTER_USE_COMPACT_STATEMENT
  linear1 = inst.ctx.lincoeff;
  linear2 = inst.ctx.lincoeff2;
  inst.ctx.lincoeff = NULL;
  inst.ctx.lincoeff2 = NULL;
#endif
  OUTER_PROFILE_MARK("coefficients");

  cnst = &st->cnst[0];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUT_R,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,inst.ctx.b_len,&open1[inst.ctx.off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,inst.ctx.b_len,&open1[inst.ctx.off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_B2,inst.ctx.b_len,&open1[inst.ctx.off[LNP22_B_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_X0,inst.ctx.x_len,&open1[inst.ctx.off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,inst.ctx.x_len,&open1[inst.ctx.off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,inst.ctx.z_len,&open1[inst.ctx.off[LNP22_Z_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Z1,inst.ctx.z_len,&open1[inst.ctx.off[LNP22_Z_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,inst.ctx.z_len,&open1[inst.ctx.off[LNP22_Z_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_Y00,inst.ctx.y_len,&open1[inst.ctx.off[LNP22_Y0_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y01,inst.ctx.y_len,&open1[inst.ctx.off[LNP22_Y0_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Y10,inst.ctx.y_len,&open1[inst.ctx.off[LNP22_Y1_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y11,inst.ctx.y_len,&open1[inst.ctx.off[LNP22_Y1_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,&open1[inst.ctx.off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,&open1[inst.ctx.off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,&open1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,&open1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,inst.ctx.m2,open2,1);
  append_cnst_view(cnst,&i,OUT_Z21,inst.ctx.m2,open2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUT_R,
                             OUTER_BZ_DIGITS*inst.ctx.b_len + 2*inst.ctx.x_len +
                             OUTER_BZ_DIGITS*inst.ctx.z_len +
                             2*inst.ctx.y_len + 2*inst.ctx.y_len + 2*tail_len +
                             2*rbit_len + 2*inst.ctx.m2,
                             1,0,0);
  i = 0;
  buf = append_cnst_block(cnst,&i,buf,OUT_B0,inst.ctx.b_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_B1,inst.ctx.b_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_B2,inst.ctx.b_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_X0,inst.ctx.x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_X1,inst.ctx.x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z0,inst.ctx.z_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z1,inst.ctx.z_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_Z2,inst.ctx.z_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_Y00,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y01,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y10,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y11,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T0,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T1,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT0,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT1,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z20,inst.ctx.m2);
  (void)append_cnst_block(cnst,&i,buf,OUT_Z21,inst.ctx.m2);
  set_scaled_polxvec(cnst->phi[OUT_B0],&open1[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,1);
  set_scaled_polxvec(cnst->phi[OUT_B1],&open1[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_B2],&open1[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_X0],&open1[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len,1);
  set_scaled_polxvec(cnst->phi[OUT_X1],&open1[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z0],&open1[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Z1],&open1[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_Z2],&open1[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_Y00],&open1[inst.ctx.off[LNP22_Y0_BLOCK]],inst.ctx.y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y01],&open1[inst.ctx.off[LNP22_Y0_BLOCK]],inst.ctx.y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Y10],&open1[inst.ctx.off[LNP22_Y1_BLOCK]],inst.ctx.y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y11],&open1[inst.ctx.off[LNP22_Y1_BLOCK]],inst.ctx.y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_T0],&open1[inst.ctx.off[LNP22_RHAT_BLOCK]],tail_len,1);
  set_scaled_polxvec(cnst->phi[OUT_T1],&open1[inst.ctx.off[LNP22_RHAT_BLOCK]],tail_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_RBIT0],&open1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],rbit_len,1);
  set_scaled_polxvec(cnst->phi[OUT_RBIT1],&open1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],rbit_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z20],open2,inst.ctx.m2,1);
  set_scaled_polxvec(cnst->phi[OUT_Z21],open2,inst.ctx.m2,OUTER_Z2_DECOMP_BASE);
#endif
#if OUTER_HASH_EXPANDED_STATEMENT
  for(i=0;i<OUTER_BZ_DIGITS*inst.ctx.b_len + 2*inst.ctx.x_len + OUTER_BZ_DIGITS*inst.ctx.z_len +
             2*inst.ctx.y_len + 2*inst.ctx.y_len + 2*tail_len +
             2*rbit_len + 2*inst.ctx.m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
#endif
  polxvec_copy(cnst->b,openrhs,1);
  OUTER_STMT_ABSORB_POLX(&shakectx,cnst->b);
  OUTER_PROFILE_MARK("constraint opening");

  cnst = &st->cnst[1];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUT_R,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,inst.ctx.b_len,&linear1[inst.ctx.off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,inst.ctx.b_len,&linear1[inst.ctx.off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_B2,inst.ctx.b_len,&linear1[inst.ctx.off[LNP22_B_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_X0,inst.ctx.x_len,&linear1[inst.ctx.off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,inst.ctx.x_len,&linear1[inst.ctx.off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,inst.ctx.z_len,&linear1[inst.ctx.off[LNP22_Z_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Z1,inst.ctx.z_len,&linear1[inst.ctx.off[LNP22_Z_BLOCK]],OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,inst.ctx.z_len,&linear1[inst.ctx.off[LNP22_Z_BLOCK]],bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_Y00,inst.ctx.y_len,&linear1[inst.ctx.off[LNP22_Y0_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y01,inst.ctx.y_len,&linear1[inst.ctx.off[LNP22_Y0_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Y10,inst.ctx.y_len,&linear1[inst.ctx.off[LNP22_Y1_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_Y11,inst.ctx.y_len,&linear1[inst.ctx.off[LNP22_Y1_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,&linear1[inst.ctx.off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,&linear1[inst.ctx.off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,&linear1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,&linear1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,inst.ctx.m2,linear2,1);
  append_cnst_view(cnst,&i,OUT_Z21,inst.ctx.m2,linear2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUT_R,
                             OUTER_BZ_DIGITS*inst.ctx.b_len + 2*inst.ctx.x_len +
                             OUTER_BZ_DIGITS*inst.ctx.z_len +
                             2*inst.ctx.y_len + 2*inst.ctx.y_len + 2*tail_len + 2*rbit_len +
                             2*inst.ctx.m2,
                             1,0,0);
  i = 0;
  buf = append_cnst_block(cnst,&i,buf,OUT_B0,inst.ctx.b_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_B1,inst.ctx.b_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_B2,inst.ctx.b_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_X0,inst.ctx.x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_X1,inst.ctx.x_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z0,inst.ctx.z_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z1,inst.ctx.z_len);
#if OUTER_BZ_DIGITS == 3
  buf = append_cnst_block(cnst,&i,buf,OUT_Z2,inst.ctx.z_len);
#endif
  buf = append_cnst_block(cnst,&i,buf,OUT_Y00,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y01,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y10,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Y11,inst.ctx.y_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T0,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_T1,tail_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT0,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_RBIT1,rbit_len);
  buf = append_cnst_block(cnst,&i,buf,OUT_Z20,inst.ctx.m2);
  (void)append_cnst_block(cnst,&i,buf,OUT_Z21,inst.ctx.m2);
  set_scaled_polxvec(cnst->phi[OUT_B0],&inst.ctx.lincoeff[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,1);
  set_scaled_polxvec(cnst->phi[OUT_B1],&inst.ctx.lincoeff[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_B2],&inst.ctx.lincoeff[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_X0],&inst.ctx.lincoeff[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len,1);
  set_scaled_polxvec(cnst->phi[OUT_X1],&inst.ctx.lincoeff[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z0],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Z1],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  set_scaled_polxvec(cnst->phi[OUT_Z2],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Z_BLOCK]],inst.ctx.z_len,bzbase2);
#endif
  set_scaled_polxvec(cnst->phi[OUT_Y00],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Y0_BLOCK]],inst.ctx.y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y01],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Y0_BLOCK]],inst.ctx.y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Y10],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Y1_BLOCK]],inst.ctx.y_len,1);
  set_scaled_polxvec(cnst->phi[OUT_Y11],&inst.ctx.lincoeff[inst.ctx.off[LNP22_Y1_BLOCK]],inst.ctx.y_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_T0],&inst.ctx.lincoeff[inst.ctx.off[LNP22_RHAT_BLOCK]],tail_len,1);
  set_scaled_polxvec(cnst->phi[OUT_T1],&inst.ctx.lincoeff[inst.ctx.off[LNP22_RHAT_BLOCK]],tail_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_RBIT0],&inst.ctx.lincoeff[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],rbit_len,1);
  set_scaled_polxvec(cnst->phi[OUT_RBIT1],&inst.ctx.lincoeff[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],rbit_len,OUTER_XY_DECOMP_BASE);
  set_scaled_polxvec(cnst->phi[OUT_Z20],inst.ctx.lincoeff2,inst.ctx.m2,1);
  set_scaled_polxvec(cnst->phi[OUT_Z21],inst.ctx.lincoeff2,inst.ctx.m2,OUTER_Z2_DECOMP_BASE);
#endif
#if OUTER_HASH_EXPANDED_STATEMENT
  for(i=0;i<OUTER_BZ_DIGITS*inst.ctx.b_len + 2*inst.ctx.x_len + OUTER_BZ_DIGITS*inst.ctx.z_len +
             2*inst.ctx.y_len + 2*inst.ctx.y_len + 2*tail_len + 2*rbit_len +
             2*inst.ctx.m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
#endif
  polxvec_copy(cnst->b,linrhs,1);
  OUTER_STMT_ABSORB_POLX(&shakectx,cnst->b);
  OUTER_PROFILE_MARK("constraint linear");

  cnst = &st->cnst[2];
#if OUTER_USE_COMPACT_STATEMENT
  init_sparsecnst_view(cnst,r,OUTER_BZ_DIGITS + 6,1,1,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_X0,inst.ctx.x_len,quadx,1);
  append_cnst_view(cnst,&i,OUT_X1,inst.ctx.x_len,quadx,OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z0,inst.ctx.z_len,quadz,1);
  append_cnst_view(cnst,&i,OUT_Z1,inst.ctx.z_len,quadz,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
  append_cnst_view(cnst,&i,OUT_Z2,inst.ctx.z_len,quadz,bzbase2);
#endif
  append_cnst_view(cnst,&i,OUT_RBIT0,rbit_len,quadbits,1);
  append_cnst_view(cnst,&i,OUT_RBIT1,rbit_len,quadbits,OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_Z20,inst.ctx.m2,quadz2,1);
  append_cnst_view(cnst,&i,OUT_Z21,inst.ctx.m2,quadz2,OUTER_Z2_DECOMP_BASE);
#else
  buf = init_sparsecnst_half(cnst,r,OUTER_BZ_DIGITS + 6,
                             2*inst.ctx.x_len + OUTER_BZ_DIGITS*inst.ctx.z_len +
                             2*rbit_len + 2*inst.ctx.m2,1,1,0);
  {
    size_t p = 0;
    const size_t q_x0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_X0,inst.ctx.x_len);
    const size_t q_x1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_X1,inst.ctx.x_len);
    const size_t q_z0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z0,inst.ctx.z_len);
    const size_t q_z1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z1,inst.ctx.z_len);
#if OUTER_BZ_DIGITS == 3
    const size_t q_z2 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z2,inst.ctx.z_len);
#endif
    const size_t q_rbit0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_RBIT0,rbit_len);
    const size_t q_rbit1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_RBIT1,rbit_len);
    const size_t q_z20 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z20,inst.ctx.m2);
    const size_t q_z21 = p; (void)append_cnst_block(cnst,&p,buf,OUT_Z21,inst.ctx.m2);
    set_scaled_polxvec(cnst->phi[q_x0],quadx,inst.ctx.x_len,1);
    set_scaled_polxvec(cnst->phi[q_x1],quadx,inst.ctx.x_len,OUTER_XY_DECOMP_BASE);
    set_scaled_polxvec(cnst->phi[q_z0],quadz,inst.ctx.z_len,1);
    set_scaled_polxvec(cnst->phi[q_z1],quadz,inst.ctx.z_len,OUTER_BZ_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 3
    set_scaled_polxvec(cnst->phi[q_z2],quadz,inst.ctx.z_len,bzbase2);
#endif
    set_scaled_polxvec(cnst->phi[q_rbit0],quadbits,rbit_len,1);
    set_scaled_polxvec(cnst->phi[q_rbit1],quadbits,rbit_len,OUTER_XY_DECOMP_BASE);
    set_scaled_polxvec(cnst->phi[q_z20],quadz2,inst.ctx.m2,1);
    set_scaled_polxvec(cnst->phi[q_z21],quadz2,inst.ctx.m2,OUTER_Z2_DECOMP_BASE);
  }
#endif
#if OUTER_HASH_EXPANDED_STATEMENT
  for(i=0;i<2*inst.ctx.x_len + OUTER_BZ_DIGITS*inst.ctx.z_len + 2*rbit_len + 2*inst.ctx.m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
#endif
  cnst->a->len = 14;
  {
    size_t idx = 0;
    int xa, yb;
    for(xa=0;xa<2;xa++) {
      for(yb=0;yb<2;yb++) {
        const int64_t scale = -outer_pow_i64(OUTER_XY_DECOMP_BASE,xa + yb);
        cnst->a->rows[idx] = (xa == 0) ? OUT_X0 : OUT_X1;
        cnst->a->cols[idx] = (yb == 0) ? OUT_Y00 : OUT_Y01;
        set_const_coeff(&cnst->a->coeffs[idx],scale);
        OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
        idx++;
      }
    }
    for(xa=0;xa<2;xa++) {
      for(yb=0;yb<2;yb++) {
        const int64_t scale = -(int64_t)Y_SPLIT_GAMMA * outer_pow_i64(OUTER_XY_DECOMP_BASE,xa + yb);
        cnst->a->rows[idx] = (xa == 0) ? OUT_X0 : OUT_X1;
        cnst->a->cols[idx] = (yb == 0) ? OUT_Y10 : OUT_Y11;
        set_const_coeff(&cnst->a->coeffs[idx],scale);
        OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
        idx++;
      }
    }
    cnst->a->rows[idx] = OUT_X0;
    cnst->a->cols[idx] = OUT_X0;
    set_const_coeff(&cnst->a->coeffs[idx],1);
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
    idx++;
    cnst->a->rows[idx] = OUT_X0;
    cnst->a->cols[idx] = OUT_X1;
    set_const_coeff(&cnst->a->coeffs[idx],2*OUTER_XY_DECOMP_BASE);
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
    idx++;
    cnst->a->rows[idx] = OUT_X1;
    cnst->a->cols[idx] = OUT_X1;
    set_const_coeff(&cnst->a->coeffs[idx],outer_pow_i64(OUTER_XY_DECOMP_BASE,2));
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
    idx++;
    cnst->a->rows[idx] = OUT_RBIT0;
    cnst->a->cols[idx] = OUT_RBIT0;
    set_const_coeff(&cnst->a->coeffs[idx],1);
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
    idx++;
    cnst->a->rows[idx] = OUT_RBIT0;
    cnst->a->cols[idx] = OUT_RBIT1;
    set_const_coeff(&cnst->a->coeffs[idx],2*OUTER_XY_DECOMP_BASE);
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
    idx++;
    cnst->a->rows[idx] = OUT_RBIT1;
    cnst->a->cols[idx] = OUT_RBIT1;
    set_const_coeff(&cnst->a->coeffs[idx],outer_pow_i64(OUTER_XY_DECOMP_BASE,2));
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->a->coeffs[idx]);
  }
  polxvec_copy(cnst->b,quadrhs,1);
  OUTER_STMT_ABSORB_POLX(&shakectx,cnst->b);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);
  OUTER_PROFILE_MARK("constraint quadratic+hash");

#if OUTER_DEBUG_PRECHECK
  {
    polx *sx[OUT_R] = {0};
    int chk0 = sparsecnst_check(&st->cnst[0],sx,wt);
    int chk1 = sparsecnst_check(&st->cnst[1],sx,wt);
    int chk2 = sparsecnst_check(&st->cnst[2],sx,wt);
    if(!chk0 || !chk1 || !chk2)
      fprintf(stderr,"ERROR: outer precheck failed: c0=%d c1=%d c2=%d\n",chk0,chk1,chk2);
    free(sx[0]);
  }
  OUTER_PROFILE_MARK("debug precheck");
#endif

  free(openfold);
#if OUTER_USE_COMPACT_STATEMENT
  {
    polx *owned[12] = {
      open1,open2,linear1,linear2,quadz,quadx,quadbits,quadz2,
      NULL,NULL,NULL,NULL
    };
    attach_outer_frontend_storage(st,owned);
  }
  open1 = open2 = linear1 = linear2 = NULL;
  quadz = quadx = quadbits = quadz2 = NULL;
#endif
  free(open1);
  free(open2);
  free(openrhs);
  free(quadz);
  free(quadx);
  free(quadbits);
  free(quadz2);
  free(linear1);
  free(linear2);
  free_outer_instance(&inst);
  clock_gettime(CLOCK_MONOTONIC,&prepare_done);
  printf("Selector outer prepare breakdown: inner-generation=%.4fs, outer-frontend=%.4fs\n",
         elapsed_sec(&prepare_start,&inner_done),elapsed_sec(&inner_done,&prepare_done));
  OUTER_PROFILE_MARK("cleanup");
#undef OUTER_PROFILE_START
#undef OUTER_PROFILE_MARK
}
#endif

static int test_outer_twolayer(uint64_t seed_offset) {
  int ret;
  prncplstmnt st0 = {};
  statement st1 = {};
  proof pi0 = {};
  witness wt0 = {}, wt1 = {};
  double size = 0.0;

  printf("Testing selector LNP22 outer Labrador compression\n\n");
  prepare_selector_lnp22_outer(&st0,&wt0,seed_offset);
  print_prncplstmnt_pp(&st0);
  ret = principle_verify(&st0,&wt0);
  if(ret) {
    fprintf(stderr,"ERROR: verification of outer instance failed: %d\n",ret);
    goto end;
  }

  ret = principle_prove(&st1,&wt1,&pi0,&st0,&wt0,0);
  if(ret) {
    fprintf(stderr,"ERROR: outer Chihuahua proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt0);
  size += print_proof_pp(&pi0);
  print_statement_pp(&st1);
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: outer Chihuahua verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st1);
  ret = principle_reduce(&st1,&pi0,&st0);
  free_prncplstmnt(&st0);
  if(ret) {
    fprintf(stderr,"ERROR: outer Labrador reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: verification of outer reduced statement failed: %d\n",ret);
    goto end;
  }

  size += print_witness_pp(&wt1);
  printf("Reduced outer proof size: %.2f KB\n\n",size);

end:
  free_prncplstmnt(&st0);
  free_statement(&st1);
  free_proof(&pi0);
  free_witness(&wt0);
  free_witness(&wt1);
  return ret;
}

static int test_outer_pack(uint64_t seed_offset) {
  int ret;
  prncplstmnt st = {};
  witness wt = {};
  composite p = {};
  outer_instance pubinst = {};
  struct timespec prep0, prep1, prove0, prove1, verify0, verify1;

  printf("Testing selector LNP22 outer composite\n\n");
  clock_gettime(CLOCK_MONOTONIC,&prep0);
#if LNP22_DIRECT_CHAIN
  prepare_selector_lnp22_outer_direct(&st,&wt,seed_offset,&pubinst);
#else
  prepare_selector_lnp22_outer(&st,&wt,seed_offset);
#endif
  clock_gettime(CLOCK_MONOTONIC,&prep1);
  print_prncplstmnt_pp(&st);
  ret = principle_verify(&st,&wt);
  if(ret) {
    fprintf(stderr,"ERROR: verification of outer composite instance failed: %d\n",ret);
    goto end;
  }

  clock_gettime(CLOCK_MONOTONIC,&prove0);
  ret = composite_prove_principle(&p,&st,&wt);
  clock_gettime(CLOCK_MONOTONIC,&prove1);
  if(ret) {
    fprintf(stderr,"ERROR: outer composite proof failed: %d\n",ret);
    goto end;
  }
  clock_gettime(CLOCK_MONOTONIC,&verify0);
  ret = composite_verify_principle(&p,&st);
  clock_gettime(CLOCK_MONOTONIC,&verify1);
  if(ret) {
    fprintf(stderr,"ERROR: outer composite verification failed: %d\n",ret);
    goto end;
  }
  printf("Selector outer composite wall time: prepare=%.4fs, prove=%.4fs, verify=%.4fs, prover-total=%.4fs\n",
         elapsed_sec(&prep0,&prep1),elapsed_sec(&prove0,&prove1),
         elapsed_sec(&verify0,&verify1),
         elapsed_sec(&prep0,&prep1) + elapsed_sec(&prove0,&prove1));

#if OUTER_VERIFY_REQUEST_BENCH
  {
    selector_lnp22_ctx vctx = {};
    prncplstmnt vst = {};
    polx vufold_public[1];
    uint8_t seeded_demo_seed[16];
    struct timespec req0, req1, base0, base1, ctx0, ctx1, stmt0, stmt1, vfy0, vfy1;
    double public_base_sec, public_ctx_sec, stmt_sec, vfy_sec;

#if !LNP22_DIRECT_CHAIN
    if(generate_inner_instance(&pubinst,seed_offset)) {
      fprintf(stderr,"ERROR: failed to regenerate public transcript for verify-request benchmark\n");
      ret = 1;
      goto end;
    }
#endif

    add_seed_offset(seeded_demo_seed,demo_fixed_seed,seed_offset);
    clock_gettime(CLOCK_MONOTONIC,&base0);
    selector_ctx_init(&vctx,seeded_demo_seed);
    prepare_selector_lnp22_verifier_static(&vctx);
    clock_gettime(CLOCK_MONOTONIC,&base1);

    clock_gettime(CLOCK_MONOTONIC,&req0);
    clock_gettime(CLOCK_MONOTONIC,&ctx0);
    if(prepare_selector_lnp22_verifier_request(&vctx,&pubinst.ctx,
                                               vufold_public,pubinst.tA)) {
      fprintf(stderr,"ERROR: invalid public coefficient-constantness transcript\n");
      ret = 1;
      selector_ctx_clear(&vctx);
      goto end;
    }
    clock_gettime(CLOCK_MONOTONIC,&ctx1);

    clock_gettime(CLOCK_MONOTONIC,&stmt0);
    build_selector_lnp22_outer_statement_public(&vst,&vctx,
                                                pubinst.tA,pubinst.w,pubinst.vlin,
                                                pubinst.tq,pubinst.vq,pubinst.c,
                                                pubinst.a3,vufold_public,st.betasq);
    clock_gettime(CLOCK_MONOTONIC,&stmt1);
    if(memcmp(vst.h,st.h,sizeof(st.h)) != 0) {
      fprintf(stderr,"WARNING: verifier rebuilt outer statement hash differs from prover statement hash\n");
    }
#if defined(DIRECT_OUTER_DEBUG) && DIRECT_OUTER_DEBUG
    {
      const outer_frontend_storage *pst = st.frontend_data;
      int raw_ret = principle_verify(&vst,&wt);
      fprintf(stderr,
              "DEBUG public coeffs: fold=%d tg=%d h=%d rho=%d A0=%d Dsel=%d directD=%d\n",
              polxvec_equal_slow(pubinst.ctx.directfold,vctx.directfold,
                                 vctx.Lchain-1),
              polxvec_equal_slow(pubinst.ctx.rbit_const_tg,
                                 vctx.rbit_const_tg,vctx.rbit_const_lambda),
              polxvec_equal_slow(pubinst.ctx.rbit_const_h,
                                 vctx.rbit_const_h,vctx.rbit_const_lambda),
              polxvec_equal_slow(pubinst.ctx.rhovec,vctx.rhovec,
                                 (vctx.Lchain-1)*vctx.m),
              polxvec_equal_slow(pubinst.ctx.A0,vctx.A0,vctx.m*vctx.gdigits),
              polxvec_equal_slow(pubinst.ctx.Dsel,vctx.Dsel,
                                 vctx.m*vctx.gdigits),
              polxvec_equal_slow(pst->owned[9],vctx.directD,
                                 (vctx.Lchain-1)*vctx.gdigits));
      fprintf(stderr,"DEBUG direct outer: rebuilt raw statement=%s (%d)\n",
              raw_ret ? "mismatch" : "ok",raw_ret);
    }
#endif

    clock_gettime(CLOCK_MONOTONIC,&vfy0);
    ret = composite_verify_principle(&p,&vst);
    clock_gettime(CLOCK_MONOTONIC,&vfy1);
    clock_gettime(CLOCK_MONOTONIC,&req1);
    if(ret) {
      fprintf(stderr,"ERROR: verify-request Labrador verification failed: %d\n",ret);
      free_prncplstmnt(&vst);
      selector_ctx_clear(&vctx);
#if !LNP22_DIRECT_CHAIN
      free_outer_instance(&pubinst);
#endif
      goto end;
    }

    public_base_sec = elapsed_sec(&base0,&base1);
    public_ctx_sec = elapsed_sec(&ctx0,&ctx1);
    stmt_sec = elapsed_sec(&stmt0,&stmt1);
    vfy_sec = elapsed_sec(&vfy0,&vfy1);
    printf("Selector verify-request wall time: public-base(one-time)=%.4fs, public-request=%.4fs, fs+ar/deser=%.4fs, outer-statement=%.4fs, outer-verify=%.4fs, total=%.4fs\n",
           public_base_sec,public_ctx_sec,0.0,stmt_sec,vfy_sec,elapsed_sec(&req0,&req1));

    free_prncplstmnt(&vst);
    selector_ctx_clear(&vctx);
#if !LNP22_DIRECT_CHAIN
    free_outer_instance(&pubinst);
#endif
  }
#endif

end:
#if LNP22_DIRECT_CHAIN
  free_outer_instance(&pubinst);
#endif
  free_prncplstmnt(&st);
  free_witness(&wt);
  free_composite(&p);
  return ret;
}

int main(int argc, char **argv) {
  int ret;
  uint64_t seed_offset = LNP22_DEFAULT_SEED_OFFSET;

  configure_openmp_threads();

  if(argc > 2) {
    fprintf(stderr,"usage: %s [seed_offset]\n",argv[0]);
    return 2;
  }
  if(argc == 2) {
    char *endptr = NULL;
    unsigned long long parsed_offset = strtoull(argv[1],&endptr,10);

    if(endptr == NULL || *endptr != '\0') {
      fprintf(stderr,"invalid seed_offset: %s\n",argv[1]);
      return 2;
    }
    seed_offset = (uint64_t)parsed_offset;
  }

#if OUTER_RUN_TWOLAYER
  ret = test_outer_twolayer(seed_offset);
  if(ret) goto end;
#endif
  ret = test_outer_pack(seed_offset);
  if(ret) goto end;

end:
  free_comkey();
  return ret;
}
