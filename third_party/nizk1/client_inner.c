#define _POSIX_C_SOURCE 200809L
#include "src/core/aesctr.h"
#ifndef DEMO_L
#define DEMO_L 80
#endif
#ifndef LOGQ
#define LOGQ 48
#endif
#ifndef LEOPARD_PROOF_FREE
#define LEOPARD_PROOF_FREE 0
#endif
#ifndef DEMO_M
#if LOGQ == 72
#define DEMO_M 42
#elif LOGQ >= 64
#define DEMO_M 37
#elif LOGQ >= 56
#define DEMO_M 34
#else
#define DEMO_M 24
#endif
#endif
#ifndef GADGET_BASE
#if LOGQ >= 64
#define GADGET_BASE 319
#else
#define GADGET_BASE 255
#endif
#endif
#ifndef Y_SPLIT_GAMMA
#define Y_SPLIT_GAMMA 200
#endif
#ifndef LNP22_KMSIS
#if LOGQ == 72
#define LNP22_KMSIS 24
#elif LOGQ >= 64
#define LNP22_KMSIS 23
#elif LOGQ >= 56
#define LNP22_KMSIS 25
#else
#define LNP22_KMSIS 30
#endif
#endif
#ifndef LNP22_M2
#if LOGQ == 72
#define LNP22_M2 74
#elif LOGQ >= 64
#define LNP22_M2 64
#elif LOGQ >= 56
#define LNP22_M2 62
#else
#define LNP22_M2 62
#endif
#endif
#ifndef UMC_ELL
#if LOGQ == 72
#define UMC_ELL 45
#elif LOGQ >= 64
#define UMC_ELL 40
#elif LOGQ >= 56
#define UMC_ELL 38
#else
#define UMC_ELL 27
#endif
#endif
#ifndef LNP22_OMEGA
#define LNP22_OMEGA 8
#endif
#ifndef LNP22_LOG2STDEV1
#if LOGQ >= 64
#define LNP22_LOG2STDEV1 27
#elif LOGQ >= 56
#define LNP22_LOG2STDEV1 26
#else
#define LNP22_LOG2STDEV1 26
#endif
#endif
#ifndef LNP22_LOG2STDEV2
#define LNP22_LOG2STDEV2 17
#endif
#ifndef LNP22_GAMMA1
#define LNP22_GAMMA1 10.0
#endif
#ifndef LNP22_GAMMA2
#define LNP22_GAMMA2 20.0
#endif
#ifndef LNP22_BATCH_ROUNDS
#define LNP22_BATCH_ROUNDS 1
#endif
#ifndef LNP22_ENABLE_REJECTION
#define LNP22_ENABLE_REJECTION 1
#endif
#ifndef LNP22_DEFAULT_SEED_OFFSET
#if LOGQ == 72
#define LNP22_DEFAULT_SEED_OFFSET 2
#elif LOGQ >= 64
#define LNP22_DEFAULT_SEED_OFFSET 11
#elif LOGQ >= 56
#define LNP22_DEFAULT_SEED_OFFSET 2
#else
#define LNP22_DEFAULT_SEED_OFFSET 11
#endif
#endif
#ifndef LNP22_RBIT_CONSTCHECK_LAMBDA
#define LNP22_RBIT_CONSTCHECK_LAMBDA (LIFTS+2)
#endif
#define RBIT_ARP_ROWS 2
#define RBIT_CONST_ROWS LIFTS
#define RBIT_ARP_LOGSIGMA 20
#define RBIT_ARP_LIMIT ((int64_t)1<<25)
#define RBIT_ARP_SIGMA (1.55*(1ULL<<RBIT_ARP_LOGSIGMA))
#ifndef LNP22_PROVE_X_BITS
#define LNP22_PROVE_X_BITS 1
#endif
#ifndef LNP22_CACHE_A1OPEN
#define LNP22_CACHE_A1OPEN 0
#endif
#ifndef LNP22_VALIDATE_BASE_WITNESS
#define LNP22_VALIDATE_BASE_WITNESS 0
#endif
#ifndef LNP22_COMBINE_TA_W
#define LNP22_COMBINE_TA_W 1
#endif
#ifndef LNP22_LINEAR_PROFILE
#define LNP22_LINEAR_PROFILE 0
#endif
#ifndef LNP22_ROW_PARALLEL
#define LNP22_ROW_PARALLEL 1
#endif
#ifndef LNP22_PREP_ROW_PARALLEL
#define LNP22_PREP_ROW_PARALLEL 1
#endif
#ifndef LNP22_DIRECT_COEFF_PARALLEL
#define LNP22_DIRECT_COEFF_PARALLEL 1
#endif
#ifndef LNP22_DEFAULT_OMP_THREADS
#define LNP22_DEFAULT_OMP_THREADS 8
#endif
#ifndef LNP22_PARALLEL_MIN_LEN
#define LNP22_PARALLEL_MIN_LEN 4096
#endif
#ifndef LNP22_COEFF_RESPONSE
#define LNP22_COEFF_RESPONSE 1
#endif
#ifndef LNP22_DIRECT_CHAIN
#define LNP22_DIRECT_CHAIN 1
#endif
#ifndef LNP22_SOUNDNESS_SELFTEST
#define LNP22_SOUNDNESS_SELFTEST 0
#endif

#define LNP22_RHAT_BINARY_BACKEND_RBIT 1
#define LNP22_RHAT_BINARY_BACKEND_NATIVE_FOLD 2
#ifndef LNP22_RHAT_BINARY_BACKEND
#define LNP22_RHAT_BINARY_BACKEND LNP22_RHAT_BINARY_BACKEND_RBIT
#endif
#if LNP22_RHAT_BINARY_BACKEND != LNP22_RHAT_BINARY_BACKEND_RBIT && \
    LNP22_RHAT_BINARY_BACKEND != LNP22_RHAT_BINARY_BACKEND_NATIVE_FOLD
#error "unsupported LNP22_RHAT_BINARY_BACKEND"
#endif
#define LNP22_USE_RBIT_BINARY \
  (LNP22_RHAT_BINARY_BACKEND == LNP22_RHAT_BINARY_BACKEND_RBIT)
#define LNP22_USE_NATIVE_FOLD_BINARY \
  (LNP22_RHAT_BINARY_BACKEND == LNP22_RHAT_BINARY_BACKEND_NATIVE_FOLD)
#if LNP22_USE_RBIT_BINARY && !LNP22_DIRECT_CHAIN
#error "Weighted RBIT requires the direct-chain inner/outer frontend"
#endif
#if !LNP22_USE_RBIT_BINARY || !LNP22_PROVE_X_BITS
#error "The binary-opening experiment requires scalar-bit and constantness proofs"
#endif
#define LNP22_ENABLE_RBIT_CONSTCHECK \
  (LNP22_RBIT_CONSTCHECK_LAMBDA > 0 && \
   (LNP22_USE_RBIT_BINARY || LNP22_USE_NATIVE_FOLD_BINARY || LNP22_PROVE_X_BITS))

#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <immintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "malloc.h"
#include "gaussian_lnp.h"

#define SELECTOR_NO_MAIN
#define SELECTOR_SINGLEPROOF_NO_MAIN
#include "src/selector/selector_singleproof_support.c"

enum {
  LNP22_B_BLOCK = 0,
  LNP22_X_BLOCK = 1,
  LNP22_Z_BLOCK = 2,
  LNP22_Y0_BLOCK = 3,
  LNP22_Y1_BLOCK = 4,
  LNP22_RHAT_BLOCK = 5,
  LNP22_UMC_EUA_BLOCK = 6,
  LNP22_UMC_EUB_BLOCK = 7,
  LNP22_UMC_R_BLOCK = 8,
  LNP22_UMC_VA_BLOCK = 9,
  LNP22_UMC_VB_BLOCK = 10,
  LNP22_RHAT_BITS_BLOCK = 11,
  LNP22_NUM_BLOCKS = 12
};

enum {
  NONCE_UMC_A_BASE = 0x01000000,
  NONCE_UMC_B_BASE = 0x02000000,
  NONCE_UMC_EUA_BASE = 1200,
  NONCE_UMC_EUB_BASE = 1300,
  NONCE_UMC_R_BASE = 1400,
  NONCE_RHAT_BASE = 1500,
  NONCE_AR_BASE = 1600,
  NONCE_LNP_A1_BASE = 0x10000000,
  NONCE_LNP_A2_BASE = 0x20000000,
  NONCE_LNP_A3_BASE = 0x30000000,
  NONCE_LINFOLD_BASE = 5000,
  NONCE_S2_BASE = 6000,
  NONCE_Y1_BASE = 7000,
  NONCE_Y2_BASE = 8000,
  NONCE_RBIT_CONST_BG_BASE = 9000,
  NONCE_RBIT_CONST_G_BASE = 10000,
  NONCE_RBIT_CONST_PHI_BASE = 11000,
  NONCE_X_CONST_PHI_BASE = 12000,
  NONCE_NATIVE_BINARY_GAMMA_BASE = 13000
};

#ifndef UMC_MODP
#if LEOPARD_PROOF_FREE
#define UMC_MODP 769
#else
#define UMC_MODP 65437
#endif
#endif

#ifndef UMC_MSG_SCALE
#define UMC_MSG_SCALE (UMC_MODP/2)
#endif

#ifndef UMC_N
#if LEOPARD_PROOF_FREE
#define UMC_N 7
#else
#define UMC_N 14
#endif
#endif

#ifndef UMC_K
#if LEOPARD_PROOF_FREE
#define UMC_K 7
#else
#define UMC_K 14
#endif
#endif

#if !LEOPARD_PROOF_FREE
_Static_assert(LOGQ == 72 && UMC_MODP == 65437 && UMC_N == 14 && UMC_K == 14,
               "compact lifting requires its audited q72 commitment parameters");
#else
_Static_assert(LOGQ == 56 || LOGQ == 64, "unsupported proof-free profile");
#endif
_Static_assert(LNP22_RBIT_CONSTCHECK_LAMBDA > 0,
               "binary opening requires masked constantness checks");

typedef struct {
  size_t Lchain;
  size_t preproc_count;
  size_t selected_j;
  size_t m;
  size_t gdigits;
  size_t b_len;
  size_t x_len;
  size_t z_len;
  size_t y_len;
  size_t rhat_len;
  size_t rhat_bits_len;
  size_t umc_n;
  size_t umc_k;
  size_t umc_v;
  size_t umc_eua_len;
  size_t umc_eub_len;
  size_t umc_r_len;
  size_t umc_va_len;
  size_t umc_vb_len;
  size_t off[LNP22_NUM_BLOCKS];
  size_t n[LNP22_NUM_BLOCKS];
  size_t m1;
  size_t m2;
  size_t m2prime;
  size_t kmsis;
  size_t cx_rows;
  size_t nlin;
  size_t rbit_const_lambda;
  int cache_a1open;
  uint8_t prechallenge_seed[16];

  poly *witness[LNP22_NUM_BLOCKS];
  polx *s1;

  polx *A0;
  polx *A1sel;
  polx *Dsel;
  polx *Aumc;
  polx *Bumc;
  polx *CX;
  polx *Arhash;
  polx *Cout;
  polx *A1open;
  polx *A2open;
  polx *kvec;
  polx *kprime;
  polx *rhovec;
  polx *alphavec;
  polx *directfold;
  polx *xbit_weight;
  polx *rbit_weight;
  polx *directD;
  polx *direct_lin_B;
  polx *g0;
  polx *g1;
  polx *gdelta;
  polx *linfold;
  polx *ulin;
  polx *lincoeff;
  polx *lincoeff2;
  polx *rbit_const_Bg;
  polx *rbit_const_g;
  polx *rbit_const_tg;
  polx *rbit_const_h;
  polx *native_binary_gamma;
  polx *native_binary_scale;
  polx *native_binary_lin1;
  polx *native_binary_lin2;
  polx *native_binary_uquad;
  polx *g1coeff;

  poly *pre_rhat;
  poly *pre_eua;
  poly *pre_eub;
  poly *pre_r;
  poly *pre_va;
  poly *pre_vb;
  polx *pre_CX;
  uint8_t *pre_ar_seed;

  int64_t umc_va_min;
  int64_t umc_va_max;
  uint64_t umc_va_maxabs;
  int64_t umc_vb_min;
  int64_t umc_vb_max;
  uint64_t umc_vb_maxabs;

  uint8_t statement_seed[16];
  uint8_t selector_chal_seed[16];
  uint8_t rbit_const_seed[16];
  uint8_t lnp_seed[16];
} selector_lnp22_ctx;

typedef struct {
  polx *y1;
  polx *y2;
  polx *w;
#if !LNP22_COEFF_RESPONSE
  polx *cs1;
  polx *cs2;
#endif
  int64_t *z1coeffs;
  int64_t *z2coeffs;
  int64_t *cs1coeffs;
  int64_t *cs2coeffs;
  int64_t ccoeffs[N];
  polx *c;
  polx *vlin;
  polx *tq;
  polx *vq;
  polx *g1;
  polx *g0;
  int rej1;
  int rej2;
} lnp22_round_ws;

static void quad_value(polx *out, const selector_lnp22_ctx *ctx, const polx *vec);
static void polxvec_centered_to_i64(int64_t *out, const polx *v, size_t len);
static void rhat_native_binary_value(polx *out, const selector_lnp22_ctx *ctx,
                                     const polx *rhat);

static size_t round_up_size(size_t n, size_t a) {
  return ((n + a - 1)/a)*a;
}

static double elapsed_sec(const struct timespec *start, const struct timespec *end) {
  return (double)(end->tv_sec - start->tv_sec) +
         1e-9 * (double)(end->tv_nsec - start->tv_nsec);
}

static size_t selector_x_pos(const selector_lnp22_ctx *ctx, size_t t) {
  (void)ctx;
#if LNP22_DIRECT_CHAIN
  return t;
#else
  return t*SMALL_BLOCK;
#endif
}

static int64_t zz_to_i64(const zz *a) {
  const __int128 r = zz_toint128(a);
  if(r < INT64_MIN || r > INT64_MAX) {
    fprintf(stderr,"centered coefficient does not fit a short int64 buffer\n");
    abort();
  }
  return (int64_t)r;
}

static polx *alloc_polxvec(size_t len) {
  if(len == 0)
    return NULL;
  return _aligned_alloc(64,round_up_size(len*sizeof(polx),64));
}

static poly *alloc_polyvec(size_t len) {
  if(len == 0)
    return NULL;
  return _aligned_alloc(64,round_up_size(len*sizeof(poly),64));
}

static void store_u64_le(uint8_t out[8], uint64_t x) {
  size_t i;
  for(i=0;i<8;i++)
    out[i] = (uint8_t)(x >> (8*i));
}

static uint64_t load_u64_le(const uint8_t in[8]) {
  uint64_t x = 0;
  size_t i;

  for(i=0;i<8;i++)
    x |= (uint64_t)in[i] << (8*i);
  return x;
}

static double sigma_from_log2stdev(unsigned int log2stdev) {
  return 1.55 * ldexp(1.0,(int)log2stdev);
}

static double rejection_logM(double gamma) {
  return sqrt(2.0 * (128.0 + 1.0) / log2(exp(1.0))) / gamma + 0.5 / (gamma * gamma);
}

static double rejection_logM_rej2(double gamma) {
  return 0.5 / (gamma * gamma);
}

static void sample_small_poly(poly *p, const uint8_t seed[32], uint64_t nonce, int bound) {
  uint8_t in[40];
  uint8_t out[2*N];
  size_t i;

  memcpy(in,seed,32);
  store_u64_le(&in[32],nonce);
  shake128(out,sizeof(out),in,sizeof(in));
  poly_set_zero(p);
  for(i=0;i<N;i++)
    p->vec->c[i] = (int16_t)((out[i] % (2*bound + 1)) - bound);
}

static void sample_small_polx(polx *r, const uint8_t seed[32], uint64_t nonce, int bound) {
  poly tmp[1];
  sample_small_poly(tmp,seed,nonce,bound);
  polx_frompoly(r,tmp);
}

static int64_t centered_div_i64(int64_t x, int64_t mod, int64_t *rem) {
  int64_t q;

  if(x >= 0)
    q = (x + mod/2) / mod;
  else
    q = -(((-x) + mod/2) / mod);
  *rem = x - mod*q;
  return q;
}

static uint64_t shake_squeeze_u64(shake128incctx *st);

static uint64_t shake_uniform_mod(shake128incctx *st, uint64_t mod) {
  const uint64_t limit = UINT64_MAX - (UINT64_MAX % mod);
  uint64_t x;

  do {
    x = shake_squeeze_u64(st);
  } while(x >= limit);
  return x % mod;
}

static void sample_autostable_polx(polx *r, const uint8_t seed[32], uint64_t nonce, int bound) {
  shake128incctx st;
  uint8_t in[40];
  poly tmp[1];
  const size_t nchal = N/2;
  const uint64_t mod = (uint64_t)(2*bound + 1);
  int64_t chal[N/2];
  size_t i;

  memcpy(in,seed,32);
  store_u64_le(&in[32],nonce);
  shake128_inc_init(&st);
  shake128_inc_absorb(&st,in,sizeof(in));
  shake128_inc_finalize(&st);

  /* Auto-stable challenges have real eigenvalues at the odd 128th roots.
     Keep an explicit numerical margin below the parameter-screen eta. */
  for(;;) {
    int accept=1;
    for(i=0;i<nchal;i++)
      chal[i]=(int64_t)shake_uniform_mod(&st,mod)-bound;
    for(size_t root=0;root<N/2;root++) {
      long double value=chal[0];
      const long double theta=acosl(-1.0L)*(2*root+1)/N;
      for(i=1;i<nchal;i++) value+=2*chal[i]*cosl(theta*i);
      if(fabsl(value)>127.0L-0.000001L) { accept=0; break; }
    }
    if(accept) break;
  }
  poly_set_zero(tmp);
  for(i=0;i<nchal;i++) tmp->vec->c[i]=(int16_t)chal[i];
  tmp->vec->c[nchal] = 0;
  for(i=nchal+1;i<N;i++)
    tmp->vec->c[i] = (int16_t)(-chal[2*nchal - i]);

  polx_frompoly(r,tmp);
}

static void sample_scalar_challenge_polxvec(polx *r, size_t len,
                                            const uint8_t seed[16],
                                            uint64_t nonce) {
  uint8_t in[24];
  uint8_t out[8*32];
  size_t i, off = sizeof(out);

  memcpy(in,seed,16);
  store_u64_le(&in[16],nonce);
  for(i=0;i<len;i++) {
    int64_t v;

    if(off + 8 > sizeof(out)) {
      store_u64_le(&in[16],nonce++);
      shake128(out,sizeof(out),in,sizeof(in));
      off = 0;
    }
    v = (int64_t)(load_u64_le(&out[off]) & ((((uint64_t)1) << 47) - 1));
    off += 8;
    v -= ((int64_t)1) << 46;
    polx_monomial(&r[i],v,0);
  }
}

static void sample_full_scalar_challenge(polx *r, const uint8_t seed[16],
                                         uint64_t nonce) {
  polx uniform[1];
  int64_t coeffs[N];

  polxvec_almostuniform(uniform,1,seed,nonce);
  polxvec_centered_to_i64(coeffs,uniform,1);
  polx_monomial(r,coeffs[0],0);
}

static void sample_cbd1_poly(poly *p, const uint8_t seed[16], uint64_t nonce) {
  uint8_t in[24];
  uint8_t out[16];
  size_t i;

  memcpy(in,seed,16);
  store_u64_le(&in[16],nonce);
  shake128(out,sizeof(out),in,sizeof(in));
  poly_set_zero(p);
  for(i=0;i<N;i++) {
    const uint8_t bits = out[i >> 2];
    const int a = (bits >> (2*(i & 3))) & 1;
    const int b = (bits >> (2*(i & 3) + 1)) & 1;
    p->vec->c[i] = (int16_t)(a - b);
  }
}

static void sample_cbd1_polxvec(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce_base) {
  poly tmp[1];
  size_t i;

  for(i=0;i<len;i++) {
    sample_cbd1_poly(tmp,seed,nonce_base + i);
    polx_frompoly(&r[i],tmp);
  }
}

static void sample_uniform_modp_poly(poly *p, const uint8_t seed[16], uint64_t nonce, int mod) {
  shake128incctx st;
  uint8_t in[24];
  size_t i;

  memcpy(in,seed,16);
  store_u64_le(&in[16],nonce);
  shake128_inc_init(&st);
  shake128_inc_absorb(&st,in,sizeof(in));
  shake128_inc_finalize(&st);
  poly_set_zero(p);
  for(i=0;i<N;i++) {
    int64_t x = (int64_t)shake_uniform_mod(&st,(uint64_t)mod);
    if(x > mod/2)
      x -= mod;
    p->vec->c[i] = (int16_t)x;
  }
}

static void sample_uniform_modp_polxvec(polx *r, size_t len, const uint8_t seed[16],
                                        uint64_t nonce_base, int mod) {
  poly tmp[1];
  size_t i;

  for(i=0;i<len;i++) {
    sample_uniform_modp_poly(tmp,seed,nonce_base + i,mod);
    polx_frompoly(&r[i],tmp);
  }
}

static void sample_binary_poly(poly *p, const uint8_t seed[16], uint64_t nonce) {
  uint8_t in[24];
  uint8_t out[8];
  size_t i;

  memcpy(in,seed,16);
  store_u64_le(&in[16],nonce);
  shake128(out,sizeof(out),in,sizeof(in));
  poly_set_zero(p);
  for(i=0;i<N;i++)
    p->vec->c[i] = (int16_t)((out[i >> 3] >> (i & 7)) & 1u);
}

static uint64_t shake_squeeze_u64(shake128incctx *st) {
  uint8_t out[8];
  shake128_inc_squeeze(out,sizeof(out),st);
  return ((uint64_t)out[0]) |
         ((uint64_t)out[1] << 8) |
         ((uint64_t)out[2] << 16) |
         ((uint64_t)out[3] << 24) |
         ((uint64_t)out[4] << 32) |
         ((uint64_t)out[5] << 40) |
         ((uint64_t)out[6] << 48) |
         ((uint64_t)out[7] << 56);
}

static void sample_zero_const_polx(polx *r, const uint8_t seed[16], uint64_t nonce) {
  polz uniform[1];

  /* The mask must cover the full zero-constant subspace, including q > 2^64. */
  polzvec_uniform(uniform,1,seed,nonce);
  polz_setcoeff_fromint64(uniform,0,0);
  polz_topolx(r,uniform);
}

static double shake_uniform01(shake128incctx *st) {
  const uint64_t x = shake_squeeze_u64(st) >> 11;
  return (((double)x) + 0.5) / 9007199254740992.0;
}

static void sample_gaussian_polxvec(polx *r, int64_t *coeffs, size_t len, unsigned int log2stdev,
                                    const uint8_t seed[16], uint64_t nonce) {
  const size_t count = len*N;

  gaussian_lnp_i64(coeffs,count,seed,nonce,log2stdev);
  polxvec_fromint64vec(r,len,1,coeffs);
}

static void sample_ternary_polxvec(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly *tmp = alloc_polyvec(len);

  polyvec_ternary(tmp,len,seed,nonce);
  polxvec_frompolyvec(r,tmp,len);
  free(tmp);
}

static void absorb_bytes(shake128incctx *st, const void *buf, size_t len) {
  shake128_inc_absorb(st,buf,len);
}

static void absorb_u64(shake128incctx *st, uint64_t x) {
  uint8_t buf[8];
  store_u64_le(buf,x);
  absorb_bytes(st,buf,sizeof(buf));
}

static void absorb_polxvec(shake128incctx *st, const polx *v, size_t len) {
  size_t i;
  for(i=0;i<len;i++)
    absorb_polx(st,&v[i]);
}

static void absorb_rbit_const_transcript(shake128incctx *st, const selector_lnp22_ctx *ctx) {
  absorb_u64(st,ctx->rbit_const_lambda);
  if(ctx->rbit_const_lambda == 0)
    return;
  absorb_polxvec(st,ctx->rbit_const_tg,ctx->rbit_const_lambda);
  absorb_polxvec(st,ctx->rbit_const_h,ctx->rbit_const_lambda);
}

static void polxvec_centered_to_i128(__int128 *, const polx *, size_t);

static int validate_rbit_const_transcript(const selector_lnp22_ctx *ctx) {
  size_t row;

  const unsigned __int128 q = ((unsigned __int128)1<<LOGQ)-QOFF;
  const unsigned __int128 b = 2*(unsigned __int128)RBIT_ARP_LIMIT;
  /* Lemma 2.7 bounds every centered scalar by 2*LIMIT. */
  if(ctx->rbit_const_lambda != RBIT_CONST_ROWS+RBIT_ARP_ROWS ||
     (unsigned __int128)ctx->rhat_bits_len*b*(b+1) >= q)
    return 1;
  for(row=RBIT_CONST_ROWS;row<ctx->rbit_const_lambda;row++) {
    __int128 coeffs[N];
    polxvec_centered_to_i128(coeffs,&ctx->rbit_const_h[row],1);
    for(size_t j=0;j<N;j++)
      if(coeffs[j] < -RBIT_ARP_LIMIT || coeffs[j] > RBIT_ARP_LIMIT)
        return 1;
  }
  for(row=0;row<RBIT_CONST_ROWS;row++) {
    if(!polx_iszero_constcoeff(&ctx->rbit_const_h[row])) {
      fprintf(stderr,
              "ERROR: coefficient-constantness response h[%zu] has non-zero constant coefficient\n",
              row);
      return 1;
    }
  }
  return 0;
}

static void bind_prechallenge_commitment(selector_lnp22_ctx *ctx, const polx *tA) {
  shake128incctx st;
  if(!LNP22_DIRECT_CHAIN)
    return;
  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-precommit-v1",21);
  absorb_bytes(&st,ctx->selector_chal_seed,16);
  absorb_polxvec(&st,tA,ctx->kmsis);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(ctx->prechallenge_seed,16,&st);
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->prechallenge_seed,64);
}

static void derive_rbit_const_seed(selector_lnp22_ctx *ctx) {
  shake128incctx st;

  if(ctx->rbit_const_lambda == 0)
    return;

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-rbit-arp-phi-v1",sizeof("selector-rbit-arp-phi-v1"));
  if(LNP22_DIRECT_CHAIN)
    absorb_bytes(&st,ctx->prechallenge_seed,16);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,ctx->rbit_const_lambda);
  absorb_bytes(&st,ctx->selector_chal_seed,sizeof(ctx->selector_chal_seed));
  absorb_polxvec(&st,ctx->rbit_const_tg,ctx->rbit_const_lambda);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(ctx->rbit_const_seed,sizeof(ctx->rbit_const_seed),&st);
}

static void sample_rbit_const_phi(polx *phi, const selector_lnp22_ctx *ctx,
                                  size_t row, size_t bit_idx) {
  if(row >= RBIT_CONST_ROWS) {
    uint8_t buf[AES128CTR_BLOCKBYTES] __attribute__((aligned(64)));
    poly tmp;
    aes128ctr_ctx aes;
    aes128ctr_init(&aes,ctx->rbit_const_seed,
                  0x4152500000000000ULL+row*ctx->rhat_bits_len+bit_idx);
    aes128ctr_squeezeblocks(buf,1,&aes);
    for(size_t j=0;j<N;j++) {
      const unsigned v=(buf[j/4]>>(2*(j%4)))&3;
      tmp.vec->c[j]=(int16_t)((int)(v&1)-(int)((v>>1)&1));
    }
    polxvec_frompolyvec(phi,&tmp,1);
    return;
  }
  sample_zero_const_polx(phi,ctx->rbit_const_seed,
                         NONCE_RBIT_CONST_PHI_BASE +
                         (uint64_t)row * (uint64_t)ctx->rhat_bits_len +
                         (uint64_t)bit_idx);
}

static void sample_x_const_phi(polx *phi, const selector_lnp22_ctx *ctx,
                               size_t row, size_t x_idx) {
  if(row >= RBIT_CONST_ROWS) {
    polxvec_setzero(phi,1);
    return;
  }
  sample_zero_const_polx(phi,ctx->rbit_const_seed,
                         NONCE_X_CONST_PHI_BASE +
                         (uint64_t)row * (uint64_t)ctx->x_len +
                         (uint64_t)x_idx);
}

static void derive_linear_fold_seed(uint8_t out[16], const selector_lnp22_ctx *ctx) {
  shake128incctx st;

  if(ctx->rbit_const_lambda == 0) {
    memcpy(out,ctx->selector_chal_seed,16);
    return;
  }

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-linear-fold",20);
  if(LNP22_DIRECT_CHAIN)
    absorb_bytes(&st,ctx->prechallenge_seed,16);
  absorb_bytes(&st,ctx->selector_chal_seed,sizeof(ctx->selector_chal_seed));
  absorb_rbit_const_transcript(&st,ctx);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(out,16,&st);
}

static void fs_challenge(polx *c, const selector_lnp22_ctx *ctx, const polx *tA, const polx *w,
                         const polx *vlin, const polx *tq, const polx *vq) {
  shake128incctx st;
  uint8_t hash[32];

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-lnp22-fs",16);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,ctx->kmsis);
  absorb_u64(&st,ctx->m2);
  absorb_u64(&st,LNP22_RHAT_BINARY_BACKEND);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  absorb_polxvec(&st,ctx->Cout,ctx->m);
  absorb_rbit_const_transcript(&st,ctx);
  absorb_polxvec(&st,tA,ctx->kmsis);
  absorb_polxvec(&st,w,ctx->kmsis);
  absorb_polx(&st,vlin);
  absorb_polx(&st,tq);
  absorb_polx(&st,vq);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(hash,sizeof(hash),&st);
  sample_autostable_polx(c,hash,0,LNP22_OMEGA);
}

static void stream_rowdot_almostuniform(polx *out, const uint8_t seed[16], uint64_t base_nonce,
                                        size_t row, size_t cols, const polx *vec) {
  const size_t chunks = (cols + 31)/32;
  const size_t row_base = row*chunks;
  size_t blk;
  polx tmp[32];
  polx prod[1];

  polxvec_setzero(out,1);
  for(blk=0;blk<chunks;blk++) {
    const size_t off = 32*blk;
    const size_t len = MIN((size_t)32,cols - off);
    polxvec_almostuniform(tmp,len,seed,base_nonce + row_base + blk);
    polxvec_sprod(prod,tmp,&vec[off],len);
    polx_add(out,out,prod);
  }
}

static void fill_almostuniform_matrix(polx *mat, size_t rows, size_t cols,
                                      const uint8_t seed[16], uint64_t base_nonce) {
  const size_t chunks = (cols + 31)/32;
  size_t row, blk;
  polx tmp[32];

  for(row=0;row<rows;row++) {
    const size_t row_base = row*chunks;
    for(blk=0;blk<chunks;blk++) {
      const size_t off = 32*blk;
      const size_t len = MIN((size_t)32,cols - off);
      polxvec_almostuniform(tmp,len,seed,base_nonce + row_base + blk);
      polxvec_copy(&mat[row*cols + off],tmp,len);
    }
  }
}

static void stream_rowdot_almostuniform_pair(polx *out0, polx *out1, const uint8_t seed[16],
                                             uint64_t base_nonce, size_t row, size_t cols,
                                             const polx *vec0, const polx *vec1) {
  const size_t chunks = (cols + 31)/32;
  const size_t row_base = row*chunks;
  size_t blk;
  polx tmp[32];
  polx prod0[1], prod1[1];

  polxvec_setzero(out0,1);
  polxvec_setzero(out1,1);
  for(blk=0;blk<chunks;blk++) {
    const size_t off = 32*blk;
    const size_t len = MIN((size_t)32,cols - off);
    polxvec_almostuniform(tmp,len,seed,base_nonce + row_base + blk);
    polxvec_sprod(prod0,tmp,&vec0[off],len);
    polxvec_sprod(prod1,tmp,&vec1[off],len);
    polx_add(out0,out0,prod0);
    polx_add(out1,out1,prod1);
  }
}

static void compute_tA_w_rows(polx *tA, polx *w, const selector_lnp22_ctx *ctx,
                              const polx *s1, const polx *s2,
                              const polx *y1, const polx *y2) {
  const polx *s21 = s2;
  const polx *s22 = &s2[ctx->m2prime];
  const polx *y21 = y2;
  const polx *y22 = &y2[ctx->m2prime];
  size_t row;
#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->kmsis;row++) {
    polx acc_tA[1], acc_w[1], tmp_tA[1], tmp_w[1];

    if(ctx->cache_a1open) {
      polxvec_sprod(acc_tA,&ctx->A1open[row*ctx->m1],s1,ctx->m1);
      polxvec_sprod(acc_w,&ctx->A1open[row*ctx->m1],y1,ctx->m1);
    } else {
      stream_rowdot_almostuniform_pair(acc_tA,acc_w,ctx->lnp_seed,NONCE_LNP_A1_BASE,
                                       row,ctx->m1,s1,y1);
    }
    polxvec_sprod(tmp_tA,&ctx->A2open[row*ctx->m2prime],s21,ctx->m2prime);
    polxvec_sprod(tmp_w,&ctx->A2open[row*ctx->m2prime],y21,ctx->m2prime);
    polx_add(acc_tA,acc_tA,tmp_tA);
    polx_add(acc_w,acc_w,tmp_w);
    polx_add(&tA[row],acc_tA,&s22[row]);
    polx_add(&w[row],acc_w,&y22[row]);
  }
}

static void compute_tA_rows(polx *tA, const selector_lnp22_ctx *ctx,
                            const polx *s1, const polx *s2) {
  const polx *s21 = s2;
  const polx *s22 = &s2[ctx->m2prime];
  size_t row;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->kmsis;row++) {
    polx acc_tA[1], tmp_tA[1];

    if(ctx->cache_a1open)
      polxvec_sprod(acc_tA,&ctx->A1open[row*ctx->m1],s1,ctx->m1);
    else
      stream_rowdot_almostuniform(acc_tA,ctx->lnp_seed,NONCE_LNP_A1_BASE,row,ctx->m1,s1);
    polxvec_sprod(tmp_tA,&ctx->A2open[row*ctx->m2prime],s21,ctx->m2prime);
    polx_add(acc_tA,acc_tA,tmp_tA);
    polx_add(&tA[row],acc_tA,&s22[row]);
  }
}

static void compute_w_rows(polx *w, const selector_lnp22_ctx *ctx,
                           const polx *y1, const polx *y2) {
  const polx *y21 = y2;
  const polx *y22 = &y2[ctx->m2prime];
  size_t row;
#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->kmsis;row++) {
    polx acc_w[1];

    if(ctx->cache_a1open)
      polxvec_sprod(acc_w,&ctx->A1open[row*ctx->m1],y1,ctx->m1);
    else
      stream_rowdot_almostuniform(acc_w,ctx->lnp_seed,NONCE_LNP_A1_BASE,row,ctx->m1,y1);
    polxvec_sprod(&w[row],&ctx->A2open[row*ctx->m2prime],y21,ctx->m2prime);
    polx_add(&w[row],&w[row],acc_w);
    polx_add(&w[row],&w[row],&y22[row]);
  }
}

static void compute_a3(polx *a3, const selector_lnp22_ctx *ctx) {
  polxvec_almostuniform(a3,ctx->m2,ctx->lnp_seed,NONCE_LNP_A3_BASE);
}

static void selector_ctx_init(selector_lnp22_ctx *ctx, const uint8_t seed_in[16]) {
  const size_t m = DEMO_M;
  const size_t gdigits = m*gadget_len_base3();
  const size_t Lchain = DEMO_L;
  const size_t rhat_len = UMC_ELL + m;
  /* Every opening coefficient is a scalar bit, packed into its source ring. */
  const size_t opening_dim = LEOPARD_PROOF_FREE ? 0 : (UMC_K + UMC_N + rhat_len)*N;
  const size_t range_rows = LEOPARD_PROOF_FREE ? 0 : UMC_K + UMC_N + rhat_len;
  const size_t rhat_bits_len = LNP22_USE_RBIT_BINARY ?
      rhat_len*N + opening_dim : 0;
  const size_t rhat_pack_rels = LNP22_USE_RBIT_BINARY ? rhat_len + range_rows : 0;
  const size_t rbit_const_lambda = LNP22_ENABLE_RBIT_CONSTCHECK ? LNP22_RBIT_CONSTCHECK_LAMBDA : 0;
  const size_t umc_n = UMC_N;
  const size_t umc_k = UMC_K;
  const size_t umc_v = rhat_len;
  const size_t umc_eua_len = umc_n;
  const size_t umc_eub_len = umc_v;
  const size_t umc_r_len = umc_k;
  const size_t umc_va_len = umc_n;
  const size_t umc_vb_len = umc_v;
  const size_t b_len = Lchain*gdigits;
  const size_t x_len = LNP22_DIRECT_CHAIN ? Lchain : Lchain*SMALL_BLOCK;
  const size_t z_len = LNP22_DIRECT_CHAIN ? 0 : (Lchain-1)*gdigits;
  const size_t y_len = LNP22_DIRECT_CHAIN ? 0 : Lchain*SMALL_BLOCK;
  const size_t m1 = b_len + x_len + z_len + y_len + y_len + rhat_len +
                    umc_eua_len + umc_eub_len + umc_r_len + umc_va_len + umc_vb_len +
                    rhat_bits_len;
  const size_t nlin = (umc_n + umc_v) + rhat_pack_rels + m + 1 +
                      (LNP22_DIRECT_CHAIN ? 0 : 2*(Lchain-1)) + rbit_const_lambda;
  size_t i;

  memset(ctx,0,sizeof(*ctx));
  ctx->Lchain = Lchain;
  ctx->preproc_count = Lchain;
  ctx->selected_j = 0;
  ctx->m = m;
  ctx->gdigits = gdigits;
  ctx->b_len = b_len;
  ctx->x_len = x_len;
  ctx->z_len = z_len;
  ctx->y_len = y_len;
  ctx->rhat_len = rhat_len;
  ctx->rhat_bits_len = rhat_bits_len;
  ctx->umc_n = umc_n;
  ctx->umc_k = umc_k;
  ctx->umc_v = umc_v;
  ctx->umc_eua_len = umc_eua_len;
  ctx->umc_eub_len = umc_eub_len;
  ctx->umc_r_len = umc_r_len;
  ctx->umc_va_len = umc_va_len;
  ctx->umc_vb_len = umc_vb_len;
  ctx->off[LNP22_B_BLOCK] = 0;
  ctx->off[LNP22_X_BLOCK] = b_len;
  ctx->off[LNP22_Z_BLOCK] = b_len + x_len;
  ctx->off[LNP22_Y0_BLOCK] = b_len + x_len + z_len;
  ctx->off[LNP22_Y1_BLOCK] = b_len + x_len + z_len + y_len;
  ctx->off[LNP22_RHAT_BLOCK] = b_len + x_len + z_len + y_len + y_len;
  ctx->off[LNP22_UMC_EUA_BLOCK] = ctx->off[LNP22_RHAT_BLOCK] + rhat_len;
  ctx->off[LNP22_UMC_EUB_BLOCK] = ctx->off[LNP22_UMC_EUA_BLOCK] + umc_eua_len;
  ctx->off[LNP22_UMC_R_BLOCK] = ctx->off[LNP22_UMC_EUB_BLOCK] + umc_eub_len;
  ctx->off[LNP22_UMC_VA_BLOCK] = ctx->off[LNP22_UMC_R_BLOCK] + umc_r_len;
  ctx->off[LNP22_UMC_VB_BLOCK] = ctx->off[LNP22_UMC_VA_BLOCK] + umc_va_len;
  ctx->off[LNP22_RHAT_BITS_BLOCK] = ctx->off[LNP22_UMC_VB_BLOCK] + umc_vb_len;
  ctx->n[LNP22_B_BLOCK] = b_len;
  ctx->n[LNP22_X_BLOCK] = x_len;
  ctx->n[LNP22_Z_BLOCK] = z_len;
  ctx->n[LNP22_Y0_BLOCK] = y_len;
  ctx->n[LNP22_Y1_BLOCK] = y_len;
  ctx->n[LNP22_RHAT_BLOCK] = rhat_len;
  ctx->n[LNP22_UMC_EUA_BLOCK] = umc_eua_len;
  ctx->n[LNP22_UMC_EUB_BLOCK] = umc_eub_len;
  ctx->n[LNP22_UMC_R_BLOCK] = umc_r_len;
  ctx->n[LNP22_UMC_VA_BLOCK] = umc_va_len;
  ctx->n[LNP22_UMC_VB_BLOCK] = umc_vb_len;
  ctx->n[LNP22_RHAT_BITS_BLOCK] = rhat_bits_len;
  ctx->m1 = m1;
  ctx->m2 = LNP22_M2;
  ctx->kmsis = LNP22_KMSIS;
  ctx->m2prime = ctx->m2 - ctx->kmsis;
  ctx->cx_rows = umc_n + umc_v;
  ctx->nlin = nlin;
  ctx->rbit_const_lambda = rbit_const_lambda;
  ctx->cache_a1open = LNP22_CACHE_A1OPEN ? 1 : 0;
  memcpy(ctx->statement_seed,seed_in,16);
  memcpy(ctx->lnp_seed,seed_in,16);

  for(i=0;i<LNP22_NUM_BLOCKS;i++) {
    ctx->witness[i] = alloc_polyvec(ctx->n[i]);
    polyvec_setzero(ctx->witness[i],ctx->n[i]);
  }
  ctx->s1 = alloc_polxvec(ctx->m1);
  ctx->A0 = alloc_polxvec(m*gdigits);
  ctx->A1sel = alloc_polxvec(m*gdigits);
  ctx->Dsel = alloc_polxvec(m*gdigits);
  ctx->Aumc = alloc_polxvec(ctx->umc_n*ctx->umc_k);
  ctx->Bumc = alloc_polxvec(ctx->umc_v*ctx->umc_k);
  ctx->CX = alloc_polxvec(ctx->cx_rows);
  ctx->Arhash = alloc_polxvec(ctx->rhat_len*ctx->m);
  ctx->Cout = alloc_polxvec(ctx->m);
  ctx->A1open = ctx->cache_a1open ? alloc_polxvec(ctx->kmsis*ctx->m1) : NULL;
  ctx->A2open = alloc_polxvec(ctx->kmsis*ctx->m2prime);
  ctx->kvec = alloc_polxvec((Lchain-1)*gdigits);
  ctx->kprime = alloc_polxvec(gdigits);
  ctx->rhovec = alloc_polxvec((Lchain-1)*m);
  ctx->alphavec = alloc_polxvec(Lchain-1);
  ctx->directfold = LNP22_DIRECT_CHAIN ? alloc_polxvec(Lchain-1) : NULL;
  ctx->xbit_weight = LNP22_DIRECT_CHAIN ? alloc_polxvec(x_len) : NULL;
  ctx->rbit_weight = (LNP22_DIRECT_CHAIN && LNP22_USE_RBIT_BINARY) ?
                     alloc_polxvec(rhat_bits_len) : NULL;
  ctx->directD = LNP22_DIRECT_CHAIN ? alloc_polxvec((Lchain-1)*gdigits) : NULL;
  ctx->direct_lin_B = LNP22_DIRECT_CHAIN ? alloc_polxvec(b_len) : NULL;
  ctx->g0 = alloc_polxvec(gdigits);
  ctx->g1 = alloc_polxvec(gdigits);
  ctx->gdelta = alloc_polxvec(gdigits);
  ctx->linfold = alloc_polxvec(nlin);
  ctx->ulin = alloc_polxvec(nlin);
  ctx->lincoeff = alloc_polxvec(ctx->m1);
  ctx->lincoeff2 = alloc_polxvec(ctx->m2);
  ctx->rbit_const_Bg = alloc_polxvec(ctx->rbit_const_lambda * ctx->m2);
  ctx->rbit_const_g = alloc_polxvec(ctx->rbit_const_lambda);
  ctx->rbit_const_tg = alloc_polxvec(ctx->rbit_const_lambda);
  ctx->rbit_const_h = alloc_polxvec(ctx->rbit_const_lambda);
  ctx->native_binary_gamma = LNP22_USE_NATIVE_FOLD_BINARY ?
                             alloc_polxvec(ctx->rbit_const_lambda) : NULL;
  ctx->native_binary_scale = LNP22_USE_NATIVE_FOLD_BINARY ? alloc_polxvec(1) : NULL;
  ctx->native_binary_lin1 = LNP22_USE_NATIVE_FOLD_BINARY ? alloc_polxvec(ctx->x_len) : NULL;
  ctx->native_binary_lin2 = LNP22_USE_NATIVE_FOLD_BINARY ? alloc_polxvec(ctx->m2) : NULL;
  ctx->native_binary_uquad = LNP22_USE_NATIVE_FOLD_BINARY ? alloc_polxvec(1) : NULL;
  ctx->g1coeff = alloc_polxvec(ctx->m1);
  ctx->pre_rhat = alloc_polyvec(ctx->preproc_count * ctx->rhat_len);
  ctx->pre_eua = alloc_polyvec(ctx->preproc_count * ctx->umc_eua_len);
  ctx->pre_eub = alloc_polyvec(ctx->preproc_count * ctx->umc_eub_len);
  ctx->pre_r = alloc_polyvec(ctx->preproc_count * ctx->umc_r_len);
  ctx->pre_va = alloc_polyvec(ctx->preproc_count * ctx->umc_va_len);
  ctx->pre_vb = alloc_polyvec(ctx->preproc_count * ctx->umc_vb_len);
  ctx->pre_CX = alloc_polxvec(ctx->preproc_count * ctx->cx_rows);
  ctx->pre_ar_seed = _malloc(16 * ctx->preproc_count);
  ctx->umc_va_min = INT64_MAX;
  ctx->umc_va_max = INT64_MIN;
  ctx->umc_va_maxabs = 0;
  ctx->umc_vb_min = INT64_MAX;
  ctx->umc_vb_max = INT64_MIN;
  ctx->umc_vb_maxabs = 0;
}

static void selector_ctx_clear(selector_lnp22_ctx *ctx) {
  size_t i;
  for(i=0;i<LNP22_NUM_BLOCKS;i++)
    free(ctx->witness[i]);
  free(ctx->s1);
  free(ctx->A0);
  free(ctx->A1sel);
  free(ctx->Dsel);
  free(ctx->Aumc);
  free(ctx->Bumc);
  free(ctx->CX);
  free(ctx->Arhash);
  free(ctx->Cout);
  free(ctx->A1open);
  free(ctx->A2open);
  free(ctx->kvec);
  free(ctx->kprime);
  free(ctx->rhovec);
  free(ctx->alphavec);
  free(ctx->directfold);
  free(ctx->xbit_weight);
  free(ctx->rbit_weight);
  free(ctx->directD);
  free(ctx->direct_lin_B);
  free(ctx->g0);
  free(ctx->g1);
  free(ctx->gdelta);
  free(ctx->linfold);
  free(ctx->ulin);
  free(ctx->lincoeff);
  free(ctx->lincoeff2);
  free(ctx->rbit_const_Bg);
  free(ctx->rbit_const_g);
  free(ctx->rbit_const_tg);
  free(ctx->rbit_const_h);
  free(ctx->native_binary_gamma);
  free(ctx->native_binary_scale);
  free(ctx->native_binary_lin1);
  free(ctx->native_binary_lin2);
  free(ctx->native_binary_uquad);
  free(ctx->g1coeff);
  free(ctx->pre_rhat);
  free(ctx->pre_eua);
  free(ctx->pre_eub);
  free(ctx->pre_r);
  free(ctx->pre_va);
  free(ctx->pre_vb);
  free(ctx->pre_CX);
  free(ctx->pre_ar_seed);
}

static void free_polxvec_null(polx **p) {
  free(*p);
  *p = NULL;
}

static void trim_ctx_for_batch(selector_lnp22_ctx *ctx) {
  free_polxvec_null(&ctx->A0);
  free_polxvec_null(&ctx->A1sel);
  free_polxvec_null(&ctx->Dsel);
  free_polxvec_null(&ctx->Aumc);
  free_polxvec_null(&ctx->Bumc);
  free_polxvec_null(&ctx->Arhash);
  free_polxvec_null(&ctx->kprime);
  free_polxvec_null(&ctx->rhovec);
  free_polxvec_null(&ctx->g0);
  free_polxvec_null(&ctx->g1);
  free_polxvec_null(&ctx->gdelta);
  free_polxvec_null(&ctx->linfold);
  free_polxvec_null(&ctx->ulin);
  free_polxvec_null(&ctx->g1coeff);
}

static void sync_witness_block_to_s1(selector_lnp22_ctx *ctx, size_t block) {
  polxvec_frompolyvec(&ctx->s1[ctx->off[block]],ctx->witness[block],ctx->n[block]);
}

static const size_t umc_range_blocks[] = {
  LNP22_UMC_R_BLOCK, LNP22_UMC_EUA_BLOCK, LNP22_UMC_EUB_BLOCK
};

static size_t umc_opening_dim(const selector_lnp22_ctx *ctx) {
  return (ctx->umc_r_len+ctx->umc_eua_len+ctx->umc_eub_len)*N;
}

/* All lifted scalars, including the opening, use omega_i*a_i*(a_i-1). */
static void umc_quad_linear(polx *out, const selector_lnp22_ctx *ctx, size_t i) {
  polx_neg(out,&ctx->rbit_weight[i]);
}

static size_t umc_range_rows(const selector_lnp22_ctx *ctx) {
  return ctx->umc_r_len + ctx->umc_eua_len + ctx->umc_eub_len;
}

static void umc_range_row(const selector_lnp22_ctx *ctx, size_t row,
                          size_t *source, size_t *bitoff, size_t *digits,
                          int64_t *bias) {
  *bitoff = ctx->rhat_len*N;
  for(size_t g=0;g<sizeof(umc_range_blocks)/sizeof(*umc_range_blocks);g++) {
    const size_t block = umc_range_blocks[g];
    if(row < ctx->n[block]) {
      *source = ctx->off[block] + row;
      *bitoff += row*N;
      *digits = 1;
      *bias = 0;
      return;
    }
    row -= ctx->n[block];
    *bitoff += ctx->n[block]*N;
  }
  abort();
}

static void umc_range_lhs(polx *out, const selector_lnp22_ctx *ctx,
                          const polx *vec, size_t row) {
  size_t source, bitoff, digits;
  int64_t bias;
  polx term[1];
  umc_range_row(ctx,row,&source,&bitoff,&digits,&bias);
  polxvec_copy(out,&vec[source],1);
  for(size_t h=0;h<N;h++) for(size_t b=0;b<digits;b++) {
    polx_monomial(term,1LL<<b,(int)h);
    polx_mul(term,term,&vec[ctx->off[LNP22_RHAT_BITS_BLOCK]+bitoff+h*digits+b]);
    polx_sub(out,out,term);
  }
}

static void populate_rhat_bits(selector_lnp22_ctx *ctx) {
  size_t row, j;

  if(!LNP22_USE_RBIT_BINARY)
    return;

  for(row=0;row<ctx->rhat_len;row++) {
    const poly *src = &ctx->witness[LNP22_RHAT_BLOCK][row];
    for(j=0;j<N;j++) {
      if(src->vec->c[j] != 0 && src->vec->c[j] != 1) abort();
      poly_set_const(&ctx->witness[LNP22_RHAT_BITS_BLOCK][row*N + j],
                     src->vec->c[j]);
    }
  }
  size_t bitoff = ctx->rhat_len*N;
  for(size_t g=0;g<sizeof(umc_range_blocks)/sizeof(*umc_range_blocks);g++) {
    const size_t block = umc_range_blocks[g];
    for(row=0;row<ctx->n[block];row++) for(j=0;j<N;j++) {
      const int v = ctx->witness[block][row].vec->c[j];
      if(v != 0 && v != 1) {
        fprintf(stderr,"UMC opening coefficient is not binary\n");
        abort();
      }
      poly_set_const(&ctx->witness[LNP22_RHAT_BITS_BLOCK][bitoff++],v);
    }
  }
  if(bitoff != ctx->rhat_bits_len) abort();
  sync_witness_block_to_s1(ctx,LNP22_RHAT_BITS_BLOCK);
}

static uint64_t preproc_nonce(uint64_t base, size_t j, size_t idx) {
  return base + 4096ULL*j + (uint64_t)idx;
}

static void sample_umc_matrices(selector_lnp22_ctx *ctx) {
  sample_uniform_modp_polxvec(ctx->Aumc,ctx->umc_n*ctx->umc_k,ctx->lnp_seed,NONCE_UMC_A_BASE,UMC_MODP);
  sample_uniform_modp_polxvec(ctx->Bumc,ctx->umc_v*ctx->umc_k,ctx->lnp_seed,NONCE_UMC_B_BASE,UMC_MODP);
}

static void compute_umc_commitment(selector_lnp22_ctx *ctx) {
  const int64_t p = UMC_MODP;
  polx prod[1], tmpmsg[1];
  int64_t coeffs[N];
  int64_t rems[N];
  size_t row, j;

  sync_witness_block_to_s1(ctx,LNP22_UMC_EUA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_EUB_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_R_BLOCK);

  ctx->umc_va_min = INT64_MAX;
  ctx->umc_va_max = INT64_MIN;
  ctx->umc_va_maxabs = 0;
  ctx->umc_vb_min = INT64_MAX;
  ctx->umc_vb_max = INT64_MIN;
  ctx->umc_vb_maxabs = 0;

  for(row=0;row<ctx->umc_n;row++) {
    polxvec_sprod(prod,&ctx->Aumc[row*ctx->umc_k],&ctx->s1[ctx->off[LNP22_UMC_R_BLOCK]],ctx->umc_k);
    polx_add(prod,prod,&ctx->s1[ctx->off[LNP22_UMC_EUA_BLOCK] + row]);
    polxvec_centered_to_i64(coeffs,prod,1);
    for(j=0;j<N;j++) {
      int64_t rem;
      const int64_t q = centered_div_i64(coeffs[j],p,&rem);
      ctx->witness[LNP22_UMC_VA_BLOCK][row].vec->c[j] = (int16_t)q;
      rems[j] = rem;
      if(q < ctx->umc_va_min) ctx->umc_va_min = q;
      if(q > ctx->umc_va_max) ctx->umc_va_max = q;
      if((uint64_t)(q < 0 ? -q : q) > ctx->umc_va_maxabs)
        ctx->umc_va_maxabs = (uint64_t)(q < 0 ? -q : q);
    }
    polxvec_fromint64vec(&ctx->CX[row],1,1,rems);
  }

  for(row=0;row<ctx->umc_v;row++) {
    const polx *msg = &ctx->s1[ctx->off[LNP22_RHAT_BLOCK] + row];
    polxvec_sprod(prod,&ctx->Bumc[row*ctx->umc_k],&ctx->s1[ctx->off[LNP22_UMC_R_BLOCK]],ctx->umc_k);
    polx_add(prod,prod,&ctx->s1[ctx->off[LNP22_UMC_EUB_BLOCK] + row]);
    polx_scale(tmpmsg,msg,UMC_MSG_SCALE);
    polx_add(prod,prod,tmpmsg);
    polxvec_centered_to_i64(coeffs,prod,1);
    for(j=0;j<N;j++) {
      int64_t rem;
      const int64_t q = centered_div_i64(coeffs[j],p,&rem);
      ctx->witness[LNP22_UMC_VB_BLOCK][row].vec->c[j] = (int16_t)q;
      rems[j] = rem;
      if(q < ctx->umc_vb_min) ctx->umc_vb_min = q;
      if(q > ctx->umc_vb_max) ctx->umc_vb_max = q;
      if((uint64_t)(q < 0 ? -q : q) > ctx->umc_vb_maxabs)
        ctx->umc_vb_maxabs = (uint64_t)(q < 0 ? -q : q);
    }
    polxvec_fromint64vec(&ctx->CX[ctx->umc_n + row],1,1,rems);
  }

  sync_witness_block_to_s1(ctx,LNP22_UMC_VA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_VB_BLOCK);
}

static void derive_ar_seed_from_commitment(uint8_t out[16], const selector_lnp22_ctx *ctx,
                                           const polx *CX) {
  shake128incctx st;

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-umc-ar",15);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,UMC_MODP);
  absorb_u64(&st,ctx->umc_n);
  absorb_u64(&st,ctx->umc_k);
  absorb_u64(&st,ctx->umc_v);
  absorb_u64(&st,ctx->selected_j);
  absorb_polxvec(&st,CX,ctx->cx_rows);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(out,16,&st);
}

static void update_umc_v_stats(selector_lnp22_ctx *ctx) {
  size_t row, j;

  ctx->umc_va_min = INT64_MAX;
  ctx->umc_va_max = INT64_MIN;
  ctx->umc_va_maxabs = 0;
  ctx->umc_vb_min = INT64_MAX;
  ctx->umc_vb_max = INT64_MIN;
  ctx->umc_vb_maxabs = 0;

  for(row=0;row<ctx->umc_va_len;row++) {
    for(j=0;j<N;j++) {
      const int64_t q = ctx->witness[LNP22_UMC_VA_BLOCK][row].vec->c[j];
      const uint64_t a = (uint64_t)(q < 0 ? -q : q);
      if(q < ctx->umc_va_min) ctx->umc_va_min = q;
      if(q > ctx->umc_va_max) ctx->umc_va_max = q;
      if(a > ctx->umc_va_maxabs) ctx->umc_va_maxabs = a;
    }
  }
  for(row=0;row<ctx->umc_vb_len;row++) {
    for(j=0;j<N;j++) {
      const int64_t q = ctx->witness[LNP22_UMC_VB_BLOCK][row].vec->c[j];
      const uint64_t a = (uint64_t)(q < 0 ? -q : q);
      if(q < ctx->umc_vb_min) ctx->umc_vb_min = q;
      if(q > ctx->umc_vb_max) ctx->umc_vb_max = q;
      if(a > ctx->umc_vb_maxabs) ctx->umc_vb_maxabs = a;
    }
  }
}

static void precompute_umc_commitments(selector_lnp22_ctx *ctx) {
  size_t j;

  sample_umc_matrices(ctx);
  for(j=0;j<ctx->preproc_count;j++) {
    size_t row;
    uint8_t *ar_seed = &ctx->pre_ar_seed[16*j];

    for(row=0;row<ctx->rhat_len;row++)
      sample_binary_poly(&ctx->pre_rhat[j*ctx->rhat_len + row],ctx->lnp_seed,preproc_nonce(NONCE_RHAT_BASE,j,row));
    for(row=0;row<ctx->umc_eua_len;row++)
      (LEOPARD_PROOF_FREE ? sample_cbd1_poly : sample_binary_poly)(&ctx->pre_eua[j*ctx->umc_eua_len + row],ctx->lnp_seed,preproc_nonce(NONCE_UMC_EUA_BASE,j,row));
    for(row=0;row<ctx->umc_eub_len;row++)
      (LEOPARD_PROOF_FREE ? sample_cbd1_poly : sample_binary_poly)(&ctx->pre_eub[j*ctx->umc_eub_len + row],ctx->lnp_seed,preproc_nonce(NONCE_UMC_EUB_BASE,j,row));
    for(row=0;row<ctx->umc_r_len;row++)
      (LEOPARD_PROOF_FREE ? sample_cbd1_poly : sample_binary_poly)(&ctx->pre_r[j*ctx->umc_r_len + row],ctx->lnp_seed,preproc_nonce(NONCE_UMC_R_BASE,j,row));

    polyvec_copy(ctx->witness[LNP22_RHAT_BLOCK],&ctx->pre_rhat[j*ctx->rhat_len],ctx->rhat_len);
    polyvec_copy(ctx->witness[LNP22_UMC_EUA_BLOCK],&ctx->pre_eua[j*ctx->umc_eua_len],ctx->umc_eua_len);
    polyvec_copy(ctx->witness[LNP22_UMC_EUB_BLOCK],&ctx->pre_eub[j*ctx->umc_eub_len],ctx->umc_eub_len);
    polyvec_copy(ctx->witness[LNP22_UMC_R_BLOCK],&ctx->pre_r[j*ctx->umc_r_len],ctx->umc_r_len);
    sync_witness_block_to_s1(ctx,LNP22_RHAT_BLOCK);
    compute_umc_commitment(ctx);
    polyvec_copy(&ctx->pre_va[j*ctx->umc_va_len],ctx->witness[LNP22_UMC_VA_BLOCK],ctx->umc_va_len);
    polyvec_copy(&ctx->pre_vb[j*ctx->umc_vb_len],ctx->witness[LNP22_UMC_VB_BLOCK],ctx->umc_vb_len);
    polxvec_copy(&ctx->pre_CX[j*ctx->cx_rows],ctx->CX,ctx->cx_rows);
    ctx->selected_j = j;
    derive_ar_seed_from_commitment(ar_seed,ctx,&ctx->pre_CX[j*ctx->cx_rows]);
  }
}

static void load_preprocessed_commitment(selector_lnp22_ctx *ctx, size_t j) {
  ctx->selected_j = j;
  polyvec_copy(ctx->witness[LNP22_RHAT_BLOCK],&ctx->pre_rhat[j*ctx->rhat_len],ctx->rhat_len);
  polyvec_copy(ctx->witness[LNP22_UMC_EUA_BLOCK],&ctx->pre_eua[j*ctx->umc_eua_len],ctx->umc_eua_len);
  polyvec_copy(ctx->witness[LNP22_UMC_EUB_BLOCK],&ctx->pre_eub[j*ctx->umc_eub_len],ctx->umc_eub_len);
  polyvec_copy(ctx->witness[LNP22_UMC_R_BLOCK],&ctx->pre_r[j*ctx->umc_r_len],ctx->umc_r_len);
  polyvec_copy(ctx->witness[LNP22_UMC_VA_BLOCK],&ctx->pre_va[j*ctx->umc_va_len],ctx->umc_va_len);
  polyvec_copy(ctx->witness[LNP22_UMC_VB_BLOCK],&ctx->pre_vb[j*ctx->umc_vb_len],ctx->umc_vb_len);
  sync_witness_block_to_s1(ctx,LNP22_RHAT_BLOCK);
  populate_rhat_bits(ctx);
  sync_witness_block_to_s1(ctx,LNP22_UMC_EUA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_EUB_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_R_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_VA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_VB_BLOCK);
  polxvec_copy(ctx->CX,&ctx->pre_CX[j*ctx->cx_rows],ctx->cx_rows);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,&ctx->pre_ar_seed[16*j],NONCE_AR_BASE);
  update_umc_v_stats(ctx);
}

static void derive_selector_chal_seed_from_public(selector_lnp22_ctx *ctx) {
  shake128incctx st;

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-nizk1-stmt",19);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,UMC_MODP);
  absorb_u64(&st,ctx->umc_n);
  absorb_u64(&st,ctx->umc_k);
  absorb_u64(&st,ctx->umc_v);
  absorb_u64(&st,ctx->selected_j);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  absorb_polxvec(&st,ctx->Cout,ctx->m);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(ctx->selector_chal_seed,sizeof(ctx->selector_chal_seed),&st);
}

static void compute_Cout(selector_lnp22_ctx *ctx) {
  polx prod[1], gb[1];
  size_t row, j;

  for(row=0;row<ctx->m;row++) {
    polxvec_sprod(prod,&ctx->Arhash[row*ctx->rhat_len],&ctx->s1[ctx->off[LNP22_RHAT_BLOCK]],ctx->rhat_len);
    polxvec_setzero(gb,1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(gb,&ctx->s1[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
    polx_add(&ctx->Cout[row],prod,gb);
  }
}

static int reject_standard_polxvec(const int64_t *, const int64_t *, size_t,
                                   double, double, const uint8_t[16], uint64_t);

static int prepare_rbit_constcheck(selector_lnp22_ctx *ctx, const polx *s2) {
  size_t row;
  polx native_value[1];
  int64_t yarp[RBIT_ARP_ROWS*N], zarp[RBIT_ARP_ROWS*N], shift[RBIT_ARP_ROWS*N];

  if(ctx->rbit_const_lambda == 0)
    return 0;

  polxvec_almostuniform(ctx->rbit_const_Bg,ctx->rbit_const_lambda * ctx->m2,
                        ctx->selector_chal_seed,NONCE_RBIT_CONST_BG_BASE);
  if(RBIT_ARP_SIGMA < 10*sqrt(RBIT_ARP_ROWS*N)*ctx->rhat_bits_len)
    return 1;
  for(size_t attempt=0;attempt<10000;attempt++) {
  sample_gaussian_polxvec(&ctx->rbit_const_g[RBIT_CONST_ROWS],yarp,
                          RBIT_ARP_ROWS,RBIT_ARP_LOGSIGMA,ctx->lnp_seed,
                          0x4152504d00000000ULL+attempt);
  for(row=0;row<ctx->rbit_const_lambda;row++) {
    if(row<RBIT_CONST_ROWS)
      sample_zero_const_polx(&ctx->rbit_const_g[row],ctx->lnp_seed,
                             NONCE_RBIT_CONST_G_BASE + row);
    polxvec_sprod(&ctx->rbit_const_tg[row],
                  &ctx->rbit_const_Bg[row*ctx->m2],s2,ctx->m2);
    polx_add(&ctx->rbit_const_tg[row],
             &ctx->rbit_const_tg[row],&ctx->rbit_const_g[row]);
  }

  derive_rbit_const_seed(ctx);
  if(LNP22_USE_NATIVE_FOLD_BINARY)
    rhat_native_binary_value(native_value,ctx,
                             &ctx->s1[ctx->off[LNP22_RHAT_BLOCK]]);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->rbit_const_lambda;row++) {
    polx phi[1], term[1];
    polxvec_copy(&ctx->rbit_const_h[row],&ctx->rbit_const_g[row],1);
    if(LNP22_USE_RBIT_BINARY) {
      for(size_t i=0;i<ctx->rhat_bits_len;i++) {
        sample_rbit_const_phi(phi,ctx,row,i);
        polx_mul(term,phi,&ctx->s1[ctx->off[LNP22_RHAT_BITS_BLOCK] + i]);
        polx_add(&ctx->rbit_const_h[row],&ctx->rbit_const_h[row],term);
      }
    }
    if(LNP22_PROVE_X_BITS) {
      for(size_t i=0;i<ctx->x_len;i++) {
        sample_x_const_phi(phi,ctx,row,i);
        polx_mul(term,phi,&ctx->s1[ctx->off[LNP22_X_BLOCK] + i]);
        polx_add(&ctx->rbit_const_h[row],&ctx->rbit_const_h[row],term);
      }
    }
    if(LNP22_USE_NATIVE_FOLD_BINARY) {
      sample_full_scalar_challenge(&ctx->native_binary_gamma[row],
                                   ctx->rbit_const_seed,
                                   NONCE_NATIVE_BINARY_GAMMA_BASE + row);
      polx_mul(term,&ctx->native_binary_gamma[row],native_value);
      polx_add(&ctx->rbit_const_h[row],&ctx->rbit_const_h[row],term);
    }
  }
  polxvec_centered_to_i64(zarp,&ctx->rbit_const_h[RBIT_CONST_ROWS],RBIT_ARP_ROWS);
  for(size_t i=0;i<RBIT_ARP_ROWS*N;i++) shift[i]=zarp[i]-yarp[i];
  if(validate_rbit_const_transcript(ctx) ||
     reject_standard_polxvec(zarp,shift,RBIT_ARP_ROWS*N,RBIT_ARP_SIGMA,10,
                             ctx->lnp_seed,0x4152505200000000ULL+attempt))
    continue;
  printf("  ARP mask rounds   = %zu\n",attempt+1);
  return 0;
  }
  return 1;
}

static void compute_linear_constants(selector_lnp22_ctx *ctx) {
  size_t idx = 0;
  size_t row, t;
  polx c0[1];

  polxvec_sprod(c0,ctx->kprime,ctx->g0,ctx->gdigits);

  for(row=0;row<ctx->cx_rows;row++) {
    polxvec_copy(&ctx->ulin[idx++],&ctx->CX[row],1);
  }
  if(LNP22_USE_RBIT_BINARY) {
    for(row=0;row<ctx->rhat_len;row++) {
      polxvec_setzero(&ctx->ulin[idx++],1);
    }
    for(row=0;row<umc_range_rows(ctx);row++) {
      size_t source, bitoff, digits;
      int64_t bias, coeffs[N];
      umc_range_row(ctx,row,&source,&bitoff,&digits,&bias);
      for(size_t h=0;h<N;h++) coeffs[h] = digits ? -bias : (h==0 ? 1 : 0);
      polxvec_fromint64vec(&ctx->ulin[idx++],1,1,coeffs);
    }
  }
  for(row=0;row<ctx->m;row++) {
    polxvec_copy(&ctx->ulin[idx++],&ctx->Cout[row],1);
  }
  polxvec_copy(&ctx->ulin[idx++],c0,1);

#if !LNP22_DIRECT_CHAIN
  for(t=0;t<ctx->Lchain-1;t++)
    polxvec_setzero(&ctx->ulin[idx++],1);

  for(t=0;t<ctx->Lchain-1;t++)
    polxvec_setzero(&ctx->ulin[idx++],1);
#else
  (void)t;
#endif

  for(row=0;row<ctx->rbit_const_lambda;row++) {
    if(LNP22_USE_NATIVE_FOLD_BINARY)
      polxvec_setzero(&ctx->ulin[idx],1);
    else
      polx_sub(&ctx->ulin[idx],&ctx->rbit_const_h[row],&ctx->rbit_const_tg[row]);
    idx++;
  }
}

static void compute_public_fold_rhs(polx *ufold, const selector_lnp22_ctx *ctx) {
  size_t idx;
  polx tmp[1];

  polxvec_setzero(ufold,1);
  for(idx=0;idx<ctx->nlin;idx++) {
    polx_mul(tmp,&ctx->linfold[idx],&ctx->ulin[idx]);
    polx_add(ufold,ufold,tmp);
  }
}

static void compute_linear_fold_constants(selector_lnp22_ctx *ctx) {
  uint8_t seed[16];

  derive_linear_fold_seed(seed,ctx);
  polxvec_challenge(ctx->linfold,ctx->nlin,seed,NONCE_LINFOLD_BASE);
}

static void compute_native_binary_fold(selector_lnp22_ctx *ctx) {
  const size_t base = ctx->nlin - ctx->rbit_const_lambda;
  size_t row, i;
  polx rho_gamma[1], factor[1], phi[1], rhs[1];

  if(!LNP22_USE_NATIVE_FOLD_BINARY)
    return;

  polxvec_setzero(ctx->native_binary_scale,1);
  polxvec_setzero(ctx->native_binary_lin1,ctx->x_len);
  polxvec_setzero(ctx->native_binary_lin2,ctx->m2);
  polxvec_setzero(ctx->native_binary_uquad,1);
  for(row=0;row<ctx->rbit_const_lambda;row++) {
    const polx *rho = &ctx->linfold[base + row];

    polx_mul(rho_gamma,rho,&ctx->native_binary_gamma[row]);
    polx_add(ctx->native_binary_scale,ctx->native_binary_scale,rho_gamma);
    if(LNP22_PROVE_X_BITS) {
      for(i=0;i<ctx->x_len;i++) {
        sample_x_const_phi(phi,ctx,row,i);
        polx_mul(factor,rho,phi);
        polx_add(&ctx->native_binary_lin1[i],&ctx->native_binary_lin1[i],factor);
      }
    }
    for(i=0;i<ctx->m2;i++) {
      polx_mul(factor,rho,&ctx->rbit_const_Bg[row*ctx->m2 + i]);
      polx_sub(&ctx->native_binary_lin2[i],&ctx->native_binary_lin2[i],factor);
    }
    polx_sub(rhs,&ctx->rbit_const_h[row],&ctx->rbit_const_tg[row]);
    polx_mul(rhs,rho,rhs);
    polx_add(ctx->native_binary_uquad,ctx->native_binary_uquad,rhs);
  }
  polx_refresh(ctx->native_binary_scale);
  polxvec_refresh(ctx->native_binary_lin1,ctx->x_len);
  polxvec_refresh(ctx->native_binary_lin2,ctx->m2);
  polx_refresh(ctx->native_binary_uquad);
}

static void compute_public_quad_rhs(polx *uquad, const selector_lnp22_ctx *ctx) {
  if(LNP22_USE_NATIVE_FOLD_BINARY)
    polxvec_copy(uquad,ctx->native_binary_uquad,1);
  else
    polxvec_setzero(uquad,1);
}

static void compute_folded_linear_coeffs(selector_lnp22_ctx *ctx) {
  size_t idx = 0;
  size_t row, t, j, i;
  polx factor[1], cdelta[1], monomial[1];
#if LNP22_LINEAR_PROFILE
  struct timespec coeff_prof_last, coeff_prof_now;
#define COEFF_PROFILE_START() clock_gettime(CLOCK_MONOTONIC,&coeff_prof_last)
#define COEFF_PROFILE_MARK(label) do {                                      \
    clock_gettime(CLOCK_MONOTONIC,&coeff_prof_now);                         \
    printf("    folded coeffs: %-22s %.6fs\n",label,                       \
           elapsed_sec(&coeff_prof_last,&coeff_prof_now));                  \
    coeff_prof_last = coeff_prof_now;                                        \
  } while(0)
#else
#define COEFF_PROFILE_START() do { } while(0)
#define COEFF_PROFILE_MARK(label) do { } while(0)
#endif
  const size_t c_idx = ctx->cx_rows;
  const size_t pack_idx = c_idx;
  const size_t cout_idx = pack_idx +
                          (LNP22_USE_RBIT_BINARY ? ctx->rhat_len+umc_range_rows(ctx) : 0);
  const size_t top_idx = cout_idx + ctx->m;
  const size_t y_idx = top_idx + 1;
  const size_t tr_idx = y_idx + (ctx->Lchain-1);

  COEFF_PROFILE_START();
  polxvec_setzero(ctx->lincoeff,ctx->m1);
  polxvec_setzero(ctx->lincoeff2,ctx->m2);
  polxvec_sprod(cdelta,ctx->kprime,ctx->gdelta,ctx->gdigits);

  for(row=0;row<ctx->umc_n;row++) {
    for(i=0;i<ctx->umc_k;i++) {
      polx_mul(factor,&ctx->linfold[idx],&ctx->Aumc[row*ctx->umc_k + i]);
      polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK] + i],
               &ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK] + i],factor);
    }
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_EUA_BLOCK] + row],
             &ctx->lincoeff[ctx->off[LNP22_UMC_EUA_BLOCK] + row],&ctx->linfold[idx]);
    polx_scale(factor,&ctx->linfold[idx],-UMC_MODP);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_VA_BLOCK] + row],
             &ctx->lincoeff[ctx->off[LNP22_UMC_VA_BLOCK] + row],factor);
    idx++;
  }

  COEFF_PROFILE_MARK("UMC A coeffs");

  for(t=0;t<ctx->umc_v;t++) {
    polx *msgcoeff = &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + t];
    for(i=0;i<ctx->umc_k;i++) {
      polx_mul(factor,&ctx->linfold[idx],&ctx->Bumc[t*ctx->umc_k + i]);
      polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK] + i],
               &ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK] + i],factor);
    }
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_EUB_BLOCK] + t],
             &ctx->lincoeff[ctx->off[LNP22_UMC_EUB_BLOCK] + t],&ctx->linfold[idx]);
    polx_scale(factor,&ctx->linfold[idx],UMC_MSG_SCALE);
    polx_add(msgcoeff,msgcoeff,factor);
    polx_scale(factor,&ctx->linfold[idx],-UMC_MODP);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_VB_BLOCK] + t],
             &ctx->lincoeff[ctx->off[LNP22_UMC_VB_BLOCK] + t],factor);
    idx++;
  }

  COEFF_PROFILE_MARK("UMC B coeffs");

  if(LNP22_USE_RBIT_BINARY) {
    for(row=0;row<ctx->rhat_len;row++) {
      polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + row],
               &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + row],
               &ctx->linfold[idx]);
      for(j=0;j<N;j++) {
        polx_monomial(monomial,-1,(int)j);
        polx_mul(factor,&ctx->linfold[idx],monomial);
        polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + row*N + j],
                 &ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + row*N + j],
                 factor);
      }
      idx++;
    }
  }

  if(LNP22_USE_RBIT_BINARY) for(row=0;row<umc_range_rows(ctx);row++) {
    size_t source, bitoff, digits;
    int64_t bias;
    umc_range_row(ctx,row,&source,&bitoff,&digits,&bias);
    polx_add(&ctx->lincoeff[source],&ctx->lincoeff[source],&ctx->linfold[idx]);
    for(j=0;j<N;j++) for(size_t b=0;b<digits;b++) {
      const size_t dest = ctx->off[LNP22_RHAT_BITS_BLOCK]+bitoff+j*digits+b;
      polx_monomial(monomial,-(1LL<<b),(int)j);
      polx_mul(factor,&ctx->linfold[idx],monomial);
      polx_add(&ctx->lincoeff[dest],&ctx->lincoeff[dest],factor);
    }
    idx++;
  }

  COEFF_PROFILE_MARK("RBIT/range packing");

  for(row=0;row<ctx->m;row++) {
    for(i=0;i<ctx->rhat_len;i++) {
      polx_mul(factor,&ctx->linfold[idx],&ctx->Arhash[row*ctx->rhat_len + i]);
      polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + i],
               &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + i],factor);
    }
    for(j=0;j<ctx->gdigits/ctx->m;j++) {
      polx_scale_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],
                     &ctx->linfold[idx],pow3_u64(j));
    }
    idx++;
  }

  COEFF_PROFILE_MARK("C wellformed");

  for(i=0;i<ctx->gdigits;i++) {
    polx_mul(factor,&ctx->linfold[top_idx],&ctx->kprime[i]);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (ctx->Lchain-1)*ctx->gdigits + i],
             &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (ctx->Lchain-1)*ctx->gdigits + i],factor);
  }
  polx_mul(factor,&ctx->linfold[top_idx],cdelta);
  polx_sub(&ctx->lincoeff[ctx->off[LNP22_X_BLOCK] +
                         selector_x_pos(ctx,ctx->Lchain-1)],
           &ctx->lincoeff[ctx->off[LNP22_X_BLOCK] +
                          selector_x_pos(ctx,ctx->Lchain-1)],factor);

  COEFF_PROFILE_MARK("top relation");

#if !LNP22_DIRECT_CHAIN
  for(t=0;t<ctx->Lchain-1;t++) {
    for(i=0;i<ctx->gdigits;i++) {
      polx_mul(factor,&ctx->alphavec[t],&ctx->kvec[t*ctx->gdigits + i]);
      polx_mul(factor,factor,&ctx->linfold[y_idx + t]);
      polx_sub(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],
               &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],factor);
    }
    polx_add(&ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK] + t*SMALL_BLOCK],
             &ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK] + t*SMALL_BLOCK],&ctx->linfold[y_idx + t]);
    polx_scale(factor,&ctx->linfold[y_idx + t],Y_SPLIT_GAMMA);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK] + t*SMALL_BLOCK],
             &ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK] + t*SMALL_BLOCK],factor);
  }

  COEFF_PROFILE_MARK("Y split coeffs");

#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
  for(size_t parity=0;parity<2;parity++) {
#pragma omp parallel for schedule(static)
    for(long long tt=(long long)parity;tt<(long long)(ctx->Lchain-1);tt+=2) {
      const size_t tloc = (size_t)tt;
      size_t rowloc, jloc, iloc;
      polx factorloc[1], tmploc[1];

      for(rowloc=0;rowloc<ctx->m;rowloc++) {
        polx_mul(factorloc,&ctx->linfold[tr_idx + tloc],&ctx->rhovec[tloc*ctx->m + rowloc]);
        for(jloc=0;jloc<ctx->gdigits/ctx->m;jloc++) {
          polx_scale_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + tloc*ctx->gdigits + rowloc*(ctx->gdigits/ctx->m) + jloc],
                         factorloc,pow3_u64(jloc));
        }
        for(iloc=0;iloc<ctx->gdigits;iloc++) {
          polx_mul(tmploc,factorloc,&ctx->A0[rowloc*ctx->gdigits + iloc]);
          polx_sub(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (tloc+1)*ctx->gdigits + iloc],
                   &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (tloc+1)*ctx->gdigits + iloc],tmploc);
          polx_mul(tmploc,factorloc,&ctx->Dsel[rowloc*ctx->gdigits + iloc]);
          polx_sub(&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + tloc*ctx->gdigits + iloc],
                   &ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + tloc*ctx->gdigits + iloc],tmploc);
        }
      }
    }
  }
#else
  for(t=0;t<ctx->Lchain-1;t++) {
    for(row=0;row<ctx->m;row++) {
      polx_mul(factor,&ctx->linfold[tr_idx + t],&ctx->rhovec[t*ctx->m + row]);
      for(j=0;j<ctx->gdigits/ctx->m;j++) {
        polx_scale_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + t*ctx->gdigits + row*(ctx->gdigits/ctx->m) + j],
                       factor,pow3_u64(j));
      }
      for(i=0;i<ctx->gdigits;i++) {
        polx tmp[1];
        polx_mul(tmp,factor,&ctx->A0[row*ctx->gdigits + i]);
        polx_sub(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],
                 &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],tmp);
        polx_mul(tmp,factor,&ctx->Dsel[row*ctx->gdigits + i]);
        polx_sub(&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits + i],
                 &ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits + i],tmp);
      }
    }
  }
#endif

  COEFF_PROFILE_MARK("transition coeffs");
#else
  COEFF_PROFILE_MARK("direct chain moved to quadratic");
#endif

  idx = LNP22_DIRECT_CHAIN ? top_idx + 1 : tr_idx + (ctx->Lchain-1);
  if(ctx->rbit_const_lambda != 0 && !LNP22_USE_NATIVE_FOLD_BINARY) {
    polx phi[1];
    /* Each scalar owns one output coefficient; fold rows locally to avoid
     * both write conflicts and a parallel region for every row. */
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
    for(size_t bit=0;bit<ctx->rhat_bits_len;bit++) {
      polx bit_phi[1], bit_factor[1];
      for(size_t crow=0;crow<ctx->rbit_const_lambda;crow++) {
        sample_rbit_const_phi(bit_phi,ctx,crow,bit);
        polx_mul(bit_factor,&ctx->linfold[idx+crow],bit_phi);
        polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK]+bit],
                 &ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK]+bit],bit_factor);
      }
    }
    for(row=0;row<ctx->rbit_const_lambda;row++) {
      if(LNP22_PROVE_X_BITS) {
        for(i=0;i<ctx->x_len;i++) {
          sample_x_const_phi(phi,ctx,row,i);
          polx_mul(factor,&ctx->linfold[idx],phi);
          polx_add(&ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + i],
                   &ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + i],factor);
        }
      }
      for(i=0;i<ctx->m2;i++) {
        polx_mul(factor,&ctx->linfold[idx],
                 &ctx->rbit_const_Bg[row*ctx->m2 + i]);
        polx_sub(&ctx->lincoeff2[i],&ctx->lincoeff2[i],factor);
      }
      idx++;
    }
  }
  COEFF_PROFILE_MARK("constcheck coeffs");
#undef COEFF_PROFILE_START
#undef COEFF_PROFILE_MARK
}

static int prepare_linear_fold_for_proof(selector_lnp22_ctx *ctx, polx *ufold_public,
                                         const polx *s2) {
#if LNP22_LINEAR_PROFILE
  struct timespec lin_prof_last, lin_prof_now;
#define LIN_PROFILE_START() clock_gettime(CLOCK_MONOTONIC,&lin_prof_last)
#define LIN_PROFILE_MARK(label) do {                                        \
    clock_gettime(CLOCK_MONOTONIC,&lin_prof_now);                           \
    printf("  linear setup profile: %-20s %.6fs\n",label,                  \
           elapsed_sec(&lin_prof_last,&lin_prof_now));                      \
    lin_prof_last = lin_prof_now;                                            \
  } while(0)
#else
#define LIN_PROFILE_START() do { } while(0)
#define LIN_PROFILE_MARK(label) do { } while(0)
#endif

  LIN_PROFILE_START();
  if(prepare_rbit_constcheck(ctx,s2))
    return 1;
  LIN_PROFILE_MARK("rbit constcheck");
  compute_linear_constants(ctx);
  LIN_PROFILE_MARK("linear constants");
  compute_linear_fold_constants(ctx);
  LIN_PROFILE_MARK("fold challenges");
  compute_native_binary_fold(ctx);
  LIN_PROFILE_MARK("native binary fold");
  compute_folded_linear_coeffs(ctx);
  LIN_PROFILE_MARK("fold coeffs");
  compute_public_fold_rhs(ufold_public,ctx);
  LIN_PROFILE_MARK("public fold rhs");
#undef LIN_PROFILE_START
#undef LIN_PROFILE_MARK
  return 0;
}

static void derive_direct_chain_fold(selector_lnp22_ctx *ctx, const polx *tA) {
  shake128incctx st;
  uint8_t seeds[48];
  size_t t;

  if(!LNP22_DIRECT_CHAIN)
    return;
  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-direct-chain-fold-v2",29);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  absorb_polxvec(&st,ctx->Cout,ctx->m);
  absorb_rbit_const_transcript(&st,ctx);
  absorb_polxvec(&st,tA,ctx->kmsis);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(seeds,sizeof(seeds),&st);
  polxvec_challenge(ctx->directfold,ctx->Lchain-1,seeds,0);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(t=0;t<ctx->Lchain-1;t++)
    polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,seeds,1 + (uint64_t)t);
  polxvec_almostuniform(ctx->xbit_weight,ctx->x_len,&seeds[16],0);
  if(LNP22_USE_RBIT_BINARY) {
    /* The separately checked ARP supplies the no-wrap bound for this sum. */
    for(size_t i=0;i<ctx->rhat_bits_len;i++)
      polx_monomial(&ctx->rbit_weight[i],1,0);
  }
}

#ifndef LNP22_DIRECT_COEFF_BLOCK
#define LNP22_DIRECT_COEFF_BLOCK ((LOGQ == 72) ? 32 : 0)
#endif
static void compute_direct_public_chain_coeffs(selector_lnp22_ctx *ctx) {
#if defined(_OPENMP) && LNP22_DIRECT_COEFF_PARALLEL
  const size_t transitions = ctx->Lchain-1;
  const size_t digits_per_row = ctx->gdigits/ctx->m;

  /* Each transition owns D[t] and the negative contribution to B[t+1].
     Add the G contribution only after the barrier, preserving the serial
     per-coefficient accumulation order without adjacent-layer races. */
#pragma omp parallel
  {
    polx factor[1];
#if LNP22_DIRECT_COEFF_BLOCK > 0
    polx products[LNP22_DIRECT_COEFF_BLOCK];
#else
    polx tmp[1];
#endif
#pragma omp for schedule(static)
    for(size_t t=0;t<ctx->Lchain;t++) {
      polxvec_setzero(&ctx->direct_lin_B[t*ctx->gdigits],ctx->gdigits);
      if(t<transitions)
        polxvec_setzero(&ctx->directD[t*ctx->gdigits],ctx->gdigits);
    }
#pragma omp for schedule(static)
    for(size_t t=0;t<transitions;t++) {
      for(size_t row=0;row<ctx->m;row++) {
        polx_mul(factor,&ctx->directfold[t],&ctx->rhovec[t*ctx->m + row]);
#if LNP22_DIRECT_COEFF_BLOCK > 0
        for(size_t off=0;off<ctx->gdigits;off+=LNP22_DIRECT_COEFF_BLOCK) {
          const size_t count=MIN(ctx->gdigits-off,LNP22_DIRECT_COEFF_BLOCK);
          polx *dst=&ctx->direct_lin_B[(t+1)*ctx->gdigits+off];
          polxvec_polx_mul(products,factor,&ctx->A0[row*ctx->gdigits+off],count);
          polxvec_sub(dst,dst,products,count);
          dst=&ctx->directD[t*ctx->gdigits+off];
          polxvec_polx_mul(products,factor,&ctx->Dsel[row*ctx->gdigits+off],count);
          polxvec_add(dst,dst,products,count);
        }
#else
        for(size_t i=0;i<ctx->gdigits;i++) {
          polx_mul(tmp,factor,&ctx->A0[row*ctx->gdigits + i]);
          polx_sub(&ctx->direct_lin_B[(t+1)*ctx->gdigits + i],
                   &ctx->direct_lin_B[(t+1)*ctx->gdigits + i],tmp);
          polx_mul(tmp,factor,&ctx->Dsel[row*ctx->gdigits + i]);
          polx_add(&ctx->directD[t*ctx->gdigits + i],
                   &ctx->directD[t*ctx->gdigits + i],tmp);
        }
#endif
      }
    }
#pragma omp for schedule(static)
    for(size_t t=0;t<transitions;t++) {
      for(size_t row=0;row<ctx->m;row++) {
        polx_mul(factor,&ctx->directfold[t],&ctx->rhovec[t*ctx->m + row]);
        for(size_t j=0;j<digits_per_row;j++)
          polx_scale_add(&ctx->direct_lin_B[t*ctx->gdigits + row*digits_per_row + j],
                         factor,pow3_u64(j));
      }
    }
#pragma omp for schedule(static)
    for(size_t i=0;i<transitions*ctx->gdigits;i++)
      polx_refresh(&ctx->directD[i]);
#pragma omp for schedule(static)
    for(size_t i=0;i<ctx->b_len;i++)
      polx_refresh(&ctx->direct_lin_B[i]);
  }
#else
  size_t t, row, i, j;
  polx factor[1], tmp[1];

  polxvec_setzero(ctx->directD,(ctx->Lchain-1)*ctx->gdigits);
  polxvec_setzero(ctx->direct_lin_B,ctx->b_len);
  for(t=0;t<ctx->Lchain-1;t++) {
    for(row=0;row<ctx->m;row++) {
      polx_mul(factor,&ctx->directfold[t],&ctx->rhovec[t*ctx->m + row]);
      for(j=0;j<ctx->gdigits/ctx->m;j++)
        polx_scale_add(&ctx->direct_lin_B[t*ctx->gdigits + row*(ctx->gdigits/ctx->m) + j],
                       factor,pow3_u64(j));
      for(i=0;i<ctx->gdigits;i++) {
        polx_mul(tmp,factor,&ctx->A0[row*ctx->gdigits + i]);
        polx_sub(&ctx->direct_lin_B[(t+1)*ctx->gdigits + i],
                 &ctx->direct_lin_B[(t+1)*ctx->gdigits + i],tmp);
        polx_mul(tmp,factor,&ctx->Dsel[row*ctx->gdigits + i]);
        polx_add(&ctx->directD[t*ctx->gdigits + i],
                 &ctx->directD[t*ctx->gdigits + i],tmp);
      }
    }
  }
  polxvec_refresh(ctx->directD,(ctx->Lchain-1)*ctx->gdigits);
  polxvec_refresh(ctx->direct_lin_B,ctx->b_len);
#endif
}

static void compute_direct_chain_coeffs(selector_lnp22_ctx *ctx) {
  size_t t, i;
  polx tmp[1], dvalue[1], one[1];
  const polx *Xs = &ctx->s1[ctx->off[LNP22_X_BLOCK]];
  const polx *RBITs = &ctx->s1[ctx->off[LNP22_RHAT_BITS_BLOCK]];

  compute_direct_public_chain_coeffs(ctx);

  polxvec_setzero(ctx->g1coeff,ctx->m1);
  polxvec_copy(&ctx->g1coeff[ctx->off[LNP22_B_BLOCK]],ctx->direct_lin_B,ctx->b_len);
  for(t=0;t<ctx->Lchain-1;t++) {
    const size_t xpos = selector_x_pos(ctx,t);
    const polx *x = &Xs[xpos];
    polxvec_sprod(dvalue,&ctx->directD[t*ctx->gdigits],
                  &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],ctx->gdigits);
    polx_sub(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + xpos],
             &ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + xpos],dvalue);
    for(i=0;i<ctx->gdigits;i++) {
      polx_mul(tmp,&ctx->directD[t*ctx->gdigits + i],x);
      polx_sub(&ctx->g1coeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],
               &ctx->g1coeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits + i],tmp);
    }
  }

  polx_monomial(one,1,0);
  if(LNP22_PROVE_X_BITS) {
    for(i=0;i<ctx->x_len;i++) {
      polx_scale(tmp,&Xs[i],2);
      polx_sub(tmp,tmp,one);
      polx_mul(tmp,&ctx->xbit_weight[i],tmp);
      polx_add(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],
               &ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],tmp);
    }
  }
  if(LNP22_USE_RBIT_BINARY) {
    for(i=0;i<ctx->rhat_bits_len;i++) {
      polx_scale(tmp,&RBITs[i],2);
      polx_sub(tmp,tmp,one);
      polx_add(&ctx->g1coeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],
               &ctx->g1coeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],tmp);
    }
  }
  if(LNP22_USE_NATIVE_FOLD_BINARY)
    polxvec_add(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK]],
                &ctx->g1coeff[ctx->off[LNP22_X_BLOCK]],
                ctx->native_binary_lin1,ctx->x_len);
  polxvec_refresh(ctx->g1coeff,ctx->m1);
}

static void compute_quad_g1_coeffs(selector_lnp22_ctx *ctx) {
  size_t t, i;
  const polx *Xs = &ctx->s1[ctx->off[LNP22_X_BLOCK]];
  const polx *Y0s = &ctx->s1[ctx->off[LNP22_Y0_BLOCK]];
  const polx *Y1s = &ctx->s1[ctx->off[LNP22_Y1_BLOCK]];
  const polx *RBITs = LNP22_USE_RBIT_BINARY ? &ctx->s1[ctx->off[LNP22_RHAT_BITS_BLOCK]] : NULL;
  polx one[1];

  if(LNP22_DIRECT_CHAIN) {
    compute_direct_chain_coeffs(ctx);
    return;
  }

  polxvec_setzero(ctx->g1coeff,ctx->m1);
  polx_monomial(one,1,0);

  for(t=0;t<ctx->Lchain-1;t++) {
    for(i=0;i<ctx->gdigits;i++) {
      polx factor[1];
      polx_mul(factor,&ctx->alphavec[t],&ctx->kvec[t*ctx->gdigits + i]);
      polx_add(&ctx->g1coeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits + i],
               &ctx->g1coeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits + i],factor);
    }
  }

  for(i=0;i<ctx->x_len;i++) {
    polx tmp[1];
    polx_scale(tmp,&Y1s[i],Y_SPLIT_GAMMA);
    polx_add(tmp,tmp,&Y0s[i]);
    polx_neg(tmp,tmp);
    polx_add(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],
             &ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],tmp);
  }

  for(i=0;i<ctx->y_len;i++) {
    polx_neg(&ctx->g1coeff[ctx->off[LNP22_Y0_BLOCK] + i],&Xs[i]);
    polx_scale(&ctx->g1coeff[ctx->off[LNP22_Y1_BLOCK] + i],&Xs[i],-Y_SPLIT_GAMMA);
  }

  if(LNP22_PROVE_X_BITS) {
    for(i=0;i<ctx->x_len;i++) {
      polx tmp[1];
      polx_scale(tmp,&Xs[i],2);
      polx_sub(tmp,tmp,one);
      polx_add(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],
               &ctx->g1coeff[ctx->off[LNP22_X_BLOCK] + i],tmp);
    }
  }

  if(LNP22_USE_RBIT_BINARY) {
    for(i=0;i<ctx->rhat_bits_len;i++) {
      polx tmp[1];
      polx_scale(tmp,&RBITs[i],2);
      polx_sub(tmp,tmp,one);
      polx_add(&ctx->g1coeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],
               &ctx->g1coeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],tmp);
    }
  }
  if(LNP22_USE_NATIVE_FOLD_BINARY)
    polxvec_add(&ctx->g1coeff[ctx->off[LNP22_X_BLOCK]],
                &ctx->g1coeff[ctx->off[LNP22_X_BLOCK]],
                ctx->native_binary_lin1,ctx->x_len);
}

static void prepare_selector_lnp22(selector_lnp22_ctx *ctx) {
  uint8_t *xbits = _malloc(ctx->Lchain);
  poly *g1poly = alloc_polyvec(ctx->gdigits);
  uint64_t nonce = 0;
  size_t row, t, j;
#if !LNP22_DIRECT_CHAIN
  polx yx[1], scaled[1];
#endif
#if defined(OUTER_PREP_PROFILE) && OUTER_PREP_PROFILE
  struct timespec prep_prof_last, prep_prof_now;
#define PREP_PROFILE_START() clock_gettime(CLOCK_MONOTONIC,&prep_prof_last)
#define PREP_PROFILE_MARK(label) do {                                       \
    clock_gettime(CLOCK_MONOTONIC,&prep_prof_now);                          \
    printf("      prepare stmt: %-24s %.6fs\n",label,                       \
           elapsed_sec(&prep_prof_last,&prep_prof_now));                    \
    prep_prof_last = prep_prof_now;                                          \
  } while(0)
#else
#define PREP_PROFILE_START() do { } while(0)
#define PREP_PROFILE_MARK(label) do { } while(0)
#endif

  PREP_PROFILE_START();
  shake128(xbits,ctx->Lchain,ctx->statement_seed,16);

  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);
  PREP_PROFILE_MARK("selector matrices");

  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits + row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);

  for(t=0;t<ctx->Lchain;t++)
    poly_set_const(&ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,t)],xbits[t] & 1);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_X_BLOCK]],ctx->witness[LNP22_X_BLOCK],ctx->x_len);
  PREP_PROFILE_MARK("top gadget+x");
  precompute_umc_commitments(ctx);
  PREP_PROFILE_MARK("precompute UMC");
  ctx->selected_j = (size_t)(xbits[0] % ctx->preproc_count);
  load_preprocessed_commitment(ctx,ctx->selected_j);
  PREP_PROFILE_MARK("load selected UMC");

  if(ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,ctx->Lchain-1)].vec->c[0]) {
    polyvec_copy(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],g1poly,ctx->gdigits);
    polxvec_copy(&ctx->s1[ctx->off[LNP22_B_BLOCK] + (ctx->Lchain-1)*ctx->gdigits],
                 ctx->g1,ctx->gdigits);
  }
  else
    polxvec_copy(&ctx->s1[ctx->off[LNP22_B_BLOCK] + (ctx->Lchain-1)*ctx->gdigits],
                 ctx->g0,ctx->gdigits);

#if LNP22_DIRECT_CHAIN && defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel private(row)
  {
    for(size_t tt=ctx->Lchain-1;tt-- > 0;) {
#pragma omp for schedule(static)
      for(row=0;row<ctx->m;row++) {
        polx rhsrow[1];

        polxvec_sprod(rhsrow,&ctx->A0[row*ctx->gdigits],
                      &ctx->s1[ctx->off[LNP22_B_BLOCK] + (tt+1)*ctx->gdigits],ctx->gdigits);
        if(ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,tt)].vec->c[0])
          polxvec_sprod_add(rhsrow,&ctx->Dsel[row*ctx->gdigits],
                            &ctx->s1[ctx->off[LNP22_B_BLOCK] + (tt+1)*ctx->gdigits],ctx->gdigits);
        gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][tt*ctx->gdigits + row*(ctx->gdigits/ctx->m)],
                                   rhsrow,ctx->gdigits/ctx->m);
      }
#pragma omp for schedule(static)
      for(j=0;j<ctx->gdigits;j++)
        polx_frompoly(&ctx->s1[ctx->off[LNP22_B_BLOCK] + tt*ctx->gdigits + j],
                      &ctx->witness[LNP22_B_BLOCK][tt*ctx->gdigits + j]);
    }
  }
#elif LNP22_DIRECT_CHAIN
  for(t=ctx->Lchain-1;t-- > 0;) {
    for(row=0;row<ctx->m;row++) {
      polx rhsrow[1];

      polxvec_sprod(rhsrow,&ctx->A0[row*ctx->gdigits],
                    &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],ctx->gdigits);
      if(ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,t)].vec->c[0])
        polxvec_sprod_add(rhsrow,&ctx->Dsel[row*ctx->gdigits],
                          &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],ctx->gdigits);
      gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits + row*(ctx->gdigits/ctx->m)],
                                 rhsrow,ctx->gdigits/ctx->m);
    }
    polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_B_BLOCK] + t*ctx->gdigits],
                        &ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits],ctx->gdigits);
  }
#elif defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel private(row)
  {
    for(size_t tt=ctx->Lchain-1;tt-- > 0;) {
#pragma omp single
      {
        if(ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,tt)].vec->c[0])
          polyvec_copy(&ctx->witness[LNP22_Z_BLOCK][tt*ctx->gdigits],
                       &ctx->witness[LNP22_B_BLOCK][(tt+1)*ctx->gdigits],ctx->gdigits);
        else
          polyvec_setzero(&ctx->witness[LNP22_Z_BLOCK][tt*ctx->gdigits],ctx->gdigits);
      }
#pragma omp for schedule(static)
      for(j=0;j<ctx->gdigits;j++)
        polx_frompoly(&ctx->s1[ctx->off[LNP22_Z_BLOCK] + tt*ctx->gdigits + j],
                      &ctx->witness[LNP22_Z_BLOCK][tt*ctx->gdigits + j]);
#pragma omp for schedule(static)
      for(row=0;row<ctx->m;row++) {
        polx rhsrow[1];

        polxvec_sprod(rhsrow,&ctx->A0[row*ctx->gdigits],
                      &ctx->s1[ctx->off[LNP22_B_BLOCK] + (tt+1)*ctx->gdigits],ctx->gdigits);
        polxvec_sprod_add(rhsrow,&ctx->Dsel[row*ctx->gdigits],
                          &ctx->s1[ctx->off[LNP22_Z_BLOCK] + tt*ctx->gdigits],ctx->gdigits);
        gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][tt*ctx->gdigits + row*(ctx->gdigits/ctx->m)],
                                   rhsrow,ctx->gdigits/ctx->m);
      }
#pragma omp for schedule(static)
      for(j=0;j<ctx->gdigits;j++)
        polx_frompoly(&ctx->s1[ctx->off[LNP22_B_BLOCK] + tt*ctx->gdigits + j],
                      &ctx->witness[LNP22_B_BLOCK][tt*ctx->gdigits + j]);
    }
  }
#else
  for(t=ctx->Lchain-1;t-- > 0;) {
    if(ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,t)].vec->c[0])
      polyvec_copy(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],
                   &ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    else
      polyvec_setzero(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],ctx->gdigits);

    polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits],
                        &ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],ctx->gdigits);
    for(row=0;row<ctx->m;row++) {
      polx rhsrow[1];

      polxvec_sprod(rhsrow,&ctx->A0[row*ctx->gdigits],
                    &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],ctx->gdigits);
      polxvec_sprod_add(rhsrow,&ctx->Dsel[row*ctx->gdigits],
                        &ctx->s1[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits],ctx->gdigits);
      gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits + row*(ctx->gdigits/ctx->m)],
                                 rhsrow,ctx->gdigits/ctx->m);
    }
    polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_B_BLOCK] + t*ctx->gdigits],
                        &ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits],ctx->gdigits);
  }
#endif

  PREP_PROFILE_MARK(LNP22_DIRECT_CHAIN ? "direct chain B" : "chain B/Z");
  PREP_PROFILE_MARK(LNP22_DIRECT_CHAIN ? "B already synced" : "B/Z already synced");
  compute_Cout(ctx);
  PREP_PROFILE_MARK("compute Cout");
  derive_selector_chal_seed_from_public(ctx);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(t=0;t<ctx->Lchain-1;t++) {
    if(!LNP22_DIRECT_CHAIN)
      polxvec_challenge(&ctx->kvec[t*ctx->gdigits],ctx->gdigits,ctx->selector_chal_seed,(uint64_t)t);
    if(!LNP22_DIRECT_CHAIN)
      polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,ctx->selector_chal_seed,32 + (uint64_t)t);
  }
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->selector_chal_seed,64);
  if(!LNP22_DIRECT_CHAIN)
    polxvec_challenge(ctx->alphavec,ctx->Lchain-1,ctx->selector_chal_seed,65);
  PREP_PROFILE_MARK("selector challenges");

#if !LNP22_DIRECT_CHAIN
  for(t=0;t<ctx->Lchain-1;t++) {
    polz yz[1];
    zz coeff;

    polxvec_sprod(yx,&ctx->kvec[t*ctx->gdigits],
                  &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],ctx->gdigits);
    polx_mul(scaled,&ctx->alphavec[t],yx);
    polz_frompolx(yz,scaled);
    polz_center(yz);
    for(j=0;j<N;j++) {
      int64_t c;
      int64_t q;
      int64_t rem;

      for(row=0;row<L;row++)
        coeff.limbs[row] = yz->limbs[row].c[j];
      c = zz_to_i64(&coeff);
      if(c >= 0)
        q = (c + Y_SPLIT_GAMMA/2) / Y_SPLIT_GAMMA;
      else
        q = -(((-c) + Y_SPLIT_GAMMA/2) / Y_SPLIT_GAMMA);
      rem = c - Y_SPLIT_GAMMA*q;
      ctx->witness[LNP22_Y0_BLOCK][t*SMALL_BLOCK].vec->c[j] = (int16_t)rem;
      ctx->witness[LNP22_Y1_BLOCK][t*SMALL_BLOCK].vec->c[j] = (int16_t)q;
    }
  }

  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Y0_BLOCK]],ctx->witness[LNP22_Y0_BLOCK],ctx->y_len);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Y1_BLOCK]],ctx->witness[LNP22_Y1_BLOCK],ctx->y_len);
#endif
  PREP_PROFILE_MARK("Y split+sync");
  if(ctx->cache_a1open)
    fill_almostuniform_matrix(ctx->A1open,ctx->kmsis,ctx->m1,ctx->lnp_seed,NONCE_LNP_A1_BASE);
  fill_almostuniform_matrix(ctx->A2open,ctx->kmsis,ctx->m2prime,ctx->lnp_seed,NONCE_LNP_A2_BASE);
  PREP_PROFILE_MARK("Aopen fill");
  if(!LNP22_DIRECT_CHAIN)
    compute_quad_g1_coeffs(ctx);
  PREP_PROFILE_MARK("quad g1 coeffs");

  free(xbits);
  free(g1poly);
  PREP_PROFILE_MARK("cleanup");
#undef PREP_PROFILE_START
#undef PREP_PROFILE_MARK
}

static void linZ_value(polx *out, const selector_lnp22_ctx *ctx, const polx *zblock) {
  size_t t;
  polx inner[1], scaled[1];
  polxvec_setzero(out,1);
  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_sprod(inner,&ctx->kvec[t*ctx->gdigits],&zblock[t*ctx->gdigits],ctx->gdigits);
    polx_mul(scaled,&ctx->alphavec[t],inner);
    polx_add(out,out,scaled);
  }
}

static void bilinear_value(polx *out, const selector_lnp22_ctx *ctx,
                           const polx *xblock, const polx *y0block, const polx *y1block) {
  polx tmp[1];
  polxvec_sprod(out,xblock,y0block,ctx->x_len);
  polxvec_sprod(tmp,xblock,y1block,ctx->x_len);
  polx_scale(tmp,tmp,Y_SPLIT_GAMMA);
  polx_add(out,out,tmp);
}

static void polxvec_sum(polx *out, const polx *vec, size_t len) {
  size_t i;

  polxvec_setzero(out,1);
  for(i=0;i<len;i++)
    polx_add(out,out,&vec[i]);
}

static void rhat_bits_quad_value(polx *out, const selector_lnp22_ctx *ctx,
                                 const polx *bits) {
  polx sum[1];
  polxvec_sprod(out,bits,bits,ctx->rhat_bits_len);
  polxvec_sum(sum,bits,ctx->rhat_bits_len);
  polx_sub(out,out,sum);
  polx_refresh(out);
}

static void rhat_bits_verifier_value(polx *out, const selector_lnp22_ctx *ctx,
                                     const polx *zbits, const polx *c) {
  polx sum[1];
  polxvec_sprod(out,zbits,zbits,ctx->rhat_bits_len);
  polxvec_sum(sum,zbits,ctx->rhat_bits_len);
  polx_mul(sum,sum,c);
  polx_sub(out,out,sum);
  polx_refresh(out);
}

static void x_bits_quad_value(polx *out, const selector_lnp22_ctx *ctx,
                              const polx *xblock) {
  polx lin[1];
  size_t i;

  if(LNP22_DIRECT_CHAIN) {
    polx term[1], weighted[1], one[1];

    polx_monomial(one,1,0);
    polxvec_setzero(out,1);
    for(i=0;i<ctx->x_len;i++) {
      polx_sub(term,&xblock[i],one);
      polx_mul(term,term,&xblock[i]);
      polx_mul(weighted,&ctx->xbit_weight[i],term);
      polx_add(out,out,weighted);
    }
    polx_refresh(out);
    return;
  }

  polxvec_sprod(out,xblock,xblock,ctx->x_len);
  polxvec_sum(lin,xblock,ctx->x_len);
  polx_sub(out,out,lin);
}

static void x_bits_verifier_value(polx *out, const selector_lnp22_ctx *ctx,
                                  const polx *zblock, const polx *c) {
  polx lin[1];
  size_t i;

  if(LNP22_DIRECT_CHAIN) {
    polx term[1], weighted[1];

    polxvec_setzero(out,1);
    for(i=0;i<ctx->x_len;i++) {
      polx_mul(term,&zblock[i],&zblock[i]);
      polx_mul(weighted,c,&zblock[i]);
      polx_sub(term,term,weighted);
      polx_mul(weighted,&ctx->xbit_weight[i],term);
      polx_add(out,out,weighted);
    }
    polx_refresh(out);
    return;
  }

  polxvec_sprod(out,zblock,zblock,ctx->x_len);
  polxvec_sum(lin,zblock,ctx->x_len);
  polx_mul(lin,c,lin);
  polx_sub(out,out,lin);
}

#if LNP22_SOUNDNESS_SELFTEST
static int run_lnp22_soundness_selftest(selector_lnp22_ctx *ctx,
                                        const polx *tA) {
  const size_t rho_len = (ctx->Lchain-1)*ctx->m;
  polx *saved_rho = alloc_polxvec(rho_len);
  polx *changed_tA = alloc_polxvec(ctx->kmsis);
  polx *bad_x = alloc_polxvec(ctx->x_len);
  polx *bad_rbit = LNP22_USE_RBIT_BINARY ?
                   alloc_polxvec(ctx->rhat_bits_len) : NULL;
  polx one[1], diff[1], bit_value[1], saved_h[1];
  size_t i;
  int rho_changed = 0;
  int ok = 1;

  if(ctx->rbit_const_lambda != 0) {
    polxvec_copy(saved_h,&ctx->rbit_const_h[0],1);
    polx_monomial(one,1,0);
    polx_add(&ctx->rbit_const_h[0],&ctx->rbit_const_h[0],one);
    if(!validate_rbit_const_transcript(ctx)) {
      fprintf(stderr,"ERROR: tampered h constant coefficient was accepted\n");
      ok = 0;
    }
    polxvec_copy(&ctx->rbit_const_h[0],saved_h,1);
  }

  polxvec_copy(saved_rho,ctx->rhovec,rho_len);
  polxvec_copy(changed_tA,tA,ctx->kmsis);
  polx_monomial(one,1,0);
  polx_add(&changed_tA[0],&changed_tA[0],one);
  derive_direct_chain_fold(ctx,changed_tA);
  for(i=0;i<rho_len;i++) {
    polx_sub(diff,&ctx->rhovec[i],&saved_rho[i]);
    if(!polx_iszero(diff)) {
      rho_changed = 1;
      break;
    }
  }
  if(!rho_changed) {
    fprintf(stderr,"ERROR: transition rho did not change after t_A tampering\n");
    ok = 0;
  }
  derive_direct_chain_fold(ctx,tA);

  polxvec_setzero(bad_x,ctx->x_len);
  polx_monomial(&bad_x[0],2,0);
  x_bits_quad_value(bit_value,ctx,bad_x);
  if(polx_iszero(bit_value)) {
    fprintf(stderr,"ERROR: weighted bit relation accepted x[0] = 2\n");
    ok = 0;
  }
  if(LNP22_USE_RBIT_BINARY) {
    polxvec_setzero(bad_rbit,ctx->rhat_bits_len);
    polx_monomial(&bad_rbit[0],2,0);
    rhat_bits_quad_value(bit_value,ctx,bad_rbit);
    if(polx_iszero(bit_value)) {
      fprintf(stderr,"ERROR: weighted RBIT relation accepted RBIT[0] = 2\n");
      ok = 0;
    }
  }

  free(saved_rho);
  free(changed_tA);
  free(bad_x);
  free(bad_rbit);
  printf("LNP22 soundness self-test: %s\n",ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
#endif

static void polx_allones(polx *out) {
  poly p[1];
  size_t j;

  poly_set_zero(p);
  for(j=0;j<N;j++)
    p->vec->c[j] = 1;
  polx_frompoly(out,p);
}

static void rhat_native_binary_value(polx *out, const selector_lnp22_ctx *ctx,
                                     const polx *rhat) {
  polx ones[1], aut[1], diff[1], term[1];
  size_t i;

  polx_allones(ones);
  polxvec_setzero(out,1);
  for(i=0;i<ctx->rhat_len;i++) {
    polx_sub(diff,&rhat[i],ones);
    polx_sigmam1(aut,&rhat[i]);
    polx_mul(term,diff,aut);
    polx_add(out,out,term);
  }
  polx_refresh(out);
}

static void rhat_native_binary_verifier_value(
    polx *out, const selector_lnp22_ctx *ctx,
    const polx *rhat_z, const polx *c) {
  polx ones[1], aut_z[1], term[1], linear[1];
  size_t i;

  polx_allones(ones);
  polxvec_setzero(out,1);
  for(i=0;i<ctx->rhat_len;i++) {
    polx_sigmam1(aut_z,&rhat_z[i]);
    polx_mul(term,&rhat_z[i],aut_z);
    polx_mul(linear,ones,aut_z);
    polx_mul(linear,c,linear);
    polx_sub(term,term,linear);
    polx_add(out,out,term);
  }
  polx_mul(out,ctx->native_binary_scale,out);
}

static void rhat_native_binary_mask_terms(
    polx *g1, polx *g0, const selector_lnp22_ctx *ctx,
    const polx *rhat_s, const polx *rhat_y) {
  polx ones[1], aut_s[1], aut_y[1], diff[1];
  polx cross[1], term[1], quadratic[1];
  size_t i;

  if(!LNP22_USE_NATIVE_FOLD_BINARY)
    return;

  polx_allones(ones);
  polxvec_setzero(cross,1);
  polxvec_setzero(quadratic,1);
  for(i=0;i<ctx->rhat_len;i++) {
    polx_sigmam1(aut_s,&rhat_s[i]);
    polx_sigmam1(aut_y,&rhat_y[i]);
    polx_mul(term,&rhat_y[i],aut_s);
    polx_add(cross,cross,term);
    polx_sub(diff,&rhat_s[i],ones);
    polx_mul(term,diff,aut_y);
    polx_add(cross,cross,term);
    polx_mul(term,&rhat_y[i],aut_y);
    polx_add(quadratic,quadratic,term);
  }
  polx_mul(cross,ctx->native_binary_scale,cross);
  polx_mul(quadratic,ctx->native_binary_scale,quadratic);
  polx_add(g1,g1,cross);
  polx_add(g0,g0,quadratic);
}

static void add_native_binary_s2_mask(polx *g1,
                                      const selector_lnp22_ctx *ctx,
                                      const polx *y2) {
  polx term[1];

  if(!LNP22_USE_NATIVE_FOLD_BINARY)
    return;
  polxvec_sprod(term,ctx->native_binary_lin2,y2,ctx->m2);
  polx_add(g1,g1,term);
}

static void quad_value(polx *out, const selector_lnp22_ctx *ctx, const polx *vec) {
  polx lin[1], bil[1], bits[1];
  const polx *xblock = &vec[ctx->off[LNP22_X_BLOCK]];

  if(LNP22_DIRECT_CHAIN) {
    const polx *bblock = &vec[ctx->off[LNP22_B_BLOCK]];
    size_t t;

    polxvec_sprod(out,ctx->direct_lin_B,bblock,ctx->b_len);
    for(t=0;t<ctx->Lchain-1;t++) {
      polxvec_sprod(lin,&ctx->directD[t*ctx->gdigits],
                    &bblock[(t+1)*ctx->gdigits],ctx->gdigits);
      polx_refresh(lin);
      polx_mul(lin,lin,&xblock[selector_x_pos(ctx,t)]);
      polx_sub(out,out,lin);
    }
  } else {
    const polx *zblock = &vec[ctx->off[LNP22_Z_BLOCK]];
    const polx *y0block = &vec[ctx->off[LNP22_Y0_BLOCK]];
    const polx *y1block = &vec[ctx->off[LNP22_Y1_BLOCK]];

    linZ_value(lin,ctx,zblock);
    bilinear_value(bil,ctx,xblock,y0block,y1block);
    polx_sub(out,lin,bil);
  }
  if(LNP22_PROVE_X_BITS) {
    x_bits_quad_value(bits,ctx,xblock);
    polx_add(out,out,bits);
  }
  if(LNP22_USE_RBIT_BINARY) {
    const polx *rbitblock = &vec[ctx->off[LNP22_RHAT_BITS_BLOCK]];
    rhat_bits_quad_value(bits,ctx,rbitblock);
    polx_add(out,out,bits);
  }
  if(LNP22_USE_NATIVE_FOLD_BINARY) {
    const polx *rhatblock = &vec[ctx->off[LNP22_RHAT_BLOCK]];

    rhat_native_binary_value(bits,ctx,rhatblock);
    polx_mul(bits,ctx->native_binary_scale,bits);
    polx_add(out,out,bits);
    polxvec_sprod(bits,ctx->native_binary_lin1,xblock,ctx->x_len);
    polx_add(out,out,bits);
  }
}

static void quad_mask_terms(polx *g1, polx *g0, const selector_lnp22_ctx *ctx,
                            const polx *s1, const polx *y1) {
  const polx *Xs = &s1[ctx->off[LNP22_X_BLOCK]];
  const polx *Y0s = &s1[ctx->off[LNP22_Y0_BLOCK]];
  const polx *Y1s = &s1[ctx->off[LNP22_Y1_BLOCK]];
  const polx *Zy = &y1[ctx->off[LNP22_Z_BLOCK]];
  const polx *Xy = &y1[ctx->off[LNP22_X_BLOCK]];
  const polx *Y0y = &y1[ctx->off[LNP22_Y0_BLOCK]];
  const polx *Y1y = &y1[ctx->off[LNP22_Y1_BLOCK]];
  const polx *RBITs = &s1[ctx->off[LNP22_RHAT_BITS_BLOCK]];
  const polx *RBITy = &y1[ctx->off[LNP22_RHAT_BITS_BLOCK]];
  const polx *RHATs = &s1[ctx->off[LNP22_RHAT_BLOCK]];
  const polx *RHATy = &y1[ctx->off[LNP22_RHAT_BLOCK]];
  polx tmp0[1], one[1];
  if(ctx->g1coeff != NULL) {
    polxvec_sprod(g1,ctx->g1coeff,y1,ctx->m1);
  } else {
    linZ_value(g1,ctx,Zy);
    polxvec_sprod(tmp0,Xy,Y0s,ctx->x_len);
    polx_sub(g1,g1,tmp0);
    polxvec_sprod(tmp0,Xs,Y0y,ctx->x_len);
    polx_sub(g1,g1,tmp0);
    polxvec_sprod(tmp0,Xy,Y1s,ctx->x_len);
    polx_scale(tmp0,tmp0,Y_SPLIT_GAMMA);
    polx_sub(g1,g1,tmp0);
    polxvec_sprod(tmp0,Xs,Y1y,ctx->x_len);
    polx_scale(tmp0,tmp0,Y_SPLIT_GAMMA);
    polx_sub(g1,g1,tmp0);
    if(LNP22_USE_RBIT_BINARY) {
      polx_monomial(one,1,0);
      for(size_t i=0;i<ctx->rhat_bits_len;i++) {
        polx coeff[1];
        polx_scale(coeff,&RBITs[i],2);
        polx_sub(coeff,coeff,one);
        polx_mul(tmp0,coeff,&RBITy[i]);
        polx_add(g1,g1,tmp0);
      }
    }
    if(LNP22_PROVE_X_BITS) {
      polx_monomial(one,1,0);
      for(size_t i=0;i<ctx->x_len;i++) {
        polx coeff[1];
        polx_scale(coeff,&Xs[i],2);
        polx_sub(coeff,coeff,one);
        polx_mul(tmp0,coeff,&Xy[i]);
        polx_add(g1,g1,tmp0);
      }
    }
  }

  if(LNP22_DIRECT_CHAIN) {
    const polx *By = &y1[ctx->off[LNP22_B_BLOCK]];
    size_t t;

    polxvec_setzero(g0,1);
    for(t=0;t<ctx->Lchain-1;t++) {
      polxvec_sprod(tmp0,&ctx->directD[t*ctx->gdigits],
                    &By[(t+1)*ctx->gdigits],ctx->gdigits);
      polx_refresh(tmp0);
      polx_mul(tmp0,tmp0,&Xy[selector_x_pos(ctx,t)]);
      polx_sub(g0,g0,tmp0);
    }
    if(LNP22_PROVE_X_BITS) {
      for(size_t i=0;i<ctx->x_len;i++) {
        polx_mul(tmp0,&Xy[i],&Xy[i]);
        polx_mul(tmp0,&ctx->xbit_weight[i],tmp0);
        polx_add(g0,g0,tmp0);
      }
      polx_refresh(g0);
    }
    if(LNP22_USE_RBIT_BINARY) {
      polxvec_sprod(tmp0,RBITy,RBITy,ctx->rhat_bits_len);
      polx_add(g0,g0,tmp0);
      polx_refresh(g0);
    }
    rhat_native_binary_mask_terms(g1,g0,ctx,RHATs,RHATy);
    return;
  }

  polxvec_sprod(g0,Xy,Y0y,ctx->x_len);
  polx_neg(g0,g0);
  polxvec_sprod(tmp0,Xy,Y1y,ctx->x_len);
  polx_scale(tmp0,tmp0,Y_SPLIT_GAMMA);
  polx_sub(g0,g0,tmp0);
  if(LNP22_PROVE_X_BITS) {
    polxvec_sprod(tmp0,Xy,Xy,ctx->x_len);
    polx_add(g0,g0,tmp0);
  }
  if(LNP22_USE_RBIT_BINARY) {
    polxvec_sprod(tmp0,RBITy,RBITy,ctx->rhat_bits_len);
    polx_add(g0,g0,tmp0);
  }
  rhat_native_binary_mask_terms(g1,g0,ctx,RHATs,RHATy);
}

static void folded_linear_eval(polx *out, const selector_lnp22_ctx *ctx,
                               const polx *vec1, const polx *vec2) {
  polx tmp[1];

  polxvec_sprod(out,ctx->lincoeff,vec1,ctx->m1);
  if(ctx->lincoeff2 != NULL && vec2 != NULL) {
    polxvec_sprod(tmp,ctx->lincoeff2,vec2,ctx->m2);
    polx_add(out,out,tmp);
  }
}

static int first_linear_relation_failure(const selector_lnp22_ctx *ctx, const polx *vec) {
  size_t row, t, j;
  const polx *B = &vec[ctx->off[LNP22_B_BLOCK]];
  const polx *X = &vec[ctx->off[LNP22_X_BLOCK]];
  const polx *Z = &vec[ctx->off[LNP22_Z_BLOCK]];
  const polx *Y0 = &vec[ctx->off[LNP22_Y0_BLOCK]];
  const polx *Y1 = &vec[ctx->off[LNP22_Y1_BLOCK]];
  const polx *RHAT = &vec[ctx->off[LNP22_RHAT_BLOCK]];
  const polx *EUA = &vec[ctx->off[LNP22_UMC_EUA_BLOCK]];
  const polx *EUB = &vec[ctx->off[LNP22_UMC_EUB_BLOCK]];
  const polx *R = &vec[ctx->off[LNP22_UMC_R_BLOCK]];
  const polx *VA = &vec[ctx->off[LNP22_UMC_VA_BLOCK]];
  const polx *VB = &vec[ctx->off[LNP22_UMC_VB_BLOCK]];
  const polx *RBITS = LNP22_USE_RBIT_BINARY ? &vec[ctx->off[LNP22_RHAT_BITS_BLOCK]] : NULL;
  polx rel[1], tmp[1], c0[1], cdelta[1], gx[1], gy[1], term[1];

  polxvec_sprod(c0,ctx->kprime,ctx->g0,ctx->gdigits);
  polxvec_sprod(cdelta,ctx->kprime,ctx->gdelta,ctx->gdigits);

  for(row=0;row<ctx->umc_n;row++) {
    polxvec_sprod(rel,&ctx->Aumc[row*ctx->umc_k],R,ctx->umc_k);
    polx_add(rel,rel,&EUA[row]);
    polx_scale(tmp,&VA[row],UMC_MODP);
    polx_sub(rel,rel,tmp);
    polx_sub(rel,rel,&ctx->CX[row]);
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in UMC c_A row %zu\n",row);
      return 1;
    }
  }

  for(t=0;t<ctx->umc_v;t++) {
    const polx *msg = &RHAT[t];
    polxvec_sprod(rel,&ctx->Bumc[t*ctx->umc_k],R,ctx->umc_k);
    polx_add(rel,rel,&EUB[t]);
    polx_scale(tmp,msg,UMC_MSG_SCALE);
    polx_add(rel,rel,tmp);
    polx_scale(tmp,&VB[t],UMC_MODP);
    polx_sub(rel,rel,tmp);
    polx_sub(rel,rel,&ctx->CX[ctx->umc_n + t]);
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in UMC c_B row %zu\n",t);
      return 1;
    }
  }

  if(LNP22_USE_RBIT_BINARY) {
    for(row=0;row<ctx->rhat_len;row++) {
      polxvec_copy(rel,&RHAT[row],1);
      for(j=0;j<N;j++) {
        polx_monomial(term,1,(int)j);
        polx_mul(term,term,&RBITS[row*N + j]);
        polx_sub(rel,rel,term);
      }
      if(!polx_iszero(rel)) {
        fprintf(stderr,"DEBUG: linear fail in RHAT bit packing row %zu\n",row);
        return 4;
      }
    }
  }

  if(LNP22_USE_RBIT_BINARY) for(row=0;row<umc_range_rows(ctx);row++) {
    umc_range_lhs(rel,ctx,vec,row);
    polx_sub(rel,rel,&ctx->ulin[ctx->cx_rows+ctx->rhat_len+row]);
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in UMC range packing row %zu\n",row);
      return 9;
    }
  }

  for(row=0;row<ctx->m;row++) {
    polxvec_sprod(rel,&ctx->Arhash[row*ctx->rhat_len],RHAT,ctx->rhat_len);
    polxvec_setzero(gx,1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(gx,&B[row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
    polx_add(rel,rel,gx);
    polx_sub(rel,rel,&ctx->Cout[row]);
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in C well-formedness row %zu\n",row);
      return 5;
    }
  }

  polxvec_sprod(rel,ctx->kprime,&B[(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polx_mul(tmp,cdelta,&X[selector_x_pos(ctx,ctx->Lchain-1)]);
  polx_sub(rel,rel,tmp);
  polx_sub(rel,rel,c0);
  if(!polx_iszero(rel)) {
    fprintf(stderr,"DEBUG: linear fail in top selector\n");
    return 6;
  }

  if(!LNP22_DIRECT_CHAIN) for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_sprod(rel,&ctx->kvec[t*ctx->gdigits],&B[(t+1)*ctx->gdigits],ctx->gdigits);
    polx_mul(rel,&ctx->alphavec[t],rel);
    polx_scale(gy,&Y1[t*SMALL_BLOCK],Y_SPLIT_GAMMA);
    polx_add(gx,&Y0[t*SMALL_BLOCK],gy);
    polx_sub(rel,gx,rel);
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in Y split relation t=%zu\n",t);
      return 7;
    }
  }

  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_setzero(rel,1);
    for(row=0;row<ctx->m;row++) {
      polx rhs[1], gb[1], dpart[1];

      polxvec_setzero(gb,1);
      for(j=0;j<ctx->gdigits/ctx->m;j++)
        polx_scale_add(gb,&B[t*ctx->gdigits + row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
      polxvec_sprod(rhs,&ctx->A0[row*ctx->gdigits],&B[(t+1)*ctx->gdigits],ctx->gdigits);
      if(LNP22_DIRECT_CHAIN) {
        polxvec_sprod(dpart,&ctx->Dsel[row*ctx->gdigits],
                      &B[(t+1)*ctx->gdigits],ctx->gdigits);
        polx_mul(dpart,dpart,&X[selector_x_pos(ctx,t)]);
      } else {
        polxvec_sprod(dpart,&ctx->Dsel[row*ctx->gdigits],&Z[t*ctx->gdigits],ctx->gdigits);
      }
      polx_add(rhs,rhs,dpart);
      polx_sub(gb,gb,rhs);
      polx_mul(gb,gb,&ctx->rhovec[t*ctx->m + row]);
      polx_add(rel,rel,gb);
    }
    if(!polx_iszero(rel)) {
      fprintf(stderr,"DEBUG: linear fail in transition t=%zu\n",t);
      return 8;
    }
  }

  return 0;
}

static void folded_linear_eval_explicit(polx *out, const selector_lnp22_ctx *ctx, const polx *vec) {
  size_t row, t, j;
  const polx *B = &vec[ctx->off[LNP22_B_BLOCK]];
  const polx *X = &vec[ctx->off[LNP22_X_BLOCK]];
  const polx *Z = &vec[ctx->off[LNP22_Z_BLOCK]];
  const polx *Y0 = &vec[ctx->off[LNP22_Y0_BLOCK]];
  const polx *Y1 = &vec[ctx->off[LNP22_Y1_BLOCK]];
  const polx *RHAT = &vec[ctx->off[LNP22_RHAT_BLOCK]];
  const polx *EUA = &vec[ctx->off[LNP22_UMC_EUA_BLOCK]];
  const polx *EUB = &vec[ctx->off[LNP22_UMC_EUB_BLOCK]];
  const polx *R = &vec[ctx->off[LNP22_UMC_R_BLOCK]];
  const polx *VA = &vec[ctx->off[LNP22_UMC_VA_BLOCK]];
  const polx *VB = &vec[ctx->off[LNP22_UMC_VB_BLOCK]];
  const polx *RBITS = LNP22_USE_RBIT_BINARY ? &vec[ctx->off[LNP22_RHAT_BITS_BLOCK]] : NULL;
  polx rel[1], tmp[1], cdelta[1], gx[1], gy[1], term[1];
  size_t idx = 0;

  polxvec_setzero(out,1);
  polxvec_sprod(cdelta,ctx->kprime,ctx->gdelta,ctx->gdigits);

  for(row=0;row<ctx->umc_n;row++) {
    polxvec_sprod(rel,&ctx->Aumc[row*ctx->umc_k],R,ctx->umc_k);
    polx_add(rel,rel,&EUA[row]);
    polx_scale(tmp,&VA[row],UMC_MODP);
    polx_sub(rel,rel,tmp);
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }

  for(t=0;t<ctx->umc_v;t++) {
    const polx *msg = &RHAT[t];
    polxvec_sprod(rel,&ctx->Bumc[t*ctx->umc_k],R,ctx->umc_k);
    polx_add(rel,rel,&EUB[t]);
    polx_scale(tmp,msg,UMC_MSG_SCALE);
    polx_add(rel,rel,tmp);
    polx_scale(tmp,&VB[t],UMC_MODP);
    polx_sub(rel,rel,tmp);
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }

  if(LNP22_USE_RBIT_BINARY) {
    for(row=0;row<ctx->rhat_len;row++) {
      polxvec_copy(rel,&RHAT[row],1);
      for(j=0;j<N;j++) {
        polx_monomial(term,1,(int)j);
        polx_mul(term,term,&RBITS[row*N + j]);
        polx_sub(rel,rel,term);
      }
      polx_mul(tmp,&ctx->linfold[idx++],rel);
      polx_add(out,out,tmp);
    }
  }

  if(LNP22_USE_RBIT_BINARY) for(row=0;row<umc_range_rows(ctx);row++) {
    umc_range_lhs(rel,ctx,vec,row);
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }

  for(row=0;row<ctx->m;row++) {
    polxvec_sprod(rel,&ctx->Arhash[row*ctx->rhat_len],RHAT,ctx->rhat_len);
    polxvec_setzero(gx,1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(gx,&B[row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
    polx_add(rel,rel,gx);
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }

  polxvec_sprod(rel,ctx->kprime,&B[(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polx_mul(tmp,cdelta,&X[selector_x_pos(ctx,ctx->Lchain-1)]);
  polx_sub(rel,rel,tmp);
  polx_mul(tmp,&ctx->linfold[idx++],rel);
  polx_add(out,out,tmp);

  if(!LNP22_DIRECT_CHAIN) for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_sprod(rel,&ctx->kvec[t*ctx->gdigits],&B[(t+1)*ctx->gdigits],ctx->gdigits);
    polx_mul(rel,&ctx->alphavec[t],rel);
    polx_scale(gy,&Y1[t*SMALL_BLOCK],Y_SPLIT_GAMMA);
    polx_add(gx,&Y0[t*SMALL_BLOCK],gy);
    polx_sub(rel,gx,rel);
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }

  if(!LNP22_DIRECT_CHAIN) for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_setzero(rel,1);
    for(row=0;row<ctx->m;row++) {
      polx rhs[1], gb[1], dpart[1];

      polxvec_setzero(gb,1);
      for(j=0;j<ctx->gdigits/ctx->m;j++)
        polx_scale_add(gb,&B[t*ctx->gdigits + row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
      polxvec_sprod(rhs,&ctx->A0[row*ctx->gdigits],&B[(t+1)*ctx->gdigits],ctx->gdigits);
      polxvec_sprod(dpart,&ctx->Dsel[row*ctx->gdigits],&Z[t*ctx->gdigits],ctx->gdigits);
      polx_add(rhs,rhs,dpart);
      polx_sub(gb,gb,rhs);
      polx_mul(gb,gb,&ctx->rhovec[t*ctx->m + row]);
      polx_add(rel,rel,gb);
    }
    polx_mul(tmp,&ctx->linfold[idx++],rel);
    polx_add(out,out,tmp);
  }
}

static void inner_product_a3(polx *out, const polx *a3, const polx *vec, size_t len) {
  polxvec_sprod(out,a3,vec,len);
}

static void print_size_breakdown(const selector_lnp22_ctx *ctx) {
  const double ring_kb = (double)(N*QBYTES) / 1024.0;
  const double challenge_kb = ((double)N * ceil(log2(2*LNP22_OMEGA + 1))) / 8192.0;
  const double cx_coeff_bits = ceil(log2((double)UMC_MODP));
  const double cx_ring_kb = ((double)N * cx_coeff_bits) / 8192.0;
  const double constcheck_kb = 2.0 * (double)ctx->rbit_const_lambda * ring_kb;

  printf("Inner LNP22 transcript size estimate (uncompressed, explicit w/v)\n");
  printf("  t_A KB          = %.3f\n",ctx->kmsis * ring_kb);
  printf("  w KB            = %.3f\n",ctx->kmsis * ring_kb);
  printf("  tg/h KB         = %.3f\n",constcheck_kb);
  printf("  v_lin KB        = %.3f\n",ring_kb);
  printf("  v_q KB          = %.3f\n",ring_kb);
  printf("  t_q KB          = %.3f\n",ring_kb);
  printf("  challenge KB    = %.6f\n",challenge_kb);
  printf("  z1 KB           = %.3f\n",ctx->m1 * ring_kb);
  printf("  z2 KB           = %.3f\n",ctx->m2 * ring_kb);
  printf("  non-z KB        = %.3f\n",(2*ctx->kmsis + 3) * ring_kb + constcheck_kb + challenge_kb);
  printf("  total KB        = %.3f\n",((double)(ctx->m1 + ctx->m2 + 2*ctx->kmsis + 3))*ring_kb + constcheck_kb + challenge_kb);
  printf("  C_X KB          = %.3f\n",ctx->cx_rows * cx_ring_kb);
  printf("  C KB            = %.3f\n",ctx->m * ring_kb);
}

static void polxvec_centered_to_i64(int64_t *out, const polx *v, size_t len) {
  const size_t chunks = (len + 15)/16;
  size_t blk;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(len >= LNP22_PARALLEL_MIN_LEN)
#endif
  for(blk=0;blk<chunks;blk++) {
    const size_t off_poly = 16*blk;
    const size_t chunk = MIN((size_t)16,len - off_poly);
    polz tmp[16];
    size_t i, j, l;

    polzvec_frompolxvec(tmp,&v[off_poly],chunk);
    polzvec_center(tmp,chunk);
    for(i=0;i<chunk;i++) {
      for(j=0;j<N;j++) {
        zz coeff;
        for(l=0;l<L;l++)
          coeff.limbs[l] = tmp[i].limbs[l].c[j];
        out[(off_poly + i)*N + j] = zz_to_i64(&coeff);
      }
    }
  }
}

static void polxvec_centered_to_i128(__int128 *out, const polx *v, size_t len) {
#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(len >= LNP22_PARALLEL_MIN_LEN)
#endif
  for(size_t i=0;i<len;i++) {
    polz z;
    polz_frompolx(&z,&v[i]);
    polz_center(&z);
    for(size_t j=0;j<N;j++) {
      zz c;
      polz_getcoeff(&c,&z,j);
      out[i*N+j] = zz_toint128(&c);
    }
  }
}

static void witness_coeffs_to_i64(int64_t *out, const selector_lnp22_ctx *ctx) {
  size_t block;

  for(block=0;block<LNP22_NUM_BLOCKS;block++) {
    size_t i;
    for(i=0;i<ctx->n[block];i++) {
      size_t j;
      const poly *src = &ctx->witness[block][i];
      int64_t *dst = &out[(ctx->off[block] + i)*N];
      for(j=0;j<N;j++)
        dst[j] = src->vec->c[j];
    }
  }
}

static inline void coefficient_mul_add(int64_t *dst, const int64_t *src,
                                       size_t len, int64_t scalar, int subtract) {
  size_t i = 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__)
  const __m512i scale = _mm512_set1_epi64(scalar);
  for(;i + 8 <= len;i += 8) {
    const __m512i d = _mm512_loadu_si512((const void *)&dst[i]);
    const __m512i s = _mm512_loadu_si512((const void *)&src[i]);
    const __m512i p = _mm512_mullo_epi64(scale,s);
    const __m512i r = subtract ? _mm512_sub_epi64(d,p) : _mm512_add_epi64(d,p);
    _mm512_storeu_si512((void *)&dst[i],r);
  }
#endif
  for(;i<len;i++) {
    if(subtract)
      dst[i] -= scalar*src[i];
    else
      dst[i] += scalar*src[i];
  }
}

static inline void negacyclic_mul_coeffs(int64_t out[N], const int64_t c[N],
                                         const int64_t s[N]) {
  size_t j;

  memset(out,0,N*sizeof(*out));
  for(j=0;j<N;j++) {
    if(c[j] == 0)
      continue;
    coefficient_mul_add(&out[j],s,N-j,c[j],0);
    if(j != 0)
      coefficient_mul_add(out,&s[N-j],j,c[j],1);
  }
}

static void challenge_mul_coeffs(int64_t *out, const int64_t c[N],
                                 const int64_t *s, size_t len) {
  size_t i;

#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(len >= LNP22_PARALLEL_MIN_LEN && !omp_in_parallel())
#endif
  for(i=0;i<len;i++)
    negacyclic_mul_coeffs(&out[i*N],c,&s[i*N]);
}

static void centered_coeff_stats(const int64_t *coeffs, size_t count,
                                 int64_t *min_coeff, int64_t *max_coeff,
                                 uint64_t *max_abs_coeff) {
  size_t i;
  int64_t mn = INT64_MAX;
  int64_t mx = INT64_MIN;
  uint64_t maxabs = 0;

  for(i=0;i<count;i++) {
    int64_t c = coeffs[i];
    uint64_t a = (c < 0) ? (uint64_t)(-c) : (uint64_t)c;
    if(c < mn) mn = c;
    if(c > mx) mx = c;
    if(a > maxabs) maxabs = a;
  }

  *min_coeff = mn;
  *max_coeff = mx;
  *max_abs_coeff = maxabs;
}

static long double centered_dot_product_i64(const int64_t *a, const int64_t *b, size_t count) {
  size_t i;
  long double dot = 0.0L;

  for(i=0;i<count;i++)
    dot += (long double)a[i] * (long double)b[i];

  return dot;
}

static double shake_uniform01_seeded(const uint8_t seed[16], uint64_t nonce) {
  shake128incctx st;
  uint8_t in[24];

  memcpy(in,seed,16);
  store_u64_le(&in[16],nonce);
  shake128_inc_init(&st);
  shake128_inc_absorb(&st,in,sizeof(in));
  shake128_inc_finalize(&st);
  return shake_uniform01(&st);
}

static void add_seed_offset(uint8_t out[16], const uint8_t base[16], uint64_t offset) {
  int carry = 0;
  int i;

  memcpy(out,base,16);
  for(i=0;i<8;i++) {
    const uint16_t sum = (uint16_t)out[i] + (uint16_t)(offset & 0xffu) + (uint16_t)carry;
    out[i] = (uint8_t)sum;
    carry = (sum >> 8) & 1;
    offset >>= 8;
  }
  for(i=8;i<16 && carry;i++) {
    const uint16_t sum = (uint16_t)out[i] + (uint16_t)carry;
    out[i] = (uint8_t)sum;
    carry = (sum >> 8) & 1;
  }
}

static void configure_openmp_threads(void) {
#ifdef _OPENMP
  const char *env_threads = getenv("OMP_NUM_THREADS");

  if((env_threads == NULL || env_threads[0] == '\0') && LNP22_DEFAULT_OMP_THREADS > 0) {
    int target = LNP22_DEFAULT_OMP_THREADS;
    int max_threads = omp_get_max_threads();

    if(target > max_threads)
      target = max_threads;
    if(target > 0)
      omp_set_num_threads(target);
  }
  omp_set_dynamic(0);
#endif
}

static int reject_standard_polxvec(const int64_t *zcoeffs, const int64_t *cscoeffs, size_t count,
                                   double sigma, double gamma,
                                   const uint8_t seed[16], uint64_t nonce) {
  const long double logM = rejection_logM(gamma);
  const long double sigma2 = (long double)sigma * (long double)sigma;
  const long double dot_zv = centered_dot_product_i64(zcoeffs,cscoeffs,count);
  const long double dot_vv = centered_dot_product_i64(cscoeffs,cscoeffs,count);
  const long double exponent = (dot_vv - 2.0L * dot_zv) / (2.0L * sigma2);
  double u = shake_uniform01_seeded(seed,nonce);

  if(u < 1e-18)
    u = 1e-18;
  return log((long double)u) > (exponent - logM);
}

static int reject_z2_polxvec(const int64_t *zcoeffs, const int64_t *cscoeffs, size_t count,
                               double sigma, double gamma,
                               const uint8_t seed[16], uint64_t nonce) {
  /* Standard rejection avoids conditioning the transcript on sign(<z,cs>).
   * At gamma=20 its envelope differs from Rej2 by less than one percent. */
  return reject_standard_polxvec(zcoeffs,cscoeffs,count,sigma,gamma,seed,nonce);
}

static void round_ws_init(lnp22_round_ws *ws, const selector_lnp22_ctx *ctx) {
  memset(ws,0,sizeof(*ws));
  ws->y1 = alloc_polxvec(ctx->m1);
  ws->y2 = alloc_polxvec(ctx->m2);
  ws->w = alloc_polxvec(ctx->kmsis);
#if !LNP22_COEFF_RESPONSE
  ws->cs1 = alloc_polxvec(ctx->m1);
  ws->cs2 = alloc_polxvec(ctx->m2);
#endif
  ws->z1coeffs = _malloc(ctx->m1 * N * sizeof(*ws->z1coeffs));
  ws->z2coeffs = _malloc(ctx->m2 * N * sizeof(*ws->z2coeffs));
  ws->cs1coeffs = _malloc(ctx->m1 * N * sizeof(*ws->cs1coeffs));
  ws->cs2coeffs = _malloc(ctx->m2 * N * sizeof(*ws->cs2coeffs));
  ws->c = alloc_polxvec(1);
  ws->vlin = alloc_polxvec(1);
  ws->tq = alloc_polxvec(1);
  ws->vq = alloc_polxvec(1);
  ws->g1 = alloc_polxvec(1);
  ws->g0 = alloc_polxvec(1);
}

static void round_ws_clear(lnp22_round_ws *ws) {
  free(ws->y1);
  free(ws->y2);
  free(ws->w);
#if !LNP22_COEFF_RESPONSE
  free(ws->cs1);
  free(ws->cs2);
#endif
  free(ws->z1coeffs);
  free(ws->z2coeffs);
  free(ws->cs1coeffs);
  free(ws->cs2coeffs);
  free(ws->c);
  free(ws->vlin);
  free(ws->tq);
  free(ws->vq);
  free(ws->g1);
  free(ws->g0);
}

#ifndef SELECTOR_LNP22_INNER_NO_MAIN
int main(int argc, char **argv) {
  selector_lnp22_ctx ctx;
  uint8_t seeded_demo_seed[16];
  polx *s2 = NULL, *y1 = NULL, *y2 = NULL, *tA = NULL, *w = NULL, *z1 = NULL, *z2 = NULL, *a3 = NULL;
#if !LNP22_COEFF_RESPONSE
  polx *cs1 = NULL, *cs2 = NULL;
#endif
  int64_t *z1coeffs = NULL, *z2coeffs = NULL, *cs1coeffs = NULL, *cs2coeffs = NULL;
  int64_t *s1coeffs = NULL, *s2coeffs = NULL;
  int64_t ccoeffs[N];
  polx c[1], c_check[1], vlin[1], tq[1], vq[1];
  polx ufold_public[1], uquad_public[1], a3s2[1];
  polx g1[1], g0[1], lhs[1], rhs[1];
  lnp22_round_ws *batch_ws = NULL;
  struct timespec prove_start, prove_end, verify_start, verify_end, t0, t1;
  double prove_sec = 0.0;
  double verify_sec = 0.0;
  double sample_setup_sec = 0.0;
  double linear_setup_sec = 0.0;
  double open_sec = 0.0;
  double sample_y_sec = 0.0;
  double w_sec = 0.0;
  double vlin_sec = 0.0;
  double quadmask_sec = 0.0;
  double tqvq_sec = 0.0;
  double respond_sec = 0.0;
  double reject_sec = 0.0;
  double batch_eval_sec = 0.0;
  double fs_verify_sec = 0.0;
  double open_verify_sec = 0.0;
  double lin_verify_sec = 0.0;
  double quad_verify_sec = 0.0;
  int64_t z1_min = 0, z1_max = 0, z2_min = 0, z2_max = 0;
  uint64_t z1_maxabs = 0, z2_maxabs = 0;
  size_t rejection_rounds = 0;
  const size_t batch_rounds = LNP22_ENABLE_REJECTION ? LNP22_BATCH_ROUNDS : 1;
  size_t batch_workers = 1;
  uint64_t seed_offset = LNP22_DEFAULT_SEED_OFFSET;
  int ok = 1;

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

  add_seed_offset(seeded_demo_seed,demo_fixed_seed,seed_offset);

  selector_ctx_init(&ctx,seeded_demo_seed);
  prepare_selector_lnp22(&ctx);
#if LNP22_VALIDATE_BASE_WITNESS
  if(first_linear_relation_failure(&ctx,ctx.s1)) {
    fprintf(stderr,"ERROR: base witness does not satisfy the intended linear relations\n");
    ok = 0;
  }
#endif

  s2 = alloc_polxvec(ctx.m2);
  tA = alloc_polxvec(ctx.kmsis);
  a3 = alloc_polxvec(ctx.m2);
  if(batch_rounds <= 1) {
    y1 = alloc_polxvec(ctx.m1);
    y2 = alloc_polxvec(ctx.m2);
    w = alloc_polxvec(ctx.kmsis);
    z1 = alloc_polxvec(ctx.m1);
    z2 = alloc_polxvec(ctx.m2);
#if !LNP22_COEFF_RESPONSE
    cs1 = alloc_polxvec(ctx.m1);
    cs2 = alloc_polxvec(ctx.m2);
#endif
    z1coeffs = _malloc(ctx.m1 * N * sizeof(*z1coeffs));
    z2coeffs = _malloc(ctx.m2 * N * sizeof(*z2coeffs));
    cs1coeffs = _malloc(ctx.m1 * N * sizeof(*cs1coeffs));
    cs2coeffs = _malloc(ctx.m2 * N * sizeof(*cs2coeffs));
  }
#if LNP22_COEFF_RESPONSE
  s1coeffs = _malloc(ctx.m1 * N * sizeof(*s1coeffs));
  s2coeffs = _malloc(ctx.m2 * N * sizeof(*s2coeffs));
  witness_coeffs_to_i64(s1coeffs,&ctx);
#endif
  if(batch_rounds > 1) {
    size_t bi;

    #ifdef _OPENMP
    batch_workers = (size_t)omp_get_max_threads();
    #endif
    if(batch_workers > batch_rounds)
      batch_workers = batch_rounds;
    batch_ws = _malloc(batch_rounds * sizeof(*batch_ws));
    for(bi=0;bi<batch_rounds;bi++)
      round_ws_init(&batch_ws[bi],&ctx);
  }

  clock_gettime(CLOCK_MONOTONIC,&prove_start);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  sample_ternary_polxvec(s2,ctx.m2,ctx.lnp_seed,NONCE_S2_BASE);
#if LNP22_COEFF_RESPONSE
  polxvec_centered_to_i64(s2coeffs,s2,ctx.m2);
#endif
  compute_a3(a3,&ctx);
  inner_product_a3(a3s2,a3,s2,ctx.m2);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  sample_setup_sec = elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(LNP22_DIRECT_CHAIN) {
    compute_tA_rows(tA,&ctx,ctx.s1,s2);
    bind_prechallenge_commitment(&ctx,tA);
  }
  if(prepare_linear_fold_for_proof(&ctx,ufold_public,s2))
    ok = 0;
  compute_public_quad_rhs(uquad_public,&ctx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  linear_setup_sec = elapsed_sec(&t0,&t1);

  if(batch_rounds > 1 && !LNP22_DIRECT_CHAIN)
    trim_ctx_for_batch(&ctx);

  clock_gettime(CLOCK_MONOTONIC,&t0);
#if LNP22_COMBINE_TA_W
  if(!LNP22_DIRECT_CHAIN && (batch_rounds > 1 || batch_ws != NULL))
    compute_tA_rows(tA,&ctx,ctx.s1,s2);
#else
  if(!LNP22_DIRECT_CHAIN)
    compute_tA_rows(tA,&ctx,ctx.s1,s2);
#endif
  if(LNP22_DIRECT_CHAIN) {
    polx direct_check[1];

    derive_direct_chain_fold(&ctx,tA);
#if LNP22_SOUNDNESS_SELFTEST
    if(run_lnp22_soundness_selftest(&ctx,tA))
      ok = 0;
#endif
    compute_quad_g1_coeffs(&ctx);
    quad_value(direct_check,&ctx,ctx.s1);
    if(LNP22_USE_NATIVE_FOLD_BINARY) {
      polx native_term[1];

      polxvec_sprod(native_term,ctx.native_binary_lin2,s2,ctx.m2);
      polx_add(direct_check,direct_check,native_term);
      polx_sub(direct_check,direct_check,uquad_public);
    }
    if(!polx_iszero(direct_check)) {
      fprintf(stderr,"ERROR: direct-chain base quadratic relation failed\n");
      ok = 0;
    }
    if(batch_rounds > 1)
      trim_ctx_for_batch(&ctx);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  open_sec = elapsed_sec(&t0,&t1);

  if(batch_rounds <= 1 || batch_ws == NULL) {
    for(rejection_rounds = 0;;rejection_rounds++) {
      int rej1, rej2;
      size_t i;

      clock_gettime(CLOCK_MONOTONIC,&t0);
      sample_gaussian_polxvec(y1,z1coeffs,ctx.m1,LNP22_LOG2STDEV1,ctx.lnp_seed,NONCE_Y1_BASE + rejection_rounds);
      sample_gaussian_polxvec(y2,z2coeffs,ctx.m2,LNP22_LOG2STDEV2,ctx.lnp_seed,NONCE_Y2_BASE + rejection_rounds);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      sample_y_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
#if LNP22_COMBINE_TA_W
      if(LNP22_DIRECT_CHAIN)
        compute_w_rows(w,&ctx,y1,y2);
      else
        compute_tA_w_rows(tA,w,&ctx,ctx.s1,s2,y1,y2);
#else
      compute_w_rows(w,&ctx,y1,y2);
#endif
      clock_gettime(CLOCK_MONOTONIC,&t1);
      w_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
      folded_linear_eval(vlin,&ctx,y1,y2);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      vlin_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
      quad_mask_terms(g1,g0,&ctx,ctx.s1,y1);
      add_native_binary_s2_mask(g1,&ctx,y2);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      quadmask_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
      polx_add(tq,a3s2,g1);
      inner_product_a3(vq,a3,y2,ctx.m2);
      polx_add(vq,vq,g0);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      tqvq_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
      fs_challenge(c,&ctx,tA,w,vlin,tq,vq);

#if LNP22_COEFF_RESPONSE
      polxvec_centered_to_i64(ccoeffs,c,1);
      challenge_mul_coeffs(cs1coeffs,ccoeffs,s1coeffs,ctx.m1);
      challenge_mul_coeffs(cs2coeffs,ccoeffs,s2coeffs,ctx.m2);
      for(i=0;i<ctx.m1*N;i++)
        z1coeffs[i] += cs1coeffs[i];
      for(i=0;i<ctx.m2*N;i++)
        z2coeffs[i] += cs2coeffs[i];
#else
      polxvec_polx_mul(cs1,c,ctx.s1,ctx.m1);
      polxvec_add(z1,cs1,y1,ctx.m1);
      polxvec_polx_mul(cs2,c,s2,ctx.m2);
      polxvec_add(z2,cs2,y2,ctx.m2);
#endif
      clock_gettime(CLOCK_MONOTONIC,&t1);
      respond_sec += elapsed_sec(&t0,&t1);

      clock_gettime(CLOCK_MONOTONIC,&t0);
      if(LNP22_ENABLE_REJECTION) {
#if !LNP22_COEFF_RESPONSE
        polxvec_centered_to_i64(cs1coeffs,cs1,ctx.m1);
        polxvec_centered_to_i64(cs2coeffs,cs2,ctx.m2);
        for(i=0;i<ctx.m1*N;i++)
          z1coeffs[i] += cs1coeffs[i];
        for(i=0;i<ctx.m2*N;i++)
          z2coeffs[i] += cs2coeffs[i];
#endif
        rej1 = reject_standard_polxvec(z1coeffs,cs1coeffs,ctx.m1*N,
                                       sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                       LNP22_GAMMA1,
                                       ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + rejection_rounds);
        rej2 = reject_z2_polxvec(z2coeffs,cs2coeffs,ctx.m2*N,
                                   sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                   LNP22_GAMMA2,
                                   ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + rejection_rounds);
      } else {
        rej1 = 0;
        rej2 = 0;
      }
      clock_gettime(CLOCK_MONOTONIC,&t1);
      reject_sec += elapsed_sec(&t0,&t1);
      if(!rej1 && !rej2) {
#if LNP22_COEFF_RESPONSE
        clock_gettime(CLOCK_MONOTONIC,&t0);
        polxvec_fromint64vec(z1,ctx.m1,1,z1coeffs);
        polxvec_fromint64vec(z2,ctx.m2,1,z2coeffs);
        clock_gettime(CLOCK_MONOTONIC,&t1);
        respond_sec += elapsed_sec(&t0,&t1);
#endif
        break;
      }
    }
  } else {
    size_t batch_base = 0;
    int found = 0;
    uint8_t *accept = _malloc(batch_rounds);
    size_t found_bi = 0;

    for(;;) {
      size_t bi;

      clock_gettime(CLOCK_MONOTONIC,&t0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(batch_workers)
#endif
      for(bi=0;bi<batch_rounds;bi++) {
        lnp22_round_ws *ws = &batch_ws[bi];
        sample_gaussian_polxvec(ws->y1,ws->z1coeffs,ctx.m1,LNP22_LOG2STDEV1,
                                ctx.lnp_seed,NONCE_Y1_BASE + batch_base + bi);
        sample_gaussian_polxvec(ws->y2,ws->z2coeffs,ctx.m2,LNP22_LOG2STDEV2,
                                ctx.lnp_seed,NONCE_Y2_BASE + batch_base + bi);
        compute_w_rows(ws->w,&ctx,ws->y1,ws->y2);
        folded_linear_eval(ws->vlin,&ctx,ws->y1,ws->y2);
        quad_mask_terms(ws->g1,ws->g0,&ctx,ctx.s1,ws->y1);
        add_native_binary_s2_mask(ws->g1,&ctx,ws->y2);
        polx_add(ws->tq,a3s2,ws->g1);
        inner_product_a3(ws->vq,a3,ws->y2,ctx.m2);
        polx_add(ws->vq,ws->vq,ws->g0);
        fs_challenge(ws->c,&ctx,tA,ws->w,ws->vlin,ws->tq,ws->vq);
#if LNP22_COEFF_RESPONSE
        polxvec_centered_to_i64(ws->ccoeffs,ws->c,1);
        challenge_mul_coeffs(ws->cs1coeffs,ws->ccoeffs,s1coeffs,ctx.m1);
        challenge_mul_coeffs(ws->cs2coeffs,ws->ccoeffs,s2coeffs,ctx.m2);
#else
        polxvec_polx_mul(ws->cs1,ws->c,ctx.s1,ctx.m1);
        polxvec_polx_mul(ws->cs2,ws->c,s2,ctx.m2);
#endif
        size_t i;

#if !LNP22_COEFF_RESPONSE
        polxvec_centered_to_i64(ws->cs1coeffs,ws->cs1,ctx.m1);
        polxvec_centered_to_i64(ws->cs2coeffs,ws->cs2,ctx.m2);
#endif
        for(i=0;i<ctx.m1*N;i++)
          ws->z1coeffs[i] += ws->cs1coeffs[i];
        for(i=0;i<ctx.m2*N;i++)
          ws->z2coeffs[i] += ws->cs2coeffs[i];
        ws->rej1 = reject_standard_polxvec(ws->z1coeffs,ws->cs1coeffs,ctx.m1*N,
                                           sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                           LNP22_GAMMA1,
                                           ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + batch_base + bi);
        ws->rej2 = reject_z2_polxvec(ws->z2coeffs,ws->cs2coeffs,ctx.m2*N,
                                       sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                       LNP22_GAMMA2,
                                       ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + batch_base + bi);
        accept[bi] = (uint8_t)(!ws->rej1 && !ws->rej2);
      }
      clock_gettime(CLOCK_MONOTONIC,&t1);
      batch_eval_sec += elapsed_sec(&t0,&t1);

      for(bi=0;bi<batch_rounds;bi++) {
        if(accept[bi]) {
          rejection_rounds = batch_base + bi;
          found_bi = bi;
          found = 1;
          break;
        }
      }
      if(found)
        break;
      batch_base += batch_rounds;
    }

    y1 = alloc_polxvec(ctx.m1);
    y2 = alloc_polxvec(ctx.m2);
    w = alloc_polxvec(ctx.kmsis);
    z1 = alloc_polxvec(ctx.m1);
    z2 = alloc_polxvec(ctx.m2);
#if !LNP22_COEFF_RESPONSE
    cs1 = alloc_polxvec(ctx.m1);
    cs2 = alloc_polxvec(ctx.m2);
#endif
    z1coeffs = _malloc(ctx.m1 * N * sizeof(*z1coeffs));
    z2coeffs = _malloc(ctx.m2 * N * sizeof(*z2coeffs));
    cs1coeffs = _malloc(ctx.m1 * N * sizeof(*cs1coeffs));
    cs2coeffs = _malloc(ctx.m2 * N * sizeof(*cs2coeffs));
    {
      const lnp22_round_ws *winner = &batch_ws[found_bi];
      size_t i;

      polxvec_copy(y1,winner->y1,ctx.m1);
      polxvec_copy(y2,winner->y2,ctx.m2);
      polxvec_copy(w,winner->w,ctx.kmsis);
#if !LNP22_COEFF_RESPONSE
      polxvec_copy(cs1,winner->cs1,ctx.m1);
      polxvec_copy(cs2,winner->cs2,ctx.m2);
#endif
      polxvec_copy(c,winner->c,1);
      polxvec_copy(vlin,winner->vlin,1);
      polxvec_copy(tq,winner->tq,1);
      polxvec_copy(vq,winner->vq,1);
      polxvec_copy(g1,winner->g1,1);
      polxvec_copy(g0,winner->g0,1);
      for(i=0;i<ctx.m1*N;i++) {
        z1coeffs[i] = winner->z1coeffs[i];
        cs1coeffs[i] = winner->cs1coeffs[i];
      }
      for(i=0;i<ctx.m2*N;i++) {
        z2coeffs[i] = winner->z2coeffs[i];
        cs2coeffs[i] = winner->cs2coeffs[i];
      }
#if LNP22_COEFF_RESPONSE
      polxvec_fromint64vec(z1,ctx.m1,1,z1coeffs);
      polxvec_fromint64vec(z2,ctx.m2,1,z2coeffs);
#else
      polxvec_add(z1,cs1,y1,ctx.m1);
      polxvec_add(z2,cs2,y2,ctx.m2);
#endif
    }

    free(accept);
  }
  clock_gettime(CLOCK_MONOTONIC,&prove_end);
  prove_sec = elapsed_sec(&prove_start,&prove_end);

  if(!LNP22_ENABLE_REJECTION) {
    polxvec_centered_to_i64(z1coeffs,z1,ctx.m1);
    polxvec_centered_to_i64(z2coeffs,z2,ctx.m2);
  }
  centered_coeff_stats(z1coeffs,ctx.m1*N,&z1_min,&z1_max,&z1_maxabs);
  centered_coeff_stats(z2coeffs,ctx.m2*N,&z2_min,&z2_max,&z2_maxabs);

  clock_gettime(CLOCK_MONOTONIC,&verify_start);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(validate_rbit_const_transcript(&ctx))
    ok = 0;
  fs_challenge(c_check,&ctx,tA,w,vlin,tq,vq);
  if(memcmp(c,c_check,sizeof(polx)) != 0) {
    fprintf(stderr,"ERROR: Fiat-Shamir challenge mismatch\n");
    ok = 0;
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  fs_verify_sec = elapsed_sec(&t0,&t1);

  if(ok) {
    int opening_ok = 1;
    clock_gettime(CLOCK_MONOTONIC,&t0);
#if defined(_OPENMP) && LNP22_ROW_PARALLEL
#pragma omp parallel for schedule(static) reduction(&&:opening_ok)
#endif
    for(long long rr=0;rr<(long long)ctx.kmsis;rr++) {
      const size_t i = (size_t)rr;
      polx row_eval[1], tmp_row[1], rhs_row[1];
      if(ctx.cache_a1open)
        polxvec_sprod(row_eval,&ctx.A1open[i*ctx.m1],z1,ctx.m1);
      else
        stream_rowdot_almostuniform(row_eval,ctx.lnp_seed,NONCE_LNP_A1_BASE,i,ctx.m1,z1);
      polxvec_sprod(tmp_row,&ctx.A2open[i*ctx.m2prime],z2,ctx.m2prime);
      polx_add(row_eval,row_eval,tmp_row);
      polx_add(row_eval,row_eval,&z2[ctx.m2prime + i]);
      polx_mul(rhs_row,c,&tA[i]);
      polx_sub(row_eval,row_eval,rhs_row);
      polx_sub(row_eval,row_eval,&w[i]);
      if(!polx_iszero(row_eval))
        opening_ok = 0;
    }
    if(!opening_ok) {
      fprintf(stderr,"ERROR: opening check failed\n");
      ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    open_verify_sec = elapsed_sec(&t0,&t1);
  }

  if(ok) {
    clock_gettime(CLOCK_MONOTONIC,&t0);
    folded_linear_eval(lhs,&ctx,z1,z2);
    polx_mul(rhs,c,ufold_public);
    polx_sub(lhs,lhs,rhs);
    polx_sub(lhs,lhs,vlin);
    if(!polx_iszero(lhs)) {
      fprintf(stderr,"ERROR: folded linear check failed\n");
      ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    lin_verify_sec = elapsed_sec(&t0,&t1);
  }

  if(ok) {
    polx lin_term[1], bil_term[1], bits_term[1], mask_term[1], a3z2[1], c2[1], c2u[1];
    clock_gettime(CLOCK_MONOTONIC,&t0);
    if(LNP22_DIRECT_CHAIN) {
      size_t t;

      polxvec_sprod(lin_term,ctx.direct_lin_B,&z1[ctx.off[LNP22_B_BLOCK]],ctx.b_len);
      polx_refresh(lin_term);
      polx_mul(lhs,c,lin_term);
      for(t=0;t<ctx.Lchain-1;t++) {
        polxvec_sprod(bil_term,&ctx.directD[t*ctx.gdigits],
                      &z1[ctx.off[LNP22_B_BLOCK] + (t+1)*ctx.gdigits],ctx.gdigits);
        polx_refresh(bil_term);
        polx_mul(bil_term,bil_term,
                 &z1[ctx.off[LNP22_X_BLOCK] + selector_x_pos(&ctx,t)]);
        polx_sub(lhs,lhs,bil_term);
      }
    } else {
      linZ_value(lin_term,&ctx,&z1[ctx.off[LNP22_Z_BLOCK]]);
      polx_mul(lin_term,c,lin_term);
      bilinear_value(bil_term,&ctx,&z1[ctx.off[LNP22_X_BLOCK]],
                     &z1[ctx.off[LNP22_Y0_BLOCK]],&z1[ctx.off[LNP22_Y1_BLOCK]]);
      polx_sub(lhs,lin_term,bil_term);
    }
    if(LNP22_PROVE_X_BITS) {
      x_bits_verifier_value(bits_term,&ctx,&z1[ctx.off[LNP22_X_BLOCK]],c);
      polx_add(lhs,lhs,bits_term);
    }
    if(LNP22_USE_RBIT_BINARY) {
      rhat_bits_verifier_value(bits_term,&ctx,&z1[ctx.off[LNP22_RHAT_BITS_BLOCK]],c);
      polx_add(lhs,lhs,bits_term);
    }
    if(LNP22_USE_NATIVE_FOLD_BINARY) {
      rhat_native_binary_verifier_value(
          bits_term,&ctx,&z1[ctx.off[LNP22_RHAT_BLOCK]],c);
      polx_add(lhs,lhs,bits_term);
      polxvec_sprod(bits_term,ctx.native_binary_lin1,
                    &z1[ctx.off[LNP22_X_BLOCK]],ctx.x_len);
      polx_mul(bits_term,c,bits_term);
      polx_add(lhs,lhs,bits_term);
      polxvec_sprod(bits_term,ctx.native_binary_lin2,z2,ctx.m2);
      polx_mul(bits_term,c,bits_term);
      polx_add(lhs,lhs,bits_term);
    }
    polx_mul(c2,c,c);
    polx_mul(c2u,c2,uquad_public);
    polx_sub(lhs,lhs,c2u);
    inner_product_a3(a3z2,a3,z2,ctx.m2);
    polx_mul(mask_term,c,tq);
    polx_sub(mask_term,mask_term,a3z2);
    polx_sub(lhs,lhs,mask_term);
    polx_sub(lhs,lhs,vq);
    if(!polx_iszero(lhs)) {
      fprintf(stderr,"ERROR: quadratic check failed\n");
      ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    quad_verify_sec = elapsed_sec(&t0,&t1);
  }
  clock_gettime(CLOCK_MONOTONIC,&verify_end);
  verify_sec = elapsed_sec(&verify_start,&verify_end);

  printf("Testing selector LNP22 inner proof (uncompressed, explicit w/v)\n\n");
  printf("Parameters:\n");
  printf("  LOGQ              = %d\n",LOGQ);
  printf("  L                 = %d\n",DEMO_L);
  printf("  GADGET_BASE       = %d\n",GADGET_BASE);
  printf("  Y split gamma     = %d\n",Y_SPLIT_GAMMA);
  printf("  kmsis             = %d\n",LNP22_KMSIS);
  printf("  m2                = %d\n",LNP22_M2);
  printf("  omega             = %d\n",LNP22_OMEGA);
  printf("  RHAT binary       = %s\n",
         LNP22_USE_NATIVE_FOLD_BINARY ? "native automorphism fold" :
                                        "RBIT coefficient lift");
  printf("  X bit proof       = %s\n",LNP22_PROVE_X_BITS ? "yes" : "no");
  printf("  combine tA/w      = %s\n",LNP22_COMBINE_TA_W ? "yes" : "no");
  printf("  const-check folds = %zu\n",ctx.rbit_const_lambda);
  printf("  rejection enabled = %s\n",LNP22_ENABLE_REJECTION ? "yes" : "no");
  printf("  stdev1            = %.1f\n",sigma_from_log2stdev(LNP22_LOG2STDEV1));
  printf("  stdev2            = %.1f\n",sigma_from_log2stdev(LNP22_LOG2STDEV2));
  printf("  witness ranks     = B=%zu X=%zu Z=%zu Y0=%zu Y1=%zu RHAT=%zu RBIT=%zu\n",
         ctx.b_len,ctx.x_len,ctx.z_len,ctx.y_len,ctx.y_len,ctx.rhat_len,ctx.rhat_bits_len);
  printf("  UMC witness ranks = EUA=%zu EUB=%zu R=%zu VA=%zu VB=%zu\n",
         ctx.umc_eua_len,ctx.umc_eub_len,ctx.umc_r_len,ctx.umc_va_len,ctx.umc_vb_len);
  printf("  s1 length         = %zu\n",ctx.m1);
  printf("  C_X rows          = %zu\n",ctx.cx_rows);
  printf("  preprocessing R_j = %zu\n",ctx.preproc_count);
  printf("  selected J        = %zu\n",ctx.selected_j + 1);
  printf("  UMC v_A range     = [%lld, %lld], maxabs=%llu\n",
         (long long)ctx.umc_va_min,(long long)ctx.umc_va_max,
         (unsigned long long)ctx.umc_va_maxabs);
  printf("  UMC v_B range     = [%lld, %lld], maxabs=%llu\n",
         (long long)ctx.umc_vb_min,(long long)ctx.umc_vb_max,
         (unsigned long long)ctx.umc_vb_maxabs);
  printf("  folded linear eqs = %zu\n\n",ctx.nlin);
  printf("  seed offset       = %llu\n\n",(unsigned long long)seed_offset);
  printf("Timing:\n");
  printf("  proving time      = %.6f s\n",prove_sec);
  printf("  verification time = %.6f s\n\n",verify_sec);
  printf("Proving breakdown:\n");
  printf("  sample s2/a3      = %.6f s\n",sample_setup_sec);
  printf("  linear setup      = %.6f s\n",linear_setup_sec);
  printf("  compute tA        = %.6f s\n",open_sec);
  printf("  sample y1/y2      = %.6f s\n",sample_y_sec);
  printf("  compute w         = %.6f s\n",w_sec);
  printf("  compute v_lin     = %.6f s\n",vlin_sec);
  printf("  compute g1/g0     = %.6f s\n",quadmask_sec);
  printf("  compute tq/vq     = %.6f s\n",tqvq_sec);
  printf("  FS + z1/z2        = %.6f s\n",respond_sec);
  printf("  rejection checks  = %.6f s\n",reject_sec);
  if(batch_rounds > 1)
    printf("  batch eval total  = %.6f s\n",batch_eval_sec);
  printf("\n");
  printf("  rejection rounds  = %zu\n\n",rejection_rounds + 1);
  printf("Verification breakdown:\n");
  printf("  recompute FS      = %.6f s\n",fs_verify_sec);
  printf("  opening check     = %.6f s\n",open_verify_sec);
  printf("  linear check      = %.6f s\n",lin_verify_sec);
  printf("  quadratic check   = %.6f s\n\n",quad_verify_sec);
  printf("Centered coefficient ranges:\n");
  printf("  z1 coeff range    = [%lld, %lld], max|coeff| = %llu\n",
         (long long)z1_min,(long long)z1_max,(unsigned long long)z1_maxabs);
  printf("  z2 coeff range    = [%lld, %lld], max|coeff| = %llu\n\n",
         (long long)z2_min,(long long)z2_max,(unsigned long long)z2_maxabs);
  print_size_breakdown(&ctx);
  printf("\nVerifier result: %s\n",ok ? "PASS" : "FAIL");
  printf("Note: this driver keeps w/v explicit and skips D/gamma compression.\n");
  printf("Note: selector challenges are derived from the public statement (C_X, C), and folded/public RHS values are recomputed from public data.\n");
  if(!LNP22_ENABLE_REJECTION)
    printf("Note: rejection sampling is disabled in this build; this is suitable for profiling/debugging, not final ZK measurements.\n");

  free(s2);
  free(y1);
  free(y2);
  free(tA);
  free(w);
  free(z1);
  free(z2);
  free(a3);
#if !LNP22_COEFF_RESPONSE
  free(cs1);
  free(cs2);
#endif
  free(z1coeffs);
  free(z2coeffs);
  free(cs1coeffs);
  free(cs2coeffs);
  free(s1coeffs);
  free(s2coeffs);
  if(batch_ws != NULL) {
    size_t bi;
    for(bi=0;bi<batch_rounds;bi++)
      round_ws_clear(&batch_ws[bi]);
    free(batch_ws);
  }
  selector_ctx_clear(&ctx);
  free_comkey();
  return ok ? 0 : 1;
}
#endif
