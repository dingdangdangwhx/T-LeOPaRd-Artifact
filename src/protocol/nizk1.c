#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "nizk1.h"

#ifndef LNP22_BATCH_ROUNDS
#define LNP22_BATCH_ROUNDS 1
#endif

#ifndef LNP22_CACHE_A1OPEN
#define LNP22_CACHE_A1OPEN 0
#endif

#ifndef LEOPARD_REQUEST_SELF_CHECK
#define LEOPARD_REQUEST_SELF_CHECK 0
#endif

#ifndef LEOPARD_PUBLIC_ROW_PARALLEL
#define LEOPARD_PUBLIC_ROW_PARALLEL 0
#endif

#ifndef LEOPARD_PUBLIC_PARALLEL_MIN_ROWS
#define LEOPARD_PUBLIC_PARALLEL_MIN_ROWS 8
#endif

#ifndef LEOPARD_ROUND_P
#define LEOPARD_ROUND_P 8
#endif

int labrador_latest_outer_main(int argc, char **argv);
#define main labrador_latest_outer_main
#include "client_outer.c"
#undef main

struct leopard_nizk1_request {
  size_t selected_j;
  polx *cr;
  polx *cx;
  polx *tA;
  polx *w;
  polx *rbit_const_tg;
  polx *rbit_const_h;
  polx *vlin;
  polx *tq;
  polx *vq;
  size_t cr_len;
  size_t cx_len;
  size_t kmsis;
  size_t rbit_const_len;
  leopard_coeff *cr_coeffs;
  leopard_coeff *cx_coeffs;
  leopard_coeff *tA_coeffs;
  leopard_coeff *w_coeffs;
  leopard_coeff *rbit_const_tg_coeffs;
  leopard_coeff *rbit_const_h_coeffs;
  leopard_coeff *vlin_coeffs;
  leopard_coeff *tq_coeffs;
  leopard_coeff *vq_coeffs;
  uint8_t ar_seed[16];
  uint8_t *proof;
  size_t proof_len;
  uint64_t outer_betasq;
};

struct leopard_nizk1_state {
  size_t selected_j;
  polx *rhat;
  size_t rhat_len;
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];
  leopard_coeff *rhat_coeffs;
};

static leopard_nizk1_request_profile g_last_request_profile;
static leopard_nizk1_verify_profile g_last_verify_profile;

typedef struct {
  double fold_sec;
  double quad_sec;
  double cnst1_sec;
  double cnst2_sec;
  double cnst3_sec;
  double hash_sec;
} leopard_outer_statement_profile;

static double leopard_elapsed_sec(const struct timespec *start,
                                  const struct timespec *end) {
  return (double)(end->tv_sec - start->tv_sec) +
         1e-9 * (double)(end->tv_nsec - start->tv_nsec);
}

#include "public_fastpath.h"

static void *checked_malloc(size_t size, const char *what) {
  void *ptr = malloc(size);

  if(ptr == NULL) {
    fprintf(stderr,"allocation failed for %s (%zu bytes)\n",what,size);
    exit(1);
  }
  return ptr;
}

static void leopard_tag_hash(uint8_t out[LEOPARD_TAG_BYTES],
                             const uint8_t *tag, size_t tag_len) {
  uint64_t s0 = 0x243f6a8885a308d3ULL;
  uint64_t s1 = 0x13198a2e03707344ULL;
  size_t i;

  for(i=0;i<tag_len;i++) {
    s0 ^= (uint64_t)tag[i];
    s0 *= 0x100000001b3ULL;
    s1 ^= (uint64_t)tag[i] + 0x9e3779b97f4a7c15ULL + (s1 << 6) + (s1 >> 2);
    s1 *= 0xff51afd7ed558ccdULL;
  }
  for(i=0;i<8;i++) {
    out[i] = (uint8_t)(s0 >> (8*i));
    out[8 + i] = (uint8_t)(s1 >> (8*i));
  }
}

static void leopard_expand_input_bits(uint8_t out[LEOPARD_X_BITS],
                                      const uint8_t in[LEOPARD_X_BYTES]) {
  size_t i;

  for(i=0;i<LEOPARD_X_BITS;i++)
    out[i] = (uint8_t)((in[i >> 3] >> (i & 7)) & 1u);
}

static void export_polx_coeffs(leopard_coeff *out, const polx *src, size_t len) {
#if LOGQ > 64
  polxvec_centered_to_i128(out,src,len);
#else
  polxvec_centered_to_i64(out,src,len);
#endif
}

static leopard_coeff leopard_center_q_i64(leopard_coeff x) {
  const __int128 q = (__int128)LEOPARD_MOD_Q;
  __int128 y = (__int128)x % q;

  if(y > q / 2)
    y -= q;
  if(y < -(q / 2))
    y += q;
  return (leopard_coeff)y;
}

static void import_polx_coeffs(polx *dst, const leopard_coeff *src, size_t len) {
#if LOGQ > 64
  for(size_t i=0;i<len;i++) {
    polz z;
    for(size_t j=0;j<N;j++)
      polz_setcoeff_fromint128(&z,leopard_center_q_i64(src[i*N+j]),j);
    polz_topolx(&dst[i],&z);
  }
#else
  const size_t coeff_count = len * (size_t)N;
  leopard_coeff *centered;
  size_t i;

  centered = checked_malloc(coeff_count * sizeof(*centered),"centered polx import");
  for(i=0;i<coeff_count;i++)
    centered[i] = leopard_center_q_i64(src[i]);
  polxvec_fromint64vec(dst,len,1,centered);
  free(centered);
#endif
}

static void export_poly_coeffs(leopard_coeff *out, const poly *src, size_t len) {
  size_t i, j;

  for(i=0;i<len;i++) {
    for(j=0;j<N;j++)
      out[i * N + j] = src[i].vec->c[j];
  }
}

static void import_poly_coeffs(poly *dst, const leopard_coeff *src, size_t len) {
  size_t i, j;

  for(i=0;i<len;i++) {
    for(j=0;j<N;j++)
      dst[i].vec->c[j] = (int16_t)src[i * N + j];
  }
}


static leopard_coeff leopard_center_q_i128(__int128 x) {
  const __int128 q = (__int128)LEOPARD_MOD_Q;
  __int128 y = x % q;

  if(y > q / 2)
    y -= q;
  if(y < -(q / 2))
    y += q;
  return (leopard_coeff)y;
}

static leopard_coeff leopard_round_div_nearest_i128(__int128 num, __int128 den) {
  if(num >= 0)
    return (leopard_coeff)((num + den / 2) / den);
  return -(leopard_coeff)((-num + den / 2) / den);
}

static uint8_t leopard_round_coeff_p(leopard_coeff coeff) {
  leopard_coeff rounded = leopard_round_div_nearest_i128(
      (__int128)LEOPARD_ROUND_P * (__int128)coeff,
      (__int128)LEOPARD_MOD_Q);
  leopard_coeff reduced = rounded % LEOPARD_ROUND_P;

  if(reduced < 0)
    reduced += LEOPARD_ROUND_P;
  return (uint8_t)reduced;
}

static void leopard_round_polx_p(uint8_t out[LEOPARD_RING_DEG],
                                 const polx *z) {
  leopard_coeff coeffs[LEOPARD_RING_DEG];
  size_t i;

  export_polx_coeffs(coeffs,z,1);
  for(i=0;i<LEOPARD_RING_DEG;i++)
    out[i] = leopard_round_coeff_p(coeffs[i]);
}

static int leopard_fast_public_compute_preprocessed(
    leopard_nizk1_fast_public_profile *out_prof,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS],
    const leopard_request_data *expected_req, polx *out_cx) {
  const int public_threads = leopard_public_threads();
  const size_t m = DEMO_M;
  const size_t gdigits = m * gadget_len_base3();
  const size_t block = gdigits / m;
  const size_t rhat_len = UMC_ELL + m;
  polx *A0 = alloc_polxvec(m * gdigits);
  polx *A1sel = alloc_polxvec(m * gdigits);
  polx *Arhash = alloc_polxvec(rhat_len * m);
  polx *rhat = alloc_polxvec(rhat_len);
  polx *bnextx = alloc_polxvec(gdigits);
  polx *cx = alloc_polxvec(m);
  poly *g0poly = alloc_polyvec(gdigits);
  poly *g1poly = alloc_polyvec(gdigits);
  poly *bpoly0 = alloc_polyvec(gdigits);
  poly *bpoly1 = alloc_polyvec(gdigits);
  poly *rhatpoly = alloc_polyvec(rhat_len);
  poly *bnext;
  uint8_t statement_seed[LEOPARD_TAG_BYTES];
  uint64_t nonce = 0;
  size_t row, j;
  int rc = 0;
  struct timespec t0, t1;

  memset(out_prof,0,sizeof(*out_prof));

  clock_gettime(CLOCK_MONOTONIC,&t0);
  memcpy(statement_seed,tag_hash,sizeof(statement_seed));
  leopard_public_almostuniform(A0,m * gdigits,statement_seed,nonce,public_threads);
  nonce += (m * gdigits + 31) / 32;
  leopard_public_almostuniform(A1sel,m * gdigits,statement_seed,nonce,public_threads);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<m;row++) {
    leopard_public_decompose(&g0poly[row * block],&A0[row * gdigits],block);
    leopard_public_decompose(&g1poly[row * block],&A1sel[row * gdigits],block);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  out_prof->selector_material_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  import_poly_coeffs(rhatpoly,
                     &pre->rhat_coeffs[slot_j * pre->rhat_rows * LEOPARD_RING_DEG],
                     pre->rhat_rows);
  polxvec_frompolyvec(rhat,rhatpoly,rhat_len);
  leopard_public_almostuniform(Arhash,rhat_len * m,&pre->ar_seeds[16 * slot_j],
                               NONCE_AR_BASE,public_threads);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  out_prof->slot_load_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  bnext = leopard_public_chain(bnextx,bpoly0,bpoly1,g0poly,g1poly,A0,A1sel,
                                xbits,public_threads,1,NULL);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  out_prof->chain_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  polxvec_frompolyvec(bnextx,bnext,gdigits);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<m;row++) {
    polx prod[1], gb[1];

    polxvec_sprod(prod,&Arhash[row * rhat_len],rhat,rhat_len);
    polxvec_setzero(gb,1);
    for(j=0;j<block;j++)
      polx_scale_add(gb,&bnextx[row * block + j],pow3_u64(j));
    polx_add(&cx[row],prod,gb);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  out_prof->cx_sec = leopard_elapsed_sec(&t0,&t1);
  out_prof->total_sec = out_prof->selector_material_sec +
                        out_prof->slot_load_sec +
                        out_prof->chain_sec +
                        out_prof->cx_sec;

  if(out_cx != NULL)
    polxvec_copy(out_cx,cx,m);
  if(expected_req != NULL) {
    leopard_coeff *cx_coeffs = checked_malloc(m * LEOPARD_RING_DEG * sizeof(*cx_coeffs),
                                        "fast public cx coeffs");
    export_polx_coeffs(cx_coeffs,cx,m);
    if(expected_req->cx_rows * expected_req->cx_cols != m ||
       memcmp(cx_coeffs,expected_req->cx_coeffs,
              m * LEOPARD_RING_DEG * sizeof(*cx_coeffs)) != 0) {
      rc = 1;
    }
    free(cx_coeffs);
  }

  free(A0);
  free(A1sel);
  free(Arhash);
  free(rhat);
  free(bnextx);
  free(cx);
  free(g0poly);
  free(g1poly);
  free(bpoly0);
  free(bpoly1);
  free(rhatpoly);
  return rc;
}


static int leopard_compute_bxt_fast(polx *bxt,
                                    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                    const uint8_t xbits[LEOPARD_X_BITS]) {
  const int public_threads = leopard_public_threads();
  const size_t m = DEMO_M;
  const size_t gdigits = m * gadget_len_base3();
  const size_t block = gdigits / m;
  polx *A0 = alloc_polxvec(m * gdigits);
  polx *A1sel = alloc_polxvec(m * gdigits);
  polx *bnextx = alloc_polxvec(gdigits);
  poly *g0poly = alloc_polyvec(gdigits);
  poly *g1poly = alloc_polyvec(gdigits);
  poly *bpoly0 = alloc_polyvec(gdigits);
  poly *bpoly1 = alloc_polyvec(gdigits);
  poly *bnext;
  uint8_t statement_seed[LEOPARD_TAG_BYTES];
  uint64_t nonce = 0;
  size_t row, j;

  memcpy(statement_seed,tag_hash,sizeof(statement_seed));
  leopard_public_almostuniform(A0,m * gdigits,statement_seed,nonce,public_threads);
  nonce += (m * gdigits + 31) / 32;
  leopard_public_almostuniform(A1sel,m * gdigits,statement_seed,nonce,public_threads);

#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<m;row++) {
    leopard_public_decompose(&g0poly[row * block],&A0[row * gdigits],block);
    leopard_public_decompose(&g1poly[row * block],&A1sel[row * gdigits],block);
  }

  bnext = leopard_public_chain(bnextx,bpoly0,bpoly1,g0poly,g1poly,A0,A1sel,
                                xbits,public_threads,1,NULL);

  polxvec_frompolyvec(bnextx,bnext,gdigits);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<m;row++) {
    polxvec_setzero(&bxt[row],1);
    for(j=0;j<block;j++)
      polx_scale_add(&bxt[row],&bnextx[row * block + j],pow3_u64(j));
  }

  free(A0);
  free(A1sel);
  free(bnextx);
  free(g0poly);
  free(g1poly);
  free(bpoly0);
  free(bpoly1);
  return 0;
}

static void leopard_hash_output(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                const uint8_t xbits[LEOPARD_X_BITS],
                                const polx *z) {
  shake128incctx st;
  uint8_t z_p[LEOPARD_RING_DEG];

  leopard_round_polx_p(z_p,z);
  shake128_inc_init(&st);
  absorb_bytes(&st,"leopard-roz-p8",14);
  absorb_bytes(&st,tag_hash,LEOPARD_TAG_BYTES);
  absorb_bytes(&st,xbits,LEOPARD_X_BITS);
  absorb_bytes(&st,z_p,sizeof(z_p));
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(y,LEOPARD_OUTPUT_BYTES,&st);
}


static void leopard_hash_output_coeffs(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                       const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                       const uint8_t xbits[LEOPARD_X_BITS],
                                       const leopard_coeff z_coeffs[LEOPARD_RING_DEG]) {
  shake128incctx st;
  uint8_t z_p[LEOPARD_RING_DEG];
  size_t i;

  for(i=0;i<LEOPARD_RING_DEG;i++)
    z_p[i] = leopard_round_coeff_p(z_coeffs[i]);
  shake128_inc_init(&st);
  absorb_bytes(&st,"leopard-roz-p8",14);
  absorb_bytes(&st,tag_hash,LEOPARD_TAG_BYTES);
  absorb_bytes(&st,xbits,LEOPARD_X_BITS);
  absorb_bytes(&st,z_p,sizeof(z_p));
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(y,LEOPARD_OUTPUT_BYTES,&st);
}

static int leopard_finalize_output_binary_rhat(
    uint8_t y[LEOPARD_OUTPUT_BYTES],
    const leopard_state_data *st,
    const leopard_coeff *ux_coeffs, size_t ux_rows,
    const leopard_coeff *vk_coeffs, size_t vk_rows) {
  __int128 acc[LEOPARD_RING_DEG] = {0};
  leopard_coeff z_coeffs[LEOPARD_RING_DEG];
  size_t row, i, j;

  if(ux_rows != 1 || vk_rows != st->rhat_rows)
    return 1;

  for(row=0;row<st->rhat_rows;row++) {
    const leopard_coeff *rh = &st->rhat_coeffs[row * LEOPARD_RING_DEG];
    const leopard_coeff *vk = &vk_coeffs[row * LEOPARD_RING_DEG];

    for(i=0;i<LEOPARD_RING_DEG;i++) {
      const leopard_coeff bit = rh[i];

      if(bit == 0)
        continue;
      if(bit != 1)
        return 1;
      for(j=0;j<LEOPARD_RING_DEG;j++) {
        const size_t idx = i + j;
        const leopard_coeff v = leopard_center_q_i64(vk[j]);

        if(idx < LEOPARD_RING_DEG)
          acc[idx] += (__int128)v;
        else
          acc[idx - LEOPARD_RING_DEG] -= (__int128)v;
      }
    }
  }

  for(i=0;i<LEOPARD_RING_DEG;i++) {
    const leopard_coeff ux = leopard_center_q_i64(ux_coeffs[i]);
    z_coeffs[i] = leopard_center_q_i128((__int128)ux - acc[i]);
  }
  leopard_hash_output_coeffs(y,st->tag_hash,st->xbits,z_coeffs);
  return 0;
}

static leopard_nizk1_request *leopard_request_alloc(size_t cr_len, size_t cx_len,
                                                    size_t kmsis,
                                                    size_t rbit_const_len) {
  leopard_nizk1_request *req;

  req = checked_malloc(sizeof(*req),"request");
  memset(req,0,sizeof(*req));
  req->cr_len = cr_len;
  req->cx_len = cx_len;
  req->kmsis = kmsis;
  req->rbit_const_len = rbit_const_len;
  req->cr = alloc_polxvec(cr_len);
  req->cx = alloc_polxvec(cx_len);
  req->tA = alloc_polxvec(kmsis);
  req->w = alloc_polxvec(kmsis);
  req->rbit_const_tg = alloc_polxvec(rbit_const_len);
  req->rbit_const_h = alloc_polxvec(rbit_const_len);
  req->vlin = alloc_polxvec(1);
  req->tq = alloc_polxvec(1);
  req->vq = alloc_polxvec(1);
  req->cr_coeffs = checked_malloc(cr_len * LEOPARD_RING_DEG * sizeof(*req->cr_coeffs),
                                  "request cr coeffs");
  req->cx_coeffs = checked_malloc(cx_len * LEOPARD_RING_DEG * sizeof(*req->cx_coeffs),
                                  "request cx coeffs");
  req->tA_coeffs = checked_malloc(kmsis * LEOPARD_RING_DEG * sizeof(*req->tA_coeffs),
                                  "request tA coeffs");
  req->w_coeffs = checked_malloc(kmsis * LEOPARD_RING_DEG * sizeof(*req->w_coeffs),
                                 "request w coeffs");
  req->rbit_const_tg_coeffs =
      checked_malloc(rbit_const_len * LEOPARD_RING_DEG *
                     sizeof(*req->rbit_const_tg_coeffs),
                     "request rbit const tg coeffs");
  req->rbit_const_h_coeffs =
      checked_malloc(rbit_const_len * LEOPARD_RING_DEG *
                     sizeof(*req->rbit_const_h_coeffs),
                     "request rbit const h coeffs");
  req->vlin_coeffs = checked_malloc(LEOPARD_RING_DEG * sizeof(*req->vlin_coeffs),
                                    "request vlin coeffs");
  req->tq_coeffs = checked_malloc(LEOPARD_RING_DEG * sizeof(*req->tq_coeffs),
                                  "request tq coeffs");
  req->vq_coeffs = checked_malloc(LEOPARD_RING_DEG * sizeof(*req->vq_coeffs),
                                  "request vq coeffs");
  polxvec_setzero(req->cr,cr_len);
  polxvec_setzero(req->cx,cx_len);
  polxvec_setzero(req->tA,kmsis);
  polxvec_setzero(req->w,kmsis);
  polxvec_setzero(req->rbit_const_tg,rbit_const_len);
  polxvec_setzero(req->rbit_const_h,rbit_const_len);
  polxvec_setzero(req->vlin,1);
  polxvec_setzero(req->tq,1);
  polxvec_setzero(req->vq,1);
  return req;
}

static leopard_nizk1_state *leopard_state_alloc(size_t rhat_len) {
  leopard_nizk1_state *st;

  st = checked_malloc(sizeof(*st),"state");
  memset(st,0,sizeof(*st));
  st->rhat_len = rhat_len;
  st->rhat = alloc_polxvec(rhat_len);
  st->rhat_coeffs = checked_malloc(rhat_len * LEOPARD_RING_DEG * sizeof(*st->rhat_coeffs),
                                   "state rhat coeffs");
  polxvec_setzero(st->rhat,rhat_len);
  return st;
}

static void leopard_nizk1_request_free(leopard_nizk1_request *req) {
  if(req == NULL)
    return;
  free(req->cr);
  free(req->cx);
  free(req->tA);
  free(req->w);
  free(req->rbit_const_tg);
  free(req->rbit_const_h);
  free(req->vlin);
  free(req->tq);
  free(req->vq);
  free(req->cr_coeffs);
  free(req->cx_coeffs);
  free(req->tA_coeffs);
  free(req->w_coeffs);
  free(req->rbit_const_tg_coeffs);
  free(req->rbit_const_h_coeffs);
  free(req->vlin_coeffs);
  free(req->tq_coeffs);
  free(req->vq_coeffs);
  free(req->proof);
  free(req);
}

static void leopard_nizk1_state_free(leopard_nizk1_state *st) {
  if(st == NULL)
    return;
  free(st->rhat);
  free(st->rhat_coeffs);
  free(st);
}

static void leopard_store_u64_le(uint8_t out[8], uint64_t x) {
  size_t i;

  for(i=0;i<8;i++)
    out[i] = (uint8_t)(x >> (8*i));
}

static void leopard_store_u32_le(uint8_t out[4], uint32_t x) {
  size_t i;

  for(i=0;i<4;i++)
    out[i] = (uint8_t)(x >> (8*i));
}

static void leopard_store_u16_le(uint8_t out[2], uint16_t x) {
  out[0] = (uint8_t)x;
  out[1] = (uint8_t)(x >> 8);
}

static uint64_t leopard_load_u64_le(const uint8_t in[8]) {
  size_t i;
  uint64_t x = 0;

  for(i=0;i<8;i++)
    x |= (uint64_t)in[i] << (8*i);
  return x;
}

static uint32_t leopard_load_u32_le(const uint8_t in[4]) {
  size_t i;
  uint32_t x = 0;

  for(i=0;i<4;i++)
    x |= (uint32_t)in[i] << (8*i);
  return x;
}

static uint16_t leopard_load_u16_le(const uint8_t in[2]) {
  return (uint16_t)in[0] | ((uint16_t)in[1] << 8);
}

static void leopard_store_i32_le(uint8_t out[4], int32_t x) {
  uint32_t u = (uint32_t)x;
  size_t i;

  for(i=0;i<4;i++)
    out[i] = (uint8_t)(u >> (8*i));
}

static int32_t leopard_load_i32_le(const uint8_t in[4]) {
  uint32_t x = 0;
  size_t i;

  for(i=0;i<4;i++)
    x |= (uint32_t)in[i] << (8*i);
  return (int32_t)x;
}

static int leopard_buf_write(uint8_t *buf, size_t len, size_t *off,
                             const void *src, size_t src_len) {
  if(*off > len || src_len > len - *off)
    return 1;
  memcpy(&buf[*off],src,src_len);
  *off += src_len;
  return 0;
}

static int leopard_buf_read(const uint8_t *buf, size_t len, size_t *off,
                            void *dst, size_t dst_len) {
  if(*off > len || dst_len > len - *off)
    return 1;
  memcpy(dst,&buf[*off],dst_len);
  *off += dst_len;
  return 0;
}

static unsigned int leopard_packed_bits_u64(uint64_t span) {
  unsigned int bits = 0;

  while(span != 0) {
    bits++;
    span >>= 1;
  }
  return bits;
}

static size_t leopard_packed_payload_bytes(size_t count, unsigned int bits) {
  if(bits == 0 || count == 0)
    return 0;
  return (count * (size_t)bits + 7) / 8;
}

static void leopard_pack_u64_offsets(uint8_t *dst, const uint64_t *vals,
                                     size_t count, unsigned int bits) {
  size_t i;
  size_t bitpos = 0;
  const size_t packed_bytes = leopard_packed_payload_bytes(count,bits);

  if(packed_bytes == 0)
    return;
  memset(dst,0,packed_bytes);
  for(i=0;i<count;i++) {
    unsigned int written = 0;

    while(written < bits) {
      const size_t byte_idx = bitpos >> 3;
      const unsigned int bit_idx = (unsigned int)(bitpos & 7);
      const unsigned int room = 8 - bit_idx;
      const unsigned int take = (bits - written < room) ? (bits - written) : room;
      const uint8_t mask = (take == 8) ? 0xffu : (uint8_t)((1u << take) - 1u);

      dst[byte_idx] |= (uint8_t)(((vals[i] >> written) & mask) << bit_idx);
      written += take;
      bitpos += take;
    }
  }
}

static void leopard_unpack_u64_offsets(uint64_t *vals, size_t count,
                                       unsigned int bits, const uint8_t *src) {
  size_t i;
  size_t bitpos = 0;

  if(bits == 0) {
    for(i=0;i<count;i++)
      vals[i] = 0;
    return;
  }
  for(i=0;i<count;i++) {
    uint64_t x = 0;
    unsigned int read = 0;

    while(read < bits) {
      const size_t byte_idx = bitpos >> 3;
      const unsigned int bit_idx = (unsigned int)(bitpos & 7);
      const unsigned int room = 8 - bit_idx;
      const unsigned int take = (bits - read < room) ? (bits - read) : room;
      const uint8_t mask = (take == 8) ? 0xffu : (uint8_t)((1u << take) - 1u);

      x |= (uint64_t)((src[byte_idx] >> bit_idx) & mask) << read;
      read += take;
      bitpos += take;
    }
    vals[i] = x;
  }
}

enum {
  LEOPARD_WIT_CHUNK_POLYS = 8,
  LEOPARD_PROJ_CHUNK_COEFFS = 256
};

#include "short_codec.h"
#ifndef LEOPARD_RICE_ENCODING
#define LEOPARD_RICE_ENCODING (LOGQ == 72)
#endif
#define LEOPARD_RICE_VERSION (UINT64_C(1) << 63)
#define LEOPARD_L2_WIDE_LENGTHS (UINT64_C(1) << 62)

static size_t leopard_serialized_proof_bytes(const proof *pi) {
  size_t proj_bytes = 0;
  size_t i;
  const size_t packed_bytes =
    (pi->cpp->u1len + pi->cpp->u2len + LIFTS) *
    LEOPARD_RING_DEG * LEOPARD_Q_BYTES;

  for(i=0;i<256;i += LEOPARD_PROJ_CHUNK_COEFFS) {
    const size_t chunk_len = MIN((size_t)LEOPARD_PROJ_CHUNK_COEFFS,(size_t)256 - i);
    int32_t pmin = pi->p[i];
    int32_t pmax = pi->p[i];
    unsigned int pbits;
    size_t j;

    for(j=1;j<chunk_len;j++) {
      if(pi->p[i + j] < pmin)
        pmin = pi->p[i + j];
      if(pi->p[i + j] > pmax)
        pmax = pi->p[i + j];
    }
    pbits = leopard_packed_bits_u64((uint64_t)((int64_t)pmax - (int64_t)pmin));
    proj_bytes += sizeof(int32_t) + sizeof(uint8_t) +
                  leopard_packed_payload_bytes(chunk_len,pbits);
  }

  return 2 * sizeof(uint8_t) +
         sizeof(uint64_t) +
         sizeof(uint32_t) +
         8 * sizeof(uint8_t) +
         2 * sizeof(uint16_t) +
         2 * pi->r * sizeof(uint32_t) +
         proj_bytes +
         packed_bytes;
}

static size_t leopard_serialized_witness_bytes(const witness *wt) {
  size_t i;
  size_t total = sizeof(uint64_t) + wt->r * sizeof(uint64_t);

  for(i=0;i<wt->r;i++) {
    size_t poly_off;

    for(poly_off=0;poly_off<wt->n[i];poly_off += LEOPARD_WIT_CHUNK_POLYS) {
      const size_t chunk_polys = MIN((size_t)LEOPARD_WIT_CHUNK_POLYS,wt->n[i] - poly_off);
      const size_t coeff_count = chunk_polys * LEOPARD_RING_DEG;
      int16_t minv = wt->s[i][poly_off].vec->c[0];
      int16_t maxv = wt->s[i][poly_off].vec->c[0];
      unsigned int bits;
      size_t j, k;

      if(LEOPARD_RICE_ENCODING) {
        int16_t coeffs[LEOPARD_WIT_CHUNK_POLYS * LEOPARD_RING_DEG], lo;
        unsigned mode;
        for(j=0;j<chunk_polys;j++)
          memcpy(coeffs+j*LEOPARD_RING_DEG,wt->s[i][poly_off+j].vec->c,
                 LEOPARD_RING_DEG*sizeof(int16_t));
        total += leopard_short_plan(coeffs,coeff_count,&mode,&lo);
        continue;
      }

      for(j=0;j<chunk_polys;j++) {
        for(k=0;k<LEOPARD_RING_DEG;k++) {
          const int16_t coeff = wt->s[i][poly_off + j].vec->c[k];

          if(coeff < minv)
            minv = coeff;
          if(coeff > maxv)
            maxv = coeff;
        }
      }
      bits = leopard_packed_bits_u64((uint64_t)((int32_t)maxv - (int32_t)minv));
      total += sizeof(int16_t) + sizeof(uint8_t) +
               leopard_packed_payload_bytes(coeff_count,bits);
    }
  }
  return total;
}

static size_t leopard_serialized_composite_bytes(const composite *p) {
  size_t i;
  size_t total = sizeof(uint64_t);

  for(i=0;i<p->l;i++)
    total += leopard_serialized_proof_bytes(p->pi[i]);
  total += leopard_serialized_witness_bytes(&p->owt);
  return total;
}

static int leopard_serialize_proof(uint8_t *buf, size_t len, size_t *off,
                                   const proof *pi) {
  uint8_t tmp64[8];
  uint8_t tmp32u[4];
  uint8_t tmp16u[2];
  uint8_t tmp32[4];
  size_t i;
  const polz *packed[3] = {pi->u1, pi->u2, pi->bb};
  const size_t packed_len[3] = {pi->cpp->u1len, pi->cpp->u2len, LIFTS};

  if(pi->r > 255 || pi->jlnonce > UINT32_MAX ||
     pi->cpp->f > 255 || pi->cpp->fu > 255 || pi->cpp->fg > 255 ||
     pi->cpp->b > 255 || pi->cpp->bu > 255 || pi->cpp->bg > 255 ||
     pi->cpp->kappa > 255 || pi->cpp->kappa1 > 255 ||
     pi->cpp->u1len > UINT16_MAX || pi->cpp->u2len > UINT16_MAX) {
    return 1;
  }
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->r},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->tail},sizeof(uint8_t)) != 0) return 1;
  leopard_store_u64_le(tmp64,pi->normsq);
  if(leopard_buf_write(buf,len,off,tmp64,sizeof(tmp64)) != 0) return 1;
  leopard_store_u32_le(tmp32u,(uint32_t)pi->jlnonce);
  if(leopard_buf_write(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) return 1;

  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->f},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->fu},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->fg},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->b},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->bu},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->bg},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->kappa},sizeof(uint8_t)) != 0) return 1;
  if(leopard_buf_write(buf,len,off,&(uint8_t){(uint8_t)pi->cpp->kappa1},sizeof(uint8_t)) != 0) return 1;
  leopard_store_u16_le(tmp16u,(uint16_t)pi->cpp->u1len);
  if(leopard_buf_write(buf,len,off,tmp16u,sizeof(tmp16u)) != 0) return 1;
  leopard_store_u16_le(tmp16u,(uint16_t)pi->cpp->u2len);
  if(leopard_buf_write(buf,len,off,tmp16u,sizeof(tmp16u)) != 0) return 1;

  for(i=0;i<pi->r;i++) {
    if(pi->n[i] > UINT32_MAX || pi->nu[i] > UINT32_MAX)
      return 1;
    leopard_store_u32_le(tmp32u,(uint32_t)pi->n[i]);
    if(leopard_buf_write(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) return 1;
  }
  for(i=0;i<pi->r;i++) {
    leopard_store_u32_le(tmp32u,(uint32_t)pi->nu[i]);
    if(leopard_buf_write(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) return 1;
  }
  {
    for(i=0;i<256;i += LEOPARD_PROJ_CHUNK_COEFFS) {
      const size_t chunk_len = MIN((size_t)LEOPARD_PROJ_CHUNK_COEFFS,(size_t)256 - i);
      int32_t pmin = pi->p[i];
      int32_t pmax = pi->p[i];
      uint8_t pbits8;
      size_t payload_len;
      uint64_t *offsets;
      uint8_t *payload;
      size_t j;

      for(j=1;j<chunk_len;j++) {
        if(pi->p[i + j] < pmin)
          pmin = pi->p[i + j];
        if(pi->p[i + j] > pmax)
          pmax = pi->p[i + j];
      }
      pbits8 = (uint8_t)leopard_packed_bits_u64((uint64_t)((int64_t)pmax - (int64_t)pmin));
      payload_len = leopard_packed_payload_bytes(chunk_len,pbits8);
      leopard_store_i32_le(tmp32,pmin);
      if(leopard_buf_write(buf,len,off,tmp32,sizeof(tmp32)) != 0) return 1;
      if(leopard_buf_write(buf,len,off,&pbits8,sizeof(pbits8)) != 0) return 1;
      if(payload_len == 0)
        continue;
      offsets = checked_malloc(chunk_len * sizeof(*offsets),"serialized proof offsets");
      payload = checked_malloc(payload_len,"serialized proof projection");
      for(j=0;j<chunk_len;j++)
        offsets[j] = (uint64_t)((int64_t)pi->p[i + j] - (int64_t)pmin);
      leopard_pack_u64_offsets(payload,offsets,chunk_len,pbits8);
      free(offsets);
      if(leopard_buf_write(buf,len,off,payload,payload_len) != 0) {
        free(payload);
        return 1;
      }
      free(payload);
    }
  }

  for(i=0;i<3;i++) {
    const size_t blob_len = packed_len[i] * LEOPARD_RING_DEG * LEOPARD_Q_BYTES;
    uint8_t *blob = checked_malloc(blob_len,"serialized proof blob");

    polzvec_bitpack(blob,packed[i],packed_len[i]);
    if(leopard_buf_write(buf,len,off,blob,blob_len) != 0) {
      free(blob);
      return 1;
    }
    free(blob);
  }

  return 0;
}

static int leopard_deserialize_proof(proof *pi, const uint8_t *buf, size_t len,
                                     size_t *off) {
  uint8_t tmp64[8];
  uint8_t tmp32u[4];
  uint8_t tmp16u[2];
  uint8_t tmp32[4];
  size_t i;
  const size_t packed_len_count = 3;
  size_t packed_len[3];
  polz *packed[3];

  memset(pi,0,sizeof(*pi));
  {
    uint8_t r8, tail8, small8;

    if(leopard_buf_read(buf,len,off,&r8,sizeof(r8)) != 0) return 1;
    pi->r = (size_t)r8;
    if(leopard_buf_read(buf,len,off,&tail8,sizeof(tail8)) != 0) return 1;
    pi->tail = (int)tail8;
    if(leopard_buf_read(buf,len,off,tmp64,sizeof(tmp64)) != 0) return 1;
    pi->normsq = leopard_load_u64_le(tmp64);
    if(leopard_buf_read(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) return 1;
    pi->jlnonce = (size_t)leopard_load_u32_le(tmp32u);

    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->f = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->fu = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->fg = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->b = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->bu = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->bg = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->kappa = (size_t)small8;
    if(leopard_buf_read(buf,len,off,&small8,sizeof(small8)) != 0) return 1;
    pi->cpp->kappa1 = (size_t)small8;
    if(leopard_buf_read(buf,len,off,tmp16u,sizeof(tmp16u)) != 0) return 1;
    pi->cpp->u1len = (size_t)leopard_load_u16_le(tmp16u);
    if(leopard_buf_read(buf,len,off,tmp16u,sizeof(tmp16u)) != 0) return 1;
    pi->cpp->u2len = (size_t)leopard_load_u16_le(tmp16u);
  }

  pi->n = _malloc(2 * pi->r * sizeof(size_t));
  pi->nu = &pi->n[pi->r];
  for(i=0;i<pi->r;i++) {
    if(leopard_buf_read(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) {
      free_proof(pi);
      return 1;
    }
    pi->n[i] = (size_t)leopard_load_u32_le(tmp32u);
    if(pi->n[i]>(1u<<20)) { free_proof(pi); return 1; }
  }
  for(i=0;i<pi->r;i++) {
    if(leopard_buf_read(buf,len,off,tmp32u,sizeof(tmp32u)) != 0) {
      free_proof(pi);
      return 1;
    }
    pi->nu[i] = (size_t)leopard_load_u32_le(tmp32u);
    if(pi->nu[i]>(1u<<20)) { free_proof(pi); return 1; }
  }
  {
    for(i=0;i<256;i += LEOPARD_PROJ_CHUNK_COEFFS) {
      const size_t chunk_len = MIN((size_t)LEOPARD_PROJ_CHUNK_COEFFS,(size_t)256 - i);
      int32_t pmin;
      uint8_t pbits8;
      size_t payload_len;
      uint8_t *payload = NULL;
      uint64_t offsets[LEOPARD_PROJ_CHUNK_COEFFS];
      size_t j;

      if(leopard_buf_read(buf,len,off,tmp32,sizeof(tmp32)) != 0) {
        free_proof(pi);
        return 1;
      }
      pmin = leopard_load_i32_le(tmp32);
      if(leopard_buf_read(buf,len,off,&pbits8,sizeof(pbits8)) != 0) {
        free_proof(pi);
        return 1;
      }
      if(pbits8>32) { free_proof(pi); return 1; }
      payload_len = leopard_packed_payload_bytes(chunk_len,pbits8);
      if(payload_len != 0) {
        payload = checked_malloc(payload_len,"serialized proof projection");
        if(leopard_buf_read(buf,len,off,payload,payload_len) != 0) {
          free(payload);
          free_proof(pi);
          return 1;
        }
        leopard_unpack_u64_offsets(offsets,chunk_len,pbits8,payload);
        free(payload);
      }
      else {
        leopard_unpack_u64_offsets(offsets,chunk_len,pbits8,NULL);
      }
      for(j=0;j<chunk_len;j++) {
        int64_t value=(int64_t)pmin+(int64_t)offsets[j];
        if(value>INT32_MAX) { free_proof(pi); return 1; }
        pi->p[i + j] = (int32_t)value;
      }
    }
  }

  pi->u1 = _aligned_alloc(64,(pi->cpp->u1len + pi->cpp->u2len + LIFTS) * sizeof(polz));
  pi->u2 = &pi->u1[pi->cpp->u1len];
  pi->bb = &pi->u2[pi->cpp->u2len];
  packed_len[0] = pi->cpp->u1len;
  packed_len[1] = pi->cpp->u2len;
  packed_len[2] = LIFTS;
  packed[0] = pi->u1;
  packed[1] = pi->u2;
  packed[2] = pi->bb;
  for(i=0;i<packed_len_count;i++) {
    size_t j;
    const size_t blob_len = packed_len[i] * LEOPARD_RING_DEG * LEOPARD_Q_BYTES;
    __attribute__((aligned(64)))
    uint8_t packed_poly[LEOPARD_RING_DEG * LEOPARD_Q_BYTES];

    if(*off > len || blob_len > len - *off) {
      free_proof(pi);
      return 1;
    }
    for(j=0;j<packed_len[i];j++) {
      memcpy(packed_poly,&buf[*off + j * LEOPARD_RING_DEG * LEOPARD_Q_BYTES],sizeof(packed_poly));
      polz_bitunpack(&packed[i][j],packed_poly);
    }
    *off += blob_len;
  }

  return 0;
}

static int leopard_serialize_witness(uint8_t *buf, size_t len, size_t *off,
                                     const witness *wt) {
  uint8_t tmp64[8];
  size_t i, j, k;

  leopard_store_u64_le(tmp64,wt->r);
  if(leopard_buf_write(buf,len,off,tmp64,sizeof(tmp64)) != 0) return 1;
  for(i=0;i<wt->r;i++) {
    leopard_store_u64_le(tmp64,wt->n[i]);
    if(leopard_buf_write(buf,len,off,tmp64,sizeof(tmp64)) != 0) return 1;
  }
  for(i=0;i<wt->r;i++) {
    size_t poly_off;

    for(poly_off=0;poly_off<wt->n[i];poly_off += LEOPARD_WIT_CHUNK_POLYS) {
      const size_t chunk_polys = MIN((size_t)LEOPARD_WIT_CHUNK_POLYS,wt->n[i] - poly_off);
      const size_t coeff_count = chunk_polys * LEOPARD_RING_DEG;
      int16_t minv = wt->s[i][poly_off].vec->c[0];
      int16_t maxv = wt->s[i][poly_off].vec->c[0];
      uint8_t bits8;
      size_t payload_len;
      uint64_t *offsets;
      uint8_t *payload;

      if(LEOPARD_RICE_ENCODING) {
        int16_t coeffs[LEOPARD_WIT_CHUNK_POLYS * LEOPARD_RING_DEG];
        size_t used;
        for(j=0;j<chunk_polys;j++)
          memcpy(coeffs+j*LEOPARD_RING_DEG,wt->s[i][poly_off+j].vec->c,
                 LEOPARD_RING_DEG*sizeof(int16_t));
        if(*off>len || leopard_short_encode(buf+*off,len-*off,&used,coeffs,coeff_count))
          return 1;
        *off+=used;
        continue;
      }

      for(j=0;j<chunk_polys;j++) {
        for(k=0;k<LEOPARD_RING_DEG;k++) {
          const int16_t coeff = wt->s[i][poly_off + j].vec->c[k];

          if(coeff < minv)
            minv = coeff;
          if(coeff > maxv)
            maxv = coeff;
        }
      }
      bits8 = (uint8_t)leopard_packed_bits_u64((uint64_t)((int32_t)maxv - (int32_t)minv));
      payload_len = leopard_packed_payload_bytes(coeff_count,bits8);
      if(leopard_buf_write(buf,len,off,&minv,sizeof(minv)) != 0) return 1;
      if(leopard_buf_write(buf,len,off,&bits8,sizeof(bits8)) != 0) return 1;
      if(payload_len == 0)
        continue;
      offsets = checked_malloc(coeff_count * sizeof(*offsets),"serialized witness offsets");
      payload = checked_malloc(payload_len,"serialized witness payload");
      for(j=0;j<chunk_polys;j++) {
        for(k=0;k<LEOPARD_RING_DEG;k++) {
          const size_t idx = j * LEOPARD_RING_DEG + k;

          offsets[idx] = (uint64_t)((int32_t)wt->s[i][poly_off + j].vec->c[k] - (int32_t)minv);
        }
      }
      leopard_pack_u64_offsets(payload,offsets,coeff_count,bits8);
      free(offsets);
      if(leopard_buf_write(buf,len,off,payload,payload_len) != 0) {
        free(payload);
        return 1;
      }
      free(payload);
    }
  }
  return 0;
}

static int leopard_deserialize_witness(witness *wt, const uint8_t *buf, size_t len,
                                       size_t *off, int rice) {
  uint8_t tmp64[8];
  size_t i, j, k;
  size_t *n;
  size_t total_polys=0;

  memset(wt,0,sizeof(*wt));
  if(leopard_buf_read(buf,len,off,tmp64,sizeof(tmp64)) != 0)
    return 1;
  wt->r = (size_t)leopard_load_u64_le(tmp64);
  if(wt->r==0 || wt->r>255 || wt->r>(len-*off)/sizeof(uint64_t)) {
    memset(wt,0,sizeof(*wt)); return 1;
  }
  n = checked_malloc(wt->r * sizeof(*n),"serialized witness sizes");
  for(i=0;i<wt->r;i++) {
    if(leopard_buf_read(buf,len,off,tmp64,sizeof(tmp64)) != 0) {
      free(n);
      return 1;
    }
    n[i] = (size_t)leopard_load_u64_le(tmp64);
    if(n[i]==0 || n[i]>UINT16_MAX-total_polys) {
      free(n); memset(wt,0,sizeof(*wt)); return 1;
    }
    total_polys+=n[i];
  }
  init_witness_raw(wt,wt->r,n);
  free(n);

  for(i=0;i<wt->r;i++) {
    size_t poly_off;

    for(poly_off=0;poly_off<wt->n[i];poly_off += LEOPARD_WIT_CHUNK_POLYS) {
      const size_t chunk_polys = MIN((size_t)LEOPARD_WIT_CHUNK_POLYS,wt->n[i] - poly_off);
      const size_t coeff_count = chunk_polys * LEOPARD_RING_DEG;
      int16_t minv;
      uint8_t bits8;
      size_t payload_len;
      uint8_t *payload = NULL;
      uint64_t *offsets = NULL;

      if(rice) {
        int16_t coeffs[LEOPARD_WIT_CHUNK_POLYS * LEOPARD_RING_DEG];
        size_t used;
        if(*off>len || leopard_short_decode(coeffs,coeff_count,buf+*off,len-*off,&used)) {
          free_witness(wt); memset(wt,0,sizeof(*wt)); return 1;
        }
        *off+=used;
        for(j=0;j<chunk_polys;j++)
          memcpy(wt->s[i][poly_off+j].vec->c,coeffs+j*LEOPARD_RING_DEG,
                 LEOPARD_RING_DEG*sizeof(int16_t));
        continue;
      }

      if(leopard_buf_read(buf,len,off,&minv,sizeof(minv)) != 0 ||
         leopard_buf_read(buf,len,off,&bits8,sizeof(bits8)) != 0) {
        free_witness(wt);
        return 1;
      }
      payload_len = leopard_packed_payload_bytes(coeff_count,bits8);
      if(bits8>16) { free_witness(wt); memset(wt,0,sizeof(*wt)); return 1; }
      if(payload_len != 0) {
        payload = checked_malloc(payload_len,"serialized witness payload");
        offsets = checked_malloc(coeff_count * sizeof(*offsets),"serialized witness offsets");
        if(leopard_buf_read(buf,len,off,payload,payload_len) != 0) {
          free(payload);
          free(offsets);
          free_witness(wt);
          return 1;
        }
        leopard_unpack_u64_offsets(offsets,coeff_count,bits8,payload);
        free(payload);
      }
      for(j=0;j<chunk_polys;j++) {
        for(k=0;k<LEOPARD_RING_DEG;k++) {
          const size_t idx = j * LEOPARD_RING_DEG + k;
          const uint64_t offv = (offsets == NULL) ? 0 : offsets[idx];

          wt->s[i][poly_off + j].vec->c[k] = (int16_t)((int32_t)minv + (int32_t)offv);
        }
      }
      free(offsets);
    }
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  }
  return 0;
}

static int leopard_serialize_composite(uint8_t **out, size_t *out_len,
                                       const composite *p) {
  uint8_t *buf;
  uint8_t tmp64[8];
  size_t off = 0;
  size_t i;
  const size_t total = leopard_serialized_composite_bytes(p);

  if(getenv("LEOPARD_CODEC_PROFILE")) {
    size_t uniform=0;
    for(size_t level=0;level<p->l;level++) {
      const proof *pi=p->pi[level];
      const size_t full=(pi->cpp->u1len+pi->cpp->u2len+LIFTS)*LEOPARD_RING_DEG*LEOPARD_Q_BYTES;
      uniform+=full;
      fprintf(stderr,"codec level=%zu uniform=%zu short_and_header=%zu\n",level,full,
              leopard_serialized_proof_bytes(pi)-full);
    }
    fprintf(stderr,"codec total=%zu uniform=%zu tail=%zu\n",total,uniform,
            leopard_serialized_witness_bytes(&p->owt));
  }

  buf = checked_malloc(total,"serialized composite");
  leopard_store_u64_le(tmp64,p->l | LEOPARD_L2_WIDE_LENGTHS |
                       (LEOPARD_RICE_ENCODING ? LEOPARD_RICE_VERSION : 0));
  if(leopard_buf_write(buf,total,&off,tmp64,sizeof(tmp64)) != 0) {
    free(buf);
    return 1;
  }
  for(i=0;i<p->l;i++) {
    if(leopard_serialize_proof(buf,total,&off,p->pi[i]) != 0) {
      free(buf);
      return 1;
    }
  }
  if(leopard_serialize_witness(buf,total,&off,&p->owt) != 0) {
    free(buf);
    return 1;
  }
  if(off != total) {
    free(buf);
    return 1;
  }
  *out = buf;
  *out_len = total;
  return 0;
}

static int leopard_deserialize_composite(composite *p, const uint8_t *buf,
                                         size_t len) {
  uint8_t tmp64[8];
  size_t off = 0;
  size_t i;

  memset(p,0,sizeof(*p));
  if(leopard_buf_read(buf,len,&off,tmp64,sizeof(tmp64)) != 0)
    return 1;
  uint64_t header=leopard_load_u64_le(tmp64);
  if(!(header & LEOPARD_L2_WIDE_LENGTHS)) return 1;
  int rice=(header & LEOPARD_RICE_VERSION)!=0;
  p->l = (size_t)(header & ~(LEOPARD_RICE_VERSION | LEOPARD_L2_WIDE_LENGTHS));
  if(p->l == 0 || p->l > 16) {
    p->l=0;
    return 1;
  }
  for(i=0;i<p->l;i++) {
    p->pi[i] = _malloc(sizeof(proof));
    if(leopard_deserialize_proof(p->pi[i],buf,len,&off) != 0) {
      p->l=i+1;
      free_composite(p);
      return 1;
    }
  }
  if(leopard_deserialize_witness(&p->owt,buf,len,&off,rice) != 0) {
    free_composite(p);
    return 1;
  }
  if(off != len) {
    free_composite(p);
    return 1;
  }
  return 0;
}

static void leopard_build_public_selector_material(selector_lnp22_ctx *ctx) {
  uint64_t nonce = 0;
  size_t row, j;
  poly *g1poly = alloc_polyvec(ctx->gdigits);

  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->m * ctx->gdigits >= 4096)
#endif
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits +
                                                            row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->gdigits >= 256)
#endif
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);

  free(g1poly);
}

static void leopard_derive_selector_chal_seed_from_cr(selector_lnp22_ctx *ctx,
                                                      const uint8_t tag_hash[LEOPARD_TAG_BYTES]) {
  shake128incctx st;

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-umc-cx-tag",19);
  absorb_bytes(&st,tag_hash,LEOPARD_TAG_BYTES);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,UMC_MODP);
  absorb_u64(&st,ctx->umc_n);
  absorb_u64(&st,ctx->umc_k);
  absorb_u64(&st,ctx->umc_v);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(ctx->selector_chal_seed,sizeof(ctx->selector_chal_seed),&st);
}

static void leopard_derive_ar_seed_from_cr(uint8_t out[16],
                                           const selector_lnp22_ctx *ctx,
                                           const uint8_t tag_hash[LEOPARD_TAG_BYTES]) {
  shake128incctx st;

  shake128_inc_init(&st);
  absorb_bytes(&st,"selector-umc-ar-tag",19);
  absorb_bytes(&st,tag_hash,LEOPARD_TAG_BYTES);
  absorb_u64(&st,LOGQ);
  absorb_u64(&st,ctx->Lchain);
  absorb_u64(&st,GADGET_BASE);
  absorb_u64(&st,UMC_MODP);
  absorb_u64(&st,ctx->umc_n);
  absorb_u64(&st,ctx->umc_k);
  absorb_u64(&st,ctx->umc_v);
  absorb_polxvec(&st,ctx->CX,ctx->cx_rows);
  shake128_inc_finalize(&st);
  shake128_inc_squeeze(out,16,&st);
}

static void leopard_fill_public_umc_material(selector_lnp22_ctx *ctx) {
  sample_uniform_modp_polxvec(ctx->Aumc,ctx->umc_n*ctx->umc_k,ctx->lnp_seed,NONCE_UMC_A_BASE,UMC_MODP);
  sample_uniform_modp_polxvec(ctx->Bumc,ctx->umc_v*ctx->umc_k,ctx->lnp_seed,NONCE_UMC_B_BASE,UMC_MODP);
  if(ctx->cache_a1open)
    fill_almostuniform_matrix(ctx->A1open,ctx->kmsis,ctx->m1,ctx->lnp_seed,NONCE_LNP_A1_BASE);
  fill_almostuniform_matrix(ctx->A2open,ctx->kmsis,ctx->m2prime,ctx->lnp_seed,NONCE_LNP_A2_BASE);
}

static void leopard_fill_public_selector_challenges(selector_lnp22_ctx *ctx) {
  size_t t;

  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_challenge(&ctx->kvec[t*ctx->gdigits],ctx->gdigits,ctx->selector_chal_seed,(uint64_t)t);
    polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,ctx->selector_chal_seed,32 + (uint64_t)t);
  }
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->selector_chal_seed,64);
  polxvec_challenge(ctx->alphavec,ctx->Lchain-1,ctx->selector_chal_seed,65);
}


#ifndef LEOPARD_VERIFY_FAST_FOLDED_COEFFS
#define LEOPARD_VERIFY_FAST_FOLDED_COEFFS 0
#endif

static void leopard_compute_folded_linear_coeffs_fast(selector_lnp22_ctx *ctx) {
  size_t idx = 0;
  size_t row, t, j, i;
  polx factor[1], cdelta[1];
  polx *tmp = alloc_polxvec(ctx->gdigits);
  polx *neg_monomials = NULL;
  const size_t c_idx = ctx->cx_rows;
  const size_t pack_idx = c_idx;
  const size_t cout_idx = pack_idx + ctx->rhat_len;
  const size_t top_idx = cout_idx + ctx->m;
  const size_t y_idx = top_idx + 1;
  const size_t tr_idx = y_idx + (ctx->Lchain-1);

  polxvec_setzero(ctx->lincoeff,ctx->m1);
  polxvec_setzero(ctx->lincoeff2,ctx->m2);
  polxvec_sprod(cdelta,ctx->kprime,ctx->gdelta,ctx->gdigits);

  for(row=0;row<ctx->umc_n;row++) {
    polxvec_polx_mul_add(&ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK]],
                         &ctx->linfold[idx],&ctx->Aumc[row*ctx->umc_k],ctx->umc_k);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_EUA_BLOCK] + row],
             &ctx->lincoeff[ctx->off[LNP22_UMC_EUA_BLOCK] + row],&ctx->linfold[idx]);
    polx_scale(factor,&ctx->linfold[idx],-UMC_MODP);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_VA_BLOCK] + row],
             &ctx->lincoeff[ctx->off[LNP22_UMC_VA_BLOCK] + row],factor);
    idx++;
  }

  for(t=0;t<ctx->umc_v;t++) {
    polx *msgcoeff = &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + t];
    polxvec_polx_mul_add(&ctx->lincoeff[ctx->off[LNP22_UMC_R_BLOCK]],
                         &ctx->linfold[idx],&ctx->Bumc[t*ctx->umc_k],ctx->umc_k);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_EUB_BLOCK] + t],
             &ctx->lincoeff[ctx->off[LNP22_UMC_EUB_BLOCK] + t],&ctx->linfold[idx]);
    polx_scale(factor,&ctx->linfold[idx],UMC_MSG_SCALE);
    polx_add(msgcoeff,msgcoeff,factor);
    polx_scale(factor,&ctx->linfold[idx],-UMC_MODP);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_UMC_VB_BLOCK] + t],
             &ctx->lincoeff[ctx->off[LNP22_UMC_VB_BLOCK] + t],factor);
    idx++;
  }

  if(LNP22_USE_RBIT_BINARY) {
    neg_monomials = alloc_polxvec(N);
    for(j=0;j<N;j++)
      polx_monomial(&neg_monomials[j],-1,(int)j);
    for(row=0;row<ctx->rhat_len;row++) {
      polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + row],
               &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK] + row],
               &ctx->linfold[idx]);
      for(j=0;j<N;j++) {
        polx_mul(factor,&ctx->linfold[idx],&neg_monomials[j]);
        polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + row*N + j],
                 &ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + row*N + j],
                 factor);
      }
      idx++;
    }
    free(neg_monomials);
  }

  for(row=0;row<ctx->m;row++) {
    polxvec_polx_mul_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK]],
                         &ctx->linfold[idx],&ctx->Arhash[row*ctx->rhat_len],ctx->rhat_len);
    for(j=0;j<ctx->gdigits/ctx->m;j++) {
      polx_scale_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],
                     &ctx->linfold[idx],pow3_u64(j));
    }
    idx++;
  }

  polxvec_polx_mul_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] +
                                     (ctx->Lchain-1)*ctx->gdigits],
                       &ctx->linfold[top_idx],ctx->kprime,ctx->gdigits);
  polx_mul(factor,&ctx->linfold[top_idx],cdelta);
  polx_sub(&ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + (ctx->Lchain-1)*SMALL_BLOCK],
           &ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + (ctx->Lchain-1)*SMALL_BLOCK],factor);

  for(t=0;t<ctx->Lchain-1;t++) {
    polx_mul(factor,&ctx->alphavec[t],&ctx->linfold[y_idx + t]);
    polxvec_polx_mul(tmp,factor,&ctx->kvec[t*ctx->gdigits],ctx->gdigits);
    polxvec_sub(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],
                &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],
                tmp,ctx->gdigits);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK] + t*SMALL_BLOCK],
             &ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK] + t*SMALL_BLOCK],&ctx->linfold[y_idx + t]);
    polx_scale(factor,&ctx->linfold[y_idx + t],Y_SPLIT_GAMMA);
    polx_add(&ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK] + t*SMALL_BLOCK],
             &ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK] + t*SMALL_BLOCK],factor);
  }

  for(t=0;t<ctx->Lchain-1;t++) {
    for(row=0;row<ctx->m;row++) {
      polx_mul(factor,&ctx->linfold[tr_idx + t],&ctx->rhovec[t*ctx->m + row]);
      for(j=0;j<ctx->gdigits/ctx->m;j++) {
        polx_scale_add(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] +
                                      t*ctx->gdigits + row*(ctx->gdigits/ctx->m) + j],
                       factor,pow3_u64(j));
      }
      polxvec_polx_mul(tmp,factor,&ctx->A0[row*ctx->gdigits],ctx->gdigits);
      polxvec_sub(&ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],
                  &ctx->lincoeff[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],
                  tmp,ctx->gdigits);
      polxvec_polx_mul(tmp,factor,&ctx->Dsel[row*ctx->gdigits],ctx->gdigits);
      polxvec_sub(&ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits],
                  &ctx->lincoeff[ctx->off[LNP22_Z_BLOCK] + t*ctx->gdigits],
                  tmp,ctx->gdigits);
    }
  }

  idx = tr_idx + (ctx->Lchain-1);
  if(ctx->rbit_const_lambda != 0) {
    polx phi[1];
    for(row=0;row<ctx->rbit_const_lambda;row++) {
      for(i=0;i<ctx->rhat_bits_len;i++) {
        sample_rbit_const_phi(phi,ctx,row,i);
        polx_mul(factor,&ctx->linfold[idx],phi);
        polx_add(&ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],
                 &ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK] + i],factor);
      }
      if(LNP22_PROVE_X_BITS) {
        for(i=0;i<ctx->x_len;i++) {
          sample_x_const_phi(phi,ctx,row,i);
          polx_mul(factor,&ctx->linfold[idx],phi);
          polx_add(&ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + i],
                   &ctx->lincoeff[ctx->off[LNP22_X_BLOCK] + i],factor);
        }
      }
      polxvec_polx_mul(tmp,&ctx->linfold[idx],&ctx->rbit_const_Bg[row*ctx->m2],ctx->m2);
      polxvec_sub(ctx->lincoeff2,ctx->lincoeff2,tmp,ctx->m2);
      idx++;
    }
  }

  free(tmp);
}

static void leopard_compute_arhash(selector_lnp22_ctx *ctx,
                                   const uint8_t tag_hash[LEOPARD_TAG_BYTES]) {
  uint8_t ar_seed[16];

  leopard_derive_ar_seed_from_cr(ar_seed,ctx,tag_hash);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,ar_seed,NONCE_AR_BASE);
}

static void leopard_extract_bxt(polx *out, const selector_lnp22_ctx *ctx) {
  size_t row, j;

  for(row=0;row<ctx->m;row++) {
    polxvec_setzero(&out[row],1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(&out[row],
                     &ctx->s1[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],
                     pow3_u64(j));
  }
}

static void leopard_compute_cout(selector_lnp22_ctx *ctx,
                                 const uint8_t tag_hash[LEOPARD_TAG_BYTES]) {
  size_t row, j;

  leopard_compute_arhash(ctx,tag_hash);
  for(row=0;row<ctx->m;row++) {
    polx prod[1], gb[1];
    polxvec_sprod(prod,&ctx->Arhash[row*ctx->rhat_len],&ctx->s1[ctx->off[LNP22_RHAT_BLOCK]],ctx->rhat_len);
    polxvec_setzero(gb,1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(gb,&ctx->s1[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],pow3_u64(j));
    polx_add(&ctx->Cout[row],prod,gb);
  }
}

static void leopard_compute_cout_from_loaded_arhash(selector_lnp22_ctx *ctx) {
  size_t row, j;

  for(row=0;row<ctx->m;row++) {
    polx prod[1], gb[1];
    polxvec_sprod(prod,&ctx->Arhash[row*ctx->rhat_len],
                  &ctx->s1[ctx->off[LNP22_RHAT_BLOCK]],ctx->rhat_len);
    polxvec_setzero(gb,1);
    for(j=0;j<ctx->gdigits/ctx->m;j++)
      polx_scale_add(gb,
                     &ctx->s1[ctx->off[LNP22_B_BLOCK] + row*(ctx->gdigits/ctx->m) + j],
                     pow3_u64(j));
    polx_add(&ctx->Cout[row],prod,gb);
  }
}

static void leopard_load_preprocessed_slot_from_state(selector_lnp22_ctx *ctx,
                                                      const leopard_pre_state_data *pre,
                                                      size_t slot_j) {
  ctx->selected_j = slot_j;
  import_poly_coeffs(ctx->witness[LNP22_RHAT_BLOCK],
                     &pre->rhat_coeffs[slot_j * pre->rhat_rows * LEOPARD_RING_DEG],
                     pre->rhat_rows);
  import_poly_coeffs(ctx->witness[LNP22_UMC_EUA_BLOCK],
                     &pre->eua_coeffs[slot_j * pre->eua_rows * LEOPARD_RING_DEG],
                     pre->eua_rows);
  import_poly_coeffs(ctx->witness[LNP22_UMC_EUB_BLOCK],
                     &pre->eub_coeffs[slot_j * pre->eub_rows * LEOPARD_RING_DEG],
                     pre->eub_rows);
  import_poly_coeffs(ctx->witness[LNP22_UMC_R_BLOCK],
                     &pre->r_coeffs[slot_j * pre->r_rows * LEOPARD_RING_DEG],
                     pre->r_rows);
  import_poly_coeffs(ctx->witness[LNP22_UMC_VA_BLOCK],
                     &pre->va_coeffs[slot_j * pre->va_rows * LEOPARD_RING_DEG],
                     pre->va_rows);
  import_poly_coeffs(ctx->witness[LNP22_UMC_VB_BLOCK],
                     &pre->vb_coeffs[slot_j * pre->vb_rows * LEOPARD_RING_DEG],
                     pre->vb_rows);
  sync_witness_block_to_s1(ctx,LNP22_RHAT_BLOCK);
  populate_rhat_bits(ctx);
  sync_witness_block_to_s1(ctx,LNP22_UMC_EUA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_EUB_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_R_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_VA_BLOCK);
  sync_witness_block_to_s1(ctx,LNP22_UMC_VB_BLOCK);
  import_polx_coeffs(ctx->CX,
                     &pre->cr_coeffs[slot_j * pre->cr_rows * LEOPARD_RING_DEG],
                     pre->cr_rows);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,
                        &pre->ar_seeds[16 * slot_j],NONCE_AR_BASE);
  update_umc_v_stats(ctx);
}

#ifndef LEOPARD_FAST_WITNESS_PREPARE
#define LEOPARD_FAST_WITNESS_PREPARE 1
#endif

#if LNP22_DIRECT_CHAIN
static void leopard_prepare_decompose(poly *out, const polx *in, size_t digits) {
#if LEOPARD_FAST_WITNESS_PREPARE
  leopard_public_decompose(out,in,digits);
#else
  gadget_decompose_balanced3(out,in,digits);
#endif
}

static void leopard_prepare_selector_direct(
    selector_lnp22_ctx *ctx,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS]) {
  poly *g1poly = alloc_polyvec(ctx->gdigits);
  uint64_t nonce = 0;
  size_t row, t, j;
  struct timespec t0, t1;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);
#if LEOPARD_FAST_WITNESS_PREPARE
  leopard_public_almostuniform(ctx->A0,ctx->m*ctx->gdigits,
                              ctx->statement_seed,nonce,leopard_public_threads());
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  leopard_public_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,
                              ctx->statement_seed,nonce,leopard_public_threads());
#else
  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,
                        ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,
                        ctx->statement_seed,nonce);
#endif
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->m;row++) {
    leopard_prepare_decompose(
        &ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits +
                                      row*(ctx->gdigits/ctx->m)],
        &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    leopard_prepare_decompose(
        &g1poly[row*(ctx->gdigits/ctx->m)],
        &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(
      ctx->g0,
      &ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],
      ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_selector_material_sec =
      leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  for(t=0;t<ctx->Lchain;t++)
    poly_set_const(
        &ctx->witness[LNP22_X_BLOCK][selector_x_pos(ctx,t)],
        xbits[t] & 1);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_X_BLOCK]],
                      ctx->witness[LNP22_X_BLOCK],ctx->x_len);

  if(pre != NULL) {
    leopard_load_preprocessed_slot_from_state(ctx,pre,slot_j);
  } else {
    ctx->selected_j = slot_j;
    for(t=0;t<ctx->rhat_len;t++)
      sample_binary_poly(&ctx->witness[LNP22_RHAT_BLOCK][t],
                         ctx->lnp_seed,NONCE_RHAT_BASE + t);
    polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_RHAT_BLOCK]],
                        ctx->witness[LNP22_RHAT_BLOCK],ctx->rhat_len);
    populate_rhat_bits(ctx);
    compute_umc_commitment(ctx);
    leopard_compute_arhash(ctx,tag_hash);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_slot_load_sec =
      leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(ctx->witness[LNP22_X_BLOCK]
         [selector_x_pos(ctx,ctx->Lchain-1)].vec->c[0]) {
    polyvec_copy(
        &ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],
        g1poly,ctx->gdigits);
    polxvec_copy(
        &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                 (ctx->Lchain-1)*ctx->gdigits],
        ctx->g1,ctx->gdigits);
  } else {
    polxvec_copy(
        &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                 (ctx->Lchain-1)*ctx->gdigits],
        ctx->g0,ctx->gdigits);
  }

#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel private(row)
  {
    for(size_t tt=ctx->Lchain-1;tt-- > 0;) {
#pragma omp for schedule(static)
      for(row=0;row<ctx->m;row++) {
        polx rhsrow[1];

#if LEOPARD_FAST_WITNESS_PREPARE
        const polx *selected = xbits[tt] ? ctx->A1sel : ctx->A0;
        polxvec_sprod(rhsrow,&selected[row*ctx->gdigits],
                      &ctx->s1[ctx->off[LNP22_B_BLOCK] + (tt+1)*ctx->gdigits],
                      ctx->gdigits);
#else
        polxvec_sprod(
            rhsrow,&ctx->A0[row*ctx->gdigits],
            &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                     (tt+1)*ctx->gdigits],
            ctx->gdigits);
        if(ctx->witness[LNP22_X_BLOCK]
              [selector_x_pos(ctx,tt)].vec->c[0])
          polxvec_sprod_add(
              rhsrow,&ctx->Dsel[row*ctx->gdigits],
              &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                       (tt+1)*ctx->gdigits],
              ctx->gdigits);
#endif
        leopard_prepare_decompose(
            &ctx->witness[LNP22_B_BLOCK]
                         [tt*ctx->gdigits +
                          row*(ctx->gdigits/ctx->m)],
            rhsrow,ctx->gdigits/ctx->m);
      }
#pragma omp for schedule(static)
      for(j=0;j<ctx->gdigits;j++)
        polx_frompoly(
            &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                     tt*ctx->gdigits + j],
            &ctx->witness[LNP22_B_BLOCK][tt*ctx->gdigits + j]);
    }
  }
#else
  for(t=ctx->Lchain-1;t-- > 0;) {
    for(row=0;row<ctx->m;row++) {
      polx rhsrow[1];

#if LEOPARD_FAST_WITNESS_PREPARE
      const polx *selected = xbits[t] ? ctx->A1sel : ctx->A0;
      polxvec_sprod(rhsrow,&selected[row*ctx->gdigits],
                    &ctx->s1[ctx->off[LNP22_B_BLOCK] + (t+1)*ctx->gdigits],
                    ctx->gdigits);
#else
      polxvec_sprod(
          rhsrow,&ctx->A0[row*ctx->gdigits],
          &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                   (t+1)*ctx->gdigits],
          ctx->gdigits);
      if(ctx->witness[LNP22_X_BLOCK]
            [selector_x_pos(ctx,t)].vec->c[0])
        polxvec_sprod_add(
            rhsrow,&ctx->Dsel[row*ctx->gdigits],
            &ctx->s1[ctx->off[LNP22_B_BLOCK] +
                     (t+1)*ctx->gdigits],
            ctx->gdigits);
#endif
      leopard_prepare_decompose(
          &ctx->witness[LNP22_B_BLOCK]
                       [t*ctx->gdigits +
                        row*(ctx->gdigits/ctx->m)],
          rhsrow,ctx->gdigits/ctx->m);
    }
    polxvec_frompolyvec(
        &ctx->s1[ctx->off[LNP22_B_BLOCK] + t*ctx->gdigits],
        &ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits],
        ctx->gdigits);
  }
#endif
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_hash_eval_sec =
      leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(pre != NULL)
    leopard_compute_cout_from_loaded_arhash(ctx);
  else
    leopard_compute_cout(ctx,tag_hash);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_cx_sec =
      leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.protocol_transcript_sec =
      g_last_request_profile.protocol_selector_material_sec +
      g_last_request_profile.protocol_slot_load_sec +
      g_last_request_profile.protocol_hash_eval_sec +
      g_last_request_profile.protocol_cx_sec;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  derive_selector_chal_seed_from_public(ctx);
  polxvec_challenge(ctx->kprime,ctx->gdigits,
                    ctx->selector_chal_seed,64);
  leopard_fill_public_umc_material(ctx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.nizk_prepare_sec =
      leopard_elapsed_sec(&t0,&t1);
  free(g1poly);
}
#endif

static void leopard_prepare_selector_lnp22_preprocessed_explicit(
    selector_lnp22_ctx *ctx,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS]) {
#if LNP22_DIRECT_CHAIN
  leopard_prepare_selector_direct(ctx,pre,slot_j,tag_hash,xbits);
  return;
#else
  polx *sx = alloc_polxvec(ctx->gdigits);
  poly *g1poly = alloc_polyvec(ctx->gdigits);
  uint64_t nonce = 0;
  size_t row, t, j;
  polx yx[1], scaled[1];
  struct timespec t0, t1;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);
  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->m * ctx->gdigits >= 4096)
#endif
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits +
                                                            row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->gdigits >= 256)
#endif
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_selector_material_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  for(t=0;t<ctx->Lchain;t++)
    poly_set_const(&ctx->witness[LNP22_X_BLOCK][t*SMALL_BLOCK],xbits[t] & 1);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_X_BLOCK]],ctx->witness[LNP22_X_BLOCK],ctx->x_len);

  leopard_load_preprocessed_slot_from_state(ctx,pre,slot_j);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_slot_load_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(ctx->witness[LNP22_X_BLOCK][(ctx->Lchain-1)*SMALL_BLOCK].vec->c[0])
    polyvec_copy(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],g1poly,ctx->gdigits);

  for(t=ctx->Lchain-1;t-- > 0;) {
    const int xbit = ctx->witness[LNP22_X_BLOCK][t*SMALL_BLOCK].vec->c[0] != 0;
    const polx *Amat = xbit ? ctx->A1sel : ctx->A0;

    if(xbit)
      polyvec_copy(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],
                   &ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    else
      polyvec_setzero(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],ctx->gdigits);

    polxvec_frompolyvec(sx,&ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
#if defined(_OPENMP) && LEOPARD_PUBLIC_ROW_PARALLEL
#pragma omp parallel for schedule(static) if(ctx->m >= LEOPARD_PUBLIC_PARALLEL_MIN_ROWS)
#endif
    for(row=0;row<ctx->m;row++) {
      polx rhs[1];
      polxvec_sprod(rhs,&Amat[row*ctx->gdigits],sx,ctx->gdigits);
      gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits +
                                                              row*(ctx->gdigits/ctx->m)],
                                 rhs,ctx->gdigits/ctx->m);
    }
  }

  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_B_BLOCK]],ctx->witness[LNP22_B_BLOCK],ctx->b_len);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Z_BLOCK]],ctx->witness[LNP22_Z_BLOCK],ctx->z_len);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_hash_eval_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  leopard_compute_cout_from_loaded_arhash(ctx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_cx_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.protocol_transcript_sec =
      g_last_request_profile.protocol_selector_material_sec +
      g_last_request_profile.protocol_slot_load_sec +
      g_last_request_profile.protocol_hash_eval_sec +
      g_last_request_profile.protocol_cx_sec;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  derive_selector_chal_seed_from_public(ctx);
  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_challenge(&ctx->kvec[t*ctx->gdigits],ctx->gdigits,ctx->selector_chal_seed,(uint64_t)t);
    polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,ctx->selector_chal_seed,32 + (uint64_t)t);
  }
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->selector_chal_seed,64);
  polxvec_challenge(ctx->alphavec,ctx->Lchain-1,ctx->selector_chal_seed,65);

  for(t=0;t<ctx->Lchain-1;t++) {
    polz yz[1];
    zz coeff;

    polxvec_frompolyvec(&sx[0],&ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    polxvec_sprod(yx,&ctx->kvec[t*ctx->gdigits],&sx[0],ctx->gdigits);
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
  leopard_fill_public_umc_material(ctx);
  compute_quad_g1_coeffs(ctx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.nizk_prepare_sec = leopard_elapsed_sec(&t0,&t1);

  free(sx);
  free(g1poly);
#endif
}

static int leopard_prepare_selector_public_ctx_snow(
    selector_lnp22_ctx *ctx,
    const leopard_pp *pp,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    size_t slot_j,
    const polx *cr,
    const polx *cx,
    const polx *rbit_const_tg,
    const polx *rbit_const_h,
    size_t rbit_const_len,
    const uint8_t ar_seed[16], const polx *tA) {
#if LNP22_DIRECT_CHAIN
  selector_ctx_init(ctx,pp->nizk1_seed);
  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);
  prepare_selector_lnp22_verifier_static(ctx);
  ctx->selected_j = slot_j;
  polxvec_copy(ctx->CX,cr,ctx->cx_rows);
  polxvec_copy(ctx->Cout,cx,ctx->m);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,
                        ar_seed,NONCE_AR_BASE);

  derive_selector_chal_seed_from_public(ctx);
  bind_prechallenge_commitment(ctx,tA);
  if(rbit_const_len != ctx->rbit_const_lambda)
    return 1;
  if(ctx->rbit_const_lambda != 0) {
    polxvec_almostuniform(
        ctx->rbit_const_Bg,
        ctx->rbit_const_lambda * ctx->m2,
        ctx->selector_chal_seed,NONCE_RBIT_CONST_BG_BASE);
    polxvec_copy(ctx->rbit_const_tg,rbit_const_tg,
                 ctx->rbit_const_lambda);
    polxvec_copy(ctx->rbit_const_h,rbit_const_h,
                 ctx->rbit_const_lambda);
    if(validate_rbit_const_transcript(ctx))
      return 1;
    derive_rbit_const_seed(ctx);
  }
  compute_linear_constants(ctx);
  compute_linear_fold_constants(ctx);
  compute_native_binary_fold(ctx);
  compute_folded_linear_coeffs(ctx);
  return 0;
#else
  poly *g0poly;
  poly *g1poly;
  uint64_t nonce = 0;
  size_t row, j, t;

  selector_ctx_init(ctx,pp->nizk1_seed);
  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);
  ctx->selected_j = slot_j;

  g0poly = alloc_polyvec(ctx->gdigits);
  g1poly = alloc_polyvec(ctx->gdigits);
  if(!g0poly || !g1poly) {
    fprintf(stderr,"ERROR: verifier public allocation failure\n");
    exit(1);
  }

  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&g0poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,g0poly,ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);

  sample_umc_matrices(ctx);
  polxvec_copy(ctx->CX,cr,ctx->cx_rows);
  polxvec_copy(ctx->Cout,cx,ctx->m);
  polxvec_almostuniform(ctx->Arhash,ctx->rhat_len*ctx->m,ar_seed,NONCE_AR_BASE);

  derive_selector_chal_seed_from_public(ctx);
#if defined(_OPENMP) && LNP22_PREP_ROW_PARALLEL
#pragma omp parallel for schedule(static)
#endif
  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_challenge(&ctx->kvec[t*ctx->gdigits],ctx->gdigits,ctx->selector_chal_seed,(uint64_t)t);
    polxvec_challenge(&ctx->rhovec[t*ctx->m],ctx->m,ctx->selector_chal_seed,32 + (uint64_t)t);
  }
  polxvec_challenge(ctx->kprime,ctx->gdigits,ctx->selector_chal_seed,64);
  polxvec_challenge(ctx->alphavec,ctx->Lchain-1,ctx->selector_chal_seed,65);

  if(ctx->cache_a1open)
    fill_almostuniform_matrix(ctx->A1open,ctx->kmsis,ctx->m1,ctx->lnp_seed,NONCE_LNP_A1_BASE);
  fill_almostuniform_matrix(ctx->A2open,ctx->kmsis,ctx->m2prime,ctx->lnp_seed,NONCE_LNP_A2_BASE);

  if(rbit_const_len == ctx->rbit_const_lambda && ctx->rbit_const_lambda != 0) {
    polxvec_almostuniform(ctx->rbit_const_Bg,
                          ctx->rbit_const_lambda * ctx->m2,
                          ctx->selector_chal_seed,
                          NONCE_RBIT_CONST_BG_BASE);
    polxvec_copy(ctx->rbit_const_tg,rbit_const_tg,ctx->rbit_const_lambda);
    polxvec_copy(ctx->rbit_const_h,rbit_const_h,ctx->rbit_const_lambda);
    derive_rbit_const_seed(ctx);
  }

  compute_linear_constants(ctx);
  compute_linear_fold_constants(ctx);
  compute_folded_linear_coeffs(ctx);

  free(g0poly);
  free(g1poly);
  return 0;
#endif
}

static void leopard_prepare_selector_lnp22_explicit(selector_lnp22_ctx *ctx,
                                                    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                                    const uint8_t xbits[LEOPARD_X_BITS]) {
  polx *sx = alloc_polxvec(2*ctx->gdigits);
  poly *g1poly = alloc_polyvec(ctx->gdigits);
  uint64_t nonce = 0;
  size_t row, t, j;
  polx yx[1], scaled[1];

  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);

  polxvec_almostuniform(ctx->A0,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  nonce += (ctx->m*ctx->gdigits + 31)/32;
  polxvec_almostuniform(ctx->A1sel,ctx->m*ctx->gdigits,ctx->statement_seed,nonce);
  for(j=0;j<ctx->m*ctx->gdigits;j++)
    polx_sub(&ctx->Dsel[j],&ctx->A1sel[j],&ctx->A0[j]);

  for(row=0;row<ctx->m;row++) {
    gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits +
                                                            row*(ctx->gdigits/ctx->m)],
                               &ctx->A0[row*ctx->gdigits],ctx->gdigits/ctx->m);
    gadget_decompose_balanced3(&g1poly[row*(ctx->gdigits/ctx->m)],
                               &ctx->A1sel[row*ctx->gdigits],ctx->gdigits/ctx->m);
  }
  polxvec_frompolyvec(ctx->g0,&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],ctx->gdigits);
  polxvec_frompolyvec(ctx->g1,g1poly,ctx->gdigits);
  for(j=0;j<ctx->gdigits;j++)
    polx_sub(&ctx->gdelta[j],&ctx->g1[j],&ctx->g0[j]);

  for(t=0;t<ctx->Lchain;t++)
    poly_set_const(&ctx->witness[LNP22_X_BLOCK][t*SMALL_BLOCK],xbits[t] & 1);
  for(t=0;t<ctx->rhat_len;t++)
    sample_binary_poly(&ctx->witness[LNP22_RHAT_BLOCK][t],ctx->lnp_seed,NONCE_RHAT_BASE + t);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_X_BLOCK]],ctx->witness[LNP22_X_BLOCK],ctx->x_len);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_RHAT_BLOCK]],ctx->witness[LNP22_RHAT_BLOCK],ctx->rhat_len);
  populate_rhat_bits(ctx);
  compute_umc_commitment(ctx);
  leopard_derive_selector_chal_seed_from_cr(ctx,tag_hash);
  leopard_fill_public_selector_challenges(ctx);

  if(ctx->witness[LNP22_X_BLOCK][(ctx->Lchain-1)*SMALL_BLOCK].vec->c[0])
    polyvec_copy(&ctx->witness[LNP22_B_BLOCK][(ctx->Lchain-1)*ctx->gdigits],g1poly,ctx->gdigits);

  for(t=ctx->Lchain-1;t-- > 0;) {
    if(ctx->witness[LNP22_X_BLOCK][t*SMALL_BLOCK].vec->c[0])
      polyvec_copy(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],
                   &ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    else
      polyvec_setzero(&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],ctx->gdigits);

    polxvec_frompolyvec(&sx[0],&ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    polxvec_frompolyvec(&sx[ctx->gdigits],&ctx->witness[LNP22_Z_BLOCK][t*ctx->gdigits],ctx->gdigits);
    for(row=0;row<ctx->m;row++) {
      polx rhs[1];
      polxvec_sprod(rhs,&ctx->A0[row*ctx->gdigits],&sx[0],ctx->gdigits);
      polxvec_sprod_add(rhs,&ctx->Dsel[row*ctx->gdigits],&sx[ctx->gdigits],ctx->gdigits);
      gadget_decompose_balanced3(&ctx->witness[LNP22_B_BLOCK][t*ctx->gdigits +
                                                              row*(ctx->gdigits/ctx->m)],
                                 rhs,ctx->gdigits/ctx->m);
    }
  }

  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_B_BLOCK]],ctx->witness[LNP22_B_BLOCK],ctx->b_len);

  for(t=0;t<ctx->Lchain-1;t++) {
    polz yz[1];
    zz coeff;

    polxvec_frompolyvec(&sx[0],&ctx->witness[LNP22_B_BLOCK][(t+1)*ctx->gdigits],ctx->gdigits);
    polxvec_sprod(yx,&ctx->kvec[t*ctx->gdigits],&sx[0],ctx->gdigits);
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

  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Z_BLOCK]],ctx->witness[LNP22_Z_BLOCK],ctx->z_len);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Y0_BLOCK]],ctx->witness[LNP22_Y0_BLOCK],ctx->y_len);
  polxvec_frompolyvec(&ctx->s1[ctx->off[LNP22_Y1_BLOCK]],ctx->witness[LNP22_Y1_BLOCK],ctx->y_len);
  leopard_compute_cout(ctx,tag_hash);
  compute_quad_g1_coeffs(ctx);

  free(sx);
  free(g1poly);
}

static void leopard_prepare_selector_public_ctx(selector_lnp22_ctx *ctx,
                                                const leopard_pp *pp,
                                                const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                                const polx *cr,
                                                const polx *cx) {
  selector_ctx_init(ctx,pp->nizk1_seed);
  memcpy(ctx->statement_seed,tag_hash,LEOPARD_TAG_BYTES);
  leopard_build_public_selector_material(ctx);
  leopard_fill_public_umc_material(ctx);
  polxvec_copy(ctx->CX,cr,ctx->cx_rows);
  polxvec_copy(ctx->Cout,cx,ctx->m);
  leopard_derive_selector_chal_seed_from_cr(ctx,tag_hash);
  leopard_fill_public_selector_challenges(ctx);
  leopard_compute_arhash(ctx,tag_hash);
  compute_linear_constants(ctx);
  compute_linear_fold_constants(ctx);
  compute_folded_linear_coeffs(ctx);
}

#if !LNP22_DIRECT_CHAIN
static int leopard_generate_request_instance(outer_instance *inst,
                                             const leopard_pp *pp,
                                             const uint8_t tag_hash[LEOPARD_TAG_BYTES],
                                             const uint8_t xbits[LEOPARD_X_BITS]) {
  const size_t batch_rounds = LNP22_ENABLE_REJECTION ? LNP22_BATCH_ROUNDS : 1;
  size_t batch_workers = 1;
  lnp22_round_ws *batch_ws = NULL;
  int ok = 1;
  struct timespec t0, t1;

  memset(inst,0,sizeof(*inst));
  selector_ctx_init(&inst->ctx,pp->nizk1_seed);
  g_last_request_profile.a1open_cache = inst->ctx.cache_a1open;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  leopard_prepare_selector_lnp22_explicit(&inst->ctx,tag_hash,xbits);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_prepare_sec = leopard_elapsed_sec(&t0,&t1);
  if(first_linear_relation_failure(&inst->ctx,inst->ctx.s1)) {
    fprintf(stderr,"ERROR: request witness does not satisfy intended relations\n");
    return 1;
  }

  inst->s2 = alloc_polxvec(inst->ctx.m2);
  inst->tA = alloc_polxvec(inst->ctx.kmsis);
  inst->a3 = alloc_polxvec(inst->ctx.m2);
  inst->w = alloc_polxvec(inst->ctx.kmsis);
  inst->z1 = alloc_polxvec(inst->ctx.m1);
  inst->z2 = alloc_polxvec(inst->ctx.m2);
  inst->cs1 = alloc_polxvec(inst->ctx.m1);
  inst->cs2 = alloc_polxvec(inst->ctx.m2);
  inst->y1 = alloc_polxvec(inst->ctx.m1);
  inst->y2 = alloc_polxvec(inst->ctx.m2);
  inst->z1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->z1coeffs));
  inst->z2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->z2coeffs));
  inst->cs1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->cs1coeffs));
  inst->cs2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->cs2coeffs));
  if(batch_rounds > 1) {
    size_t bi;

#ifdef _OPENMP
    batch_workers = (size_t)omp_get_max_threads();
#endif
    if(batch_workers > batch_rounds)
      batch_workers = batch_rounds;
    batch_ws = _malloc(batch_rounds * sizeof(*batch_ws));
    for(bi=0;bi<batch_rounds;bi++)
      round_ws_init(&batch_ws[bi],&inst->ctx);
  }

  clock_gettime(CLOCK_MONOTONIC,&t0);
  compute_public_quad_rhs(inst->uquad_public,&inst->ctx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.nizk_inner_setup_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  sample_ternary_polxvec(inst->s2,inst->ctx.m2,inst->ctx.lnp_seed,NONCE_S2_BASE);
  compute_a3(inst->a3,&inst->ctx);
  inner_product_a3(inst->a3s2,inst->a3,inst->s2,inst->ctx.m2);
  {
    struct timespec lt0, lt1;
    clock_gettime(CLOCK_MONOTONIC,&lt0);
    if(prepare_linear_fold_for_proof(&inst->ctx,inst->ufold_public,inst->s2))
      return 1;
    clock_gettime(CLOCK_MONOTONIC,&lt1);
    g_last_request_profile.nizk_prepare_sec += leopard_elapsed_sec(&lt0,&lt1);
  }
  compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);

  if(batch_rounds <= 1 || batch_ws == NULL) {
    for(inst->rejection_rounds = 0;;inst->rejection_rounds++) {
      size_t i;
      int rej1, rej2;

      sample_gaussian_polxvec(inst->y1,inst->z1coeffs,inst->ctx.m1,LNP22_LOG2STDEV1,
                              inst->ctx.lnp_seed,NONCE_Y1_BASE + inst->rejection_rounds);
      sample_gaussian_polxvec(inst->y2,inst->z2coeffs,inst->ctx.m2,LNP22_LOG2STDEV2,
                              inst->ctx.lnp_seed,NONCE_Y2_BASE + inst->rejection_rounds);
      compute_w_rows(inst->w,&inst->ctx,inst->y1,inst->y2);
      folded_linear_eval(inst->vlin,&inst->ctx,inst->y1,inst->y2);
      quad_mask_terms(inst->g1,inst->g0,&inst->ctx,inst->ctx.s1,inst->y1);
      polx_add(inst->tq,inst->a3s2,inst->g1);
      inner_product_a3(inst->vq,inst->a3,inst->y2,inst->ctx.m2);
      polx_add(inst->vq,inst->vq,inst->g0);
      fs_challenge(inst->c,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
      polxvec_polx_mul(inst->cs1,inst->c,inst->ctx.s1,inst->ctx.m1);
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_polx_mul(inst->cs2,inst->c,inst->s2,inst->ctx.m2);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);

      polxvec_centered_to_i64(inst->cs1coeffs,inst->cs1,inst->ctx.m1);
      polxvec_centered_to_i64(inst->cs2coeffs,inst->cs2,inst->ctx.m2);
      for(i=0;i<inst->ctx.m1*N;i++)
        inst->z1coeffs[i] += inst->cs1coeffs[i];
      for(i=0;i<inst->ctx.m2*N;i++)
        inst->z2coeffs[i] += inst->cs2coeffs[i];

      rej1 = reject_standard_polxvec(inst->z1coeffs,inst->cs1coeffs,inst->ctx.m1*N,
                                     sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                     LNP22_GAMMA1,
                                     inst->ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + inst->rejection_rounds);
      rej2 = reject_z2_polxvec(inst->z2coeffs,inst->cs2coeffs,inst->ctx.m2*N,
                                 sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                 LNP22_GAMMA2,
                                 inst->ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + inst->rejection_rounds);
      if(!rej1 && !rej2)
        break;
    }
  } else {
    size_t batch_base = 0;
    uint8_t *accept = _malloc(batch_rounds);
    size_t found_bi = 0;
    int found = 0;

    for(;;) {
      size_t bi;

#ifdef _OPENMP
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
        polx_add(ws->tq,inst->a3s2,ws->g1);
        inner_product_a3(ws->vq,inst->a3,ws->y2,inst->ctx.m2);
        polx_add(ws->vq,ws->vq,ws->g0);
        fs_challenge(ws->c,&inst->ctx,inst->tA,ws->w,ws->vlin,ws->tq,ws->vq);
        polxvec_polx_mul(ws->cs1,ws->c,inst->ctx.s1,inst->ctx.m1);
        polxvec_polx_mul(ws->cs2,ws->c,inst->s2,inst->ctx.m2);

        polxvec_centered_to_i64(ws->cs1coeffs,ws->cs1,inst->ctx.m1);
        polxvec_centered_to_i64(ws->cs2coeffs,ws->cs2,inst->ctx.m2);
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
      if(found)
        break;
      batch_base += batch_rounds;
    }

    {
      const lnp22_round_ws *winner = &batch_ws[found_bi];
      size_t i;

      polxvec_copy(inst->y1,winner->y1,inst->ctx.m1);
      polxvec_copy(inst->y2,winner->y2,inst->ctx.m2);
      polxvec_copy(inst->w,winner->w,inst->ctx.kmsis);
      polxvec_copy(inst->cs1,winner->cs1,inst->ctx.m1);
      polxvec_copy(inst->cs2,winner->cs2,inst->ctx.m2);
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
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);
    }

    free(accept);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_prove_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.inner_rejection_rounds = inst->rejection_rounds + 1;
  g_last_request_profile.inner_batch_rounds = batch_rounds;

#if LEOPARD_REQUEST_SELF_CHECK
  {
    polx c_check[1], lhs[1], rhs[1], tmp[1];
    size_t i;

    clock_gettime(CLOCK_MONOTONIC,&t0);
    fs_challenge(c_check,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
    if(memcmp(inst->c,c_check,sizeof(polx)) != 0)
      ok = 0;

    for(i=0;i<inst->ctx.kmsis && ok;i++) {
      polx row_eval[1];
      polxvec_sprod(row_eval,&inst->ctx.A1open[i*inst->ctx.m1],inst->z1,inst->ctx.m1);
      polxvec_sprod(tmp,&inst->ctx.A2open[i*inst->ctx.m2prime],inst->z2,inst->ctx.m2prime);
      polx_add(row_eval,row_eval,tmp);
      polx_add(row_eval,row_eval,&inst->z2[inst->ctx.m2prime + i]);
      polx_mul(rhs,inst->c,&inst->tA[i]);
      polx_sub(row_eval,row_eval,rhs);
      polx_sub(row_eval,row_eval,&inst->w[i]);
      if(!polx_iszero(row_eval))
        ok = 0;
    }

    if(ok) {
      folded_linear_eval(lhs,&inst->ctx,inst->z1,inst->z2);
      polx_mul(rhs,inst->c,inst->ufold_public);
      polx_sub(lhs,lhs,rhs);
      polx_sub(lhs,lhs,inst->vlin);
      if(!polx_iszero(lhs))
        ok = 0;
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
      if(LNP22_USE_RBIT_BINARY) {
        rhat_bits_verifier_value(bits_term,&inst->ctx,
                                 &inst->z1[inst->ctx.off[LNP22_RHAT_BITS_BLOCK]],
                                 inst->c);
        polx_add(lhs,lhs,bits_term);
      }
      polx_mul(c2,inst->c,inst->c);
      polx_mul(c2u,c2,inst->uquad_public);
      polx_sub(lhs,lhs,c2u);
      inner_product_a3(a3z2,inst->a3,inst->z2,inst->ctx.m2);
      polx_mul(mask_term,inst->c,inst->tq);
      polx_sub(mask_term,mask_term,a3z2);
      polx_sub(lhs,lhs,mask_term);
      polx_sub(lhs,lhs,inst->vq);
      if(!polx_iszero(lhs))
        ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    g_last_request_profile.inner_selfcheck_sec = leopard_elapsed_sec(&t0,&t1);
  }
#else
  g_last_request_profile.inner_selfcheck_sec = 0.0;
#endif

  if(batch_ws != NULL) {
    size_t bi;

    for(bi=0;bi<batch_rounds;bi++)
      round_ws_clear(&batch_ws[bi]);
    free(batch_ws);
  }

  return ok ? 0 : 1;
}

static int leopard_generate_request_instance_preprocessed(
    outer_instance *inst,
    const leopard_pp *pp,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS]) {
  const size_t batch_rounds = LNP22_ENABLE_REJECTION ? LNP22_BATCH_ROUNDS : 1;
  size_t batch_workers = 1;
  lnp22_round_ws *batch_ws = NULL;
  int ok = 1;
  struct timespec t0, t1;

  if(pre == NULL || pre->slot_count == 0 || slot_j >= pre->slot_count)
    return 1;

  memset(inst,0,sizeof(*inst));
  selector_ctx_init(&inst->ctx,pp->nizk1_seed);
  g_last_request_profile.a1open_cache = inst->ctx.cache_a1open;
  if(pre->slot_count > inst->ctx.preproc_count)
    return 1;
  inst->ctx.preproc_count = pre->slot_count;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  leopard_prepare_selector_lnp22_preprocessed_explicit(&inst->ctx,pre,slot_j,
                                                       tag_hash,xbits);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_prepare_sec = leopard_elapsed_sec(&t0,&t1);
  if(first_linear_relation_failure(&inst->ctx,inst->ctx.s1)) {
    fprintf(stderr,"ERROR: preprocessed request witness does not satisfy intended relations\n");
    return 1;
  }

  inst->s2 = alloc_polxvec(inst->ctx.m2);
  inst->tA = alloc_polxvec(inst->ctx.kmsis);
  inst->a3 = alloc_polxvec(inst->ctx.m2);
  inst->w = alloc_polxvec(inst->ctx.kmsis);
  inst->z1 = alloc_polxvec(inst->ctx.m1);
  inst->z2 = alloc_polxvec(inst->ctx.m2);
  inst->cs1 = alloc_polxvec(inst->ctx.m1);
  inst->cs2 = alloc_polxvec(inst->ctx.m2);
  inst->y1 = alloc_polxvec(inst->ctx.m1);
  inst->y2 = alloc_polxvec(inst->ctx.m2);
  inst->z1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->z1coeffs));
  inst->z2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->z2coeffs));
  inst->cs1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*inst->cs1coeffs));
  inst->cs2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*inst->cs2coeffs));
  if(batch_rounds > 1) {
    size_t bi;

#ifdef _OPENMP
    batch_workers = (size_t)omp_get_max_threads();
#endif
    if(batch_workers > batch_rounds)
      batch_workers = batch_rounds;
    batch_ws = _malloc(batch_rounds * sizeof(*batch_ws));
    for(bi=0;bi<batch_rounds;bi++)
      round_ws_init(&batch_ws[bi],&inst->ctx);
  }

  compute_public_quad_rhs(inst->uquad_public,&inst->ctx);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  sample_ternary_polxvec(inst->s2,inst->ctx.m2,inst->ctx.lnp_seed,NONCE_S2_BASE);
  compute_a3(inst->a3,&inst->ctx);
  inner_product_a3(inst->a3s2,inst->a3,inst->s2,inst->ctx.m2);
  {
    struct timespec lt0, lt1;
    clock_gettime(CLOCK_MONOTONIC,&lt0);
    if(prepare_linear_fold_for_proof(&inst->ctx,inst->ufold_public,inst->s2))
      return 1;
    clock_gettime(CLOCK_MONOTONIC,&lt1);
    g_last_request_profile.nizk_prepare_sec += leopard_elapsed_sec(&lt0,&lt1);
  }
  compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);

  if(batch_rounds <= 1 || batch_ws == NULL) {
    for(inst->rejection_rounds = 0;;inst->rejection_rounds++) {
      size_t i;
      int rej1, rej2;

      sample_gaussian_polxvec(inst->y1,inst->z1coeffs,inst->ctx.m1,LNP22_LOG2STDEV1,
                              inst->ctx.lnp_seed,NONCE_Y1_BASE + inst->rejection_rounds);
      sample_gaussian_polxvec(inst->y2,inst->z2coeffs,inst->ctx.m2,LNP22_LOG2STDEV2,
                              inst->ctx.lnp_seed,NONCE_Y2_BASE + inst->rejection_rounds);
      compute_w_rows(inst->w,&inst->ctx,inst->y1,inst->y2);
      folded_linear_eval(inst->vlin,&inst->ctx,inst->y1,inst->y2);
      quad_mask_terms(inst->g1,inst->g0,&inst->ctx,inst->ctx.s1,inst->y1);
      polx_add(inst->tq,inst->a3s2,inst->g1);
      inner_product_a3(inst->vq,inst->a3,inst->y2,inst->ctx.m2);
      polx_add(inst->vq,inst->vq,inst->g0);
      fs_challenge(inst->c,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
      polxvec_polx_mul(inst->cs1,inst->c,inst->ctx.s1,inst->ctx.m1);
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_polx_mul(inst->cs2,inst->c,inst->s2,inst->ctx.m2);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);

      polxvec_centered_to_i64(inst->cs1coeffs,inst->cs1,inst->ctx.m1);
      polxvec_centered_to_i64(inst->cs2coeffs,inst->cs2,inst->ctx.m2);
      for(i=0;i<inst->ctx.m1*N;i++)
        inst->z1coeffs[i] += inst->cs1coeffs[i];
      for(i=0;i<inst->ctx.m2*N;i++)
        inst->z2coeffs[i] += inst->cs2coeffs[i];

      rej1 = reject_standard_polxvec(inst->z1coeffs,inst->cs1coeffs,inst->ctx.m1*N,
                                     sigma_from_log2stdev(LNP22_LOG2STDEV1),
                                     LNP22_GAMMA1,
                                     inst->ctx.lnp_seed,NONCE_Y1_BASE + 1000000 + inst->rejection_rounds);
      rej2 = reject_z2_polxvec(inst->z2coeffs,inst->cs2coeffs,inst->ctx.m2*N,
                                 sigma_from_log2stdev(LNP22_LOG2STDEV2),
                                 LNP22_GAMMA2,
                                 inst->ctx.lnp_seed,NONCE_Y2_BASE + 1000000 + inst->rejection_rounds);
      if(!rej1 && !rej2)
        break;
    }
  } else {
    size_t batch_base = 0;
    uint8_t *accept = _malloc(batch_rounds);
    size_t found_bi = 0;
    int found = 0;

    for(;;) {
      size_t bi;

#ifdef _OPENMP
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
        polx_add(ws->tq,inst->a3s2,ws->g1);
        inner_product_a3(ws->vq,inst->a3,ws->y2,inst->ctx.m2);
        polx_add(ws->vq,ws->vq,ws->g0);
        fs_challenge(ws->c,&inst->ctx,inst->tA,ws->w,ws->vlin,ws->tq,ws->vq);
        polxvec_polx_mul(ws->cs1,ws->c,inst->ctx.s1,inst->ctx.m1);
        polxvec_polx_mul(ws->cs2,ws->c,inst->s2,inst->ctx.m2);

        polxvec_centered_to_i64(ws->cs1coeffs,ws->cs1,inst->ctx.m1);
        polxvec_centered_to_i64(ws->cs2coeffs,ws->cs2,inst->ctx.m2);
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
      if(found)
        break;
      batch_base += batch_rounds;
    }

    {
      const lnp22_round_ws *winner = &batch_ws[found_bi];
      size_t i;

      polxvec_copy(inst->y1,winner->y1,inst->ctx.m1);
      polxvec_copy(inst->y2,winner->y2,inst->ctx.m2);
      polxvec_copy(inst->w,winner->w,inst->ctx.kmsis);
      polxvec_copy(inst->cs1,winner->cs1,inst->ctx.m1);
      polxvec_copy(inst->cs2,winner->cs2,inst->ctx.m2);
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
      polxvec_add(inst->z1,inst->cs1,inst->y1,inst->ctx.m1);
      polxvec_add(inst->z2,inst->cs2,inst->y2,inst->ctx.m2);
    }

    free(accept);
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_prove_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.inner_rejection_rounds = inst->rejection_rounds + 1;
  g_last_request_profile.inner_batch_rounds = batch_rounds;

#if LEOPARD_REQUEST_SELF_CHECK
  {
    polx c_check[1], lhs[1], rhs[1], tmp[1];
    size_t i;

    clock_gettime(CLOCK_MONOTONIC,&t0);
    fs_challenge(c_check,&inst->ctx,inst->tA,inst->w,inst->vlin,inst->tq,inst->vq);
    if(memcmp(inst->c,c_check,sizeof(polx)) != 0)
      ok = 0;

    for(i=0;i<inst->ctx.kmsis && ok;i++) {
      polx row_eval[1];
      polxvec_sprod(row_eval,&inst->ctx.A1open[i*inst->ctx.m1],inst->z1,inst->ctx.m1);
      polxvec_sprod(tmp,&inst->ctx.A2open[i*inst->ctx.m2prime],inst->z2,inst->ctx.m2prime);
      polx_add(row_eval,row_eval,tmp);
      polx_add(row_eval,row_eval,&inst->z2[inst->ctx.m2prime + i]);
      polx_mul(rhs,inst->c,&inst->tA[i]);
      polx_sub(row_eval,row_eval,rhs);
      polx_sub(row_eval,row_eval,&inst->w[i]);
      if(!polx_iszero(row_eval))
        ok = 0;
    }

    if(ok) {
      folded_linear_eval(lhs,&inst->ctx,inst->z1,inst->z2);
      polx_mul(rhs,inst->c,inst->ufold_public);
      polx_sub(lhs,lhs,rhs);
      polx_sub(lhs,lhs,inst->vlin);
      if(!polx_iszero(lhs))
        ok = 0;
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
      if(LNP22_USE_RBIT_BINARY) {
        rhat_bits_verifier_value(bits_term,&inst->ctx,
                                 &inst->z1[inst->ctx.off[LNP22_RHAT_BITS_BLOCK]],
                                 inst->c);
        polx_add(lhs,lhs,bits_term);
      }
      polx_mul(c2,inst->c,inst->c);
      polx_mul(c2u,c2,inst->uquad_public);
      polx_sub(lhs,lhs,c2u);
      inner_product_a3(a3z2,inst->a3,inst->z2,inst->ctx.m2);
      polx_mul(mask_term,inst->c,inst->tq);
      polx_sub(mask_term,mask_term,a3z2);
      polx_sub(lhs,lhs,mask_term);
      polx_sub(lhs,lhs,inst->vq);
      if(!polx_iszero(lhs))
        ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    g_last_request_profile.inner_selfcheck_sec = leopard_elapsed_sec(&t0,&t1);
  }
#else
  g_last_request_profile.inner_selfcheck_sec = 0.0;
#endif

  if(batch_ws != NULL) {
    size_t bi;

    for(bi=0;bi<batch_rounds;bi++)
      round_ws_clear(&batch_ws[bi]);
    free(batch_ws);
  }

  return ok ? 0 : 1;
}
#else
#ifndef LEOPARD_DIRECT_COMBINE_TA_W
#define LEOPARD_DIRECT_COMBINE_TA_W (LOGQ == 72)
#endif
static int leopard_finish_direct_inner_instance(outer_instance *inst) {
  int64_t *s1coeffs = NULL;
  int64_t *s2coeffs = NULL;
  int64_t ccoeffs[N];
  struct timespec prove0, prove1;

  clock_gettime(CLOCK_MONOTONIC,&prove0);
  inst->s2 = alloc_polxvec(inst->ctx.m2);
  inst->tA = alloc_polxvec(inst->ctx.kmsis);
  inst->a3 = alloc_polxvec(inst->ctx.m2);
  inst->w = alloc_polxvec(inst->ctx.kmsis);
  inst->z1 = alloc_polxvec(inst->ctx.m1);
  inst->z2 = alloc_polxvec(inst->ctx.m2);
  inst->y1 = alloc_polxvec(inst->ctx.m1);
  inst->y2 = alloc_polxvec(inst->ctx.m2);
  inst->z1coeffs =
      _malloc(inst->ctx.m1 * N * sizeof(*inst->z1coeffs));
  inst->z2coeffs =
      _malloc(inst->ctx.m2 * N * sizeof(*inst->z2coeffs));
  inst->cs1coeffs =
      _malloc(inst->ctx.m1 * N * sizeof(*inst->cs1coeffs));
  inst->cs2coeffs =
      _malloc(inst->ctx.m2 * N * sizeof(*inst->cs2coeffs));
  s1coeffs = _malloc(inst->ctx.m1 * N * sizeof(*s1coeffs));
  s2coeffs = _malloc(inst->ctx.m2 * N * sizeof(*s2coeffs));
  witness_coeffs_to_i64(s1coeffs,&inst->ctx);

  sample_ternary_polxvec(inst->s2,inst->ctx.m2,
                         inst->ctx.lnp_seed,NONCE_S2_BASE);
  polxvec_centered_to_i64(s2coeffs,inst->s2,inst->ctx.m2);
  compute_a3(inst->a3,&inst->ctx);
  inner_product_a3(inst->a3s2,inst->a3,inst->s2,inst->ctx.m2);
#if LEOPARD_DIRECT_COMBINE_TA_W
  /* The first masks do not depend on tA. Reuse each generated A1 row for
   * both products, retaining the same mask nonces and folding order. */
  sample_gaussian_polxvec(inst->y1,inst->z1coeffs,inst->ctx.m1,
                          LNP22_LOG2STDEV1,inst->ctx.lnp_seed,NONCE_Y1_BASE);
  sample_gaussian_polxvec(inst->y2,inst->z2coeffs,inst->ctx.m2,
                          LNP22_LOG2STDEV2,inst->ctx.lnp_seed,NONCE_Y2_BASE);
  compute_tA_w_rows(inst->tA,inst->w,&inst->ctx,inst->ctx.s1,inst->s2,
                    inst->y1,inst->y2);
#else
  compute_tA_rows(inst->tA,&inst->ctx,inst->ctx.s1,inst->s2);
#endif
  bind_prechallenge_commitment(&inst->ctx,inst->tA);
  if(prepare_linear_fold_for_proof(
         &inst->ctx,inst->ufold_public,inst->s2)) {
    free(s1coeffs);
    free(s2coeffs);
    return 1;
  }
  compute_public_quad_rhs(inst->uquad_public,&inst->ctx);
  derive_direct_chain_fold(&inst->ctx,inst->tA);
  compute_quad_g1_coeffs(&inst->ctx);
  {
    polx direct_check[1];

    quad_value(direct_check,&inst->ctx,inst->ctx.s1);
    if(!polx_iszero(direct_check)) {
      fprintf(stderr,
              "ERROR: direct-chain base quadratic relation failed\n");
      free(s1coeffs);
      free(s2coeffs);
      return 1;
    }
  }

  for(inst->rejection_rounds=0;;inst->rejection_rounds++) {
    size_t i;
    int rej1, rej2;

#if LEOPARD_DIRECT_COMBINE_TA_W
    if(inst->rejection_rounds != 0) {
#endif
    sample_gaussian_polxvec(
        inst->y1,inst->z1coeffs,inst->ctx.m1,
        LNP22_LOG2STDEV1,inst->ctx.lnp_seed,
        NONCE_Y1_BASE + inst->rejection_rounds);
    sample_gaussian_polxvec(
        inst->y2,inst->z2coeffs,inst->ctx.m2,
        LNP22_LOG2STDEV2,inst->ctx.lnp_seed,
        NONCE_Y2_BASE + inst->rejection_rounds);
    compute_w_rows(inst->w,&inst->ctx,inst->y1,inst->y2);
#if LEOPARD_DIRECT_COMBINE_TA_W
    }
#endif
    folded_linear_eval(inst->vlin,&inst->ctx,inst->y1,inst->y2);
    quad_mask_terms(inst->g1,inst->g0,
                    &inst->ctx,inst->ctx.s1,inst->y1);
    polx_add(inst->tq,inst->a3s2,inst->g1);
    inner_product_a3(inst->vq,inst->a3,inst->y2,inst->ctx.m2);
    polx_add(inst->vq,inst->vq,inst->g0);
    fs_challenge(inst->c,&inst->ctx,inst->tA,inst->w,
                 inst->vlin,inst->tq,inst->vq);
    polxvec_centered_to_i64(ccoeffs,inst->c,1);
    challenge_mul_coeffs(
        inst->cs1coeffs,ccoeffs,s1coeffs,inst->ctx.m1);
    challenge_mul_coeffs(
        inst->cs2coeffs,ccoeffs,s2coeffs,inst->ctx.m2);
    for(i=0;i<inst->ctx.m1*N;i++)
      inst->z1coeffs[i] += inst->cs1coeffs[i];
    for(i=0;i<inst->ctx.m2*N;i++)
      inst->z2coeffs[i] += inst->cs2coeffs[i];

    if(LNP22_ENABLE_REJECTION) {
      rej1 = reject_standard_polxvec(
          inst->z1coeffs,inst->cs1coeffs,inst->ctx.m1*N,
          sigma_from_log2stdev(LNP22_LOG2STDEV1),LNP22_GAMMA1,
          inst->ctx.lnp_seed,
          NONCE_Y1_BASE + 1000000 + inst->rejection_rounds);
      rej2 = reject_z2_polxvec(
          inst->z2coeffs,inst->cs2coeffs,inst->ctx.m2*N,
          sigma_from_log2stdev(LNP22_LOG2STDEV2),LNP22_GAMMA2,
          inst->ctx.lnp_seed,
          NONCE_Y2_BASE + 1000000 + inst->rejection_rounds);
    } else {
      rej1 = 0;
      rej2 = 0;
    }
    if(!rej1 && !rej2) {
      polxvec_fromint64vec(
          inst->z1,inst->ctx.m1,1,inst->z1coeffs);
      polxvec_fromint64vec(
          inst->z2,inst->ctx.m2,1,inst->z2coeffs);
      break;
    }
  }

  free(s1coeffs);
  free(s2coeffs);
  clock_gettime(CLOCK_MONOTONIC,&prove1);
  g_last_request_profile.inner_prove_sec =
      leopard_elapsed_sec(&prove0,&prove1);
  g_last_request_profile.inner_rejection_rounds =
      inst->rejection_rounds + 1;
  g_last_request_profile.inner_batch_rounds = 1;
  g_last_request_profile.inner_selfcheck_sec = 0.0;
  return 0;
}

static int leopard_generate_request_instance(
    outer_instance *inst,
    const leopard_pp *pp,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS]) {
  struct timespec prep0, prep1;

  memset(inst,0,sizeof(*inst));
  selector_ctx_init(&inst->ctx,pp->nizk1_seed);
  g_last_request_profile.a1open_cache = inst->ctx.cache_a1open;
  clock_gettime(CLOCK_MONOTONIC,&prep0);
  leopard_prepare_selector_direct(
      &inst->ctx,NULL,0,tag_hash,xbits);
  clock_gettime(CLOCK_MONOTONIC,&prep1);
  g_last_request_profile.inner_prepare_sec =
      leopard_elapsed_sec(&prep0,&prep1);
  return leopard_finish_direct_inner_instance(inst);
}

static int leopard_generate_request_instance_preprocessed(
    outer_instance *inst,
    const leopard_pp *pp,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t tag_hash[LEOPARD_TAG_BYTES],
    const uint8_t xbits[LEOPARD_X_BITS]) {
  struct timespec prep0, prep1;

  memset(inst,0,sizeof(*inst));
  selector_ctx_init(&inst->ctx,pp->nizk1_seed);
  g_last_request_profile.a1open_cache = inst->ctx.cache_a1open;
  if(pre->slot_count > inst->ctx.preproc_count)
    return 1;
  inst->ctx.preproc_count = pre->slot_count;
  clock_gettime(CLOCK_MONOTONIC,&prep0);
  leopard_prepare_selector_direct(
      &inst->ctx,pre,slot_j,tag_hash,xbits);
  clock_gettime(CLOCK_MONOTONIC,&prep1);
  g_last_request_profile.inner_prepare_sec =
      leopard_elapsed_sec(&prep0,&prep1);
  return leopard_finish_direct_inner_instance(inst);
}
#endif

#if !LNP22_DIRECT_CHAIN
static size_t leopard_outer_tail_len(const selector_lnp22_ctx *ctx) {
  return ctx->rhat_len + ctx->umc_eua_len + ctx->umc_eub_len + ctx->umc_r_len +
         ctx->umc_va_len + ctx->umc_vb_len;
}

static void leopard_fill_two_digit_views(polx *dst0, polx *dst1,
                                         const polx *src, size_t len,
                                         int64_t base) {
  size_t i;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(len > 32)
#endif
  for(i=0;i<len;i++)
    memcpy(&dst0[i],&src[i],sizeof(*dst0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(len > 32)
#endif
  for(i=0;i<len;i++)
    polx_scale(&dst1[i],&src[i],base);
}

static void leopard_fill_three_digit_views(polx *dst0, polx *dst1, polx *dst2,
                                           const polx *src, size_t len,
                                           int64_t base, int64_t base2) {
  size_t i;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(len > 32)
#endif
  for(i=0;i<len;i++) {
    memcpy(&dst0[i],&src[i],sizeof(*dst0));
    polx_scale(&dst1[i],&src[i],base);
    polx_scale(&dst2[i],&src[i],base2);
  }
}

static void leopard_outer_witness_shape(size_t n[OUT_R],
                                        const selector_lnp22_ctx *ctx) {
  const size_t tail_len = leopard_outer_tail_len(ctx);

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
  n[OUT_RBIT0] = ctx->rhat_bits_len; n[OUT_RBIT1] = ctx->rhat_bits_len;
  n[OUT_Z20] = ctx->m2;    n[OUT_Z21] = ctx->m2;
}

static int leopard_prepare_outer_witness(witness *wt, const outer_instance *inst,
                                         uint64_t *betasq_out) {
  size_t n[OUT_R];
  size_t i;
  size_t maxcount;
  int64_t *b_d[3] = {NULL,NULL,NULL};
  int64_t *x_d[2] = {NULL,NULL};
  int64_t *z_d[3] = {NULL,NULL,NULL};
  int64_t *y0_d[2] = {NULL,NULL};
  int64_t *y1_d[2] = {NULL,NULL};
  int64_t *t_d[2] = {NULL,NULL};
  int64_t *rbit_d[2] = {NULL,NULL};
  int64_t *z2_d[2] = {NULL,NULL};
  uint64_t betasq = 0;
  const selector_lnp22_ctx *ctx = &inst->ctx;

  leopard_outer_witness_shape(n,ctx);
  init_witness_raw(wt,OUT_R,n);
  maxcount = MAX(ctx->b_len,
             MAX(ctx->x_len,
             MAX(ctx->z_len,
             MAX(ctx->y_len,
             MAX(MAX(leopard_outer_tail_len(ctx),ctx->rhat_bits_len),ctx->m2)))))*N;
  for(i=0;i<3;i++) {
    b_d[i] = _malloc(maxcount*sizeof(*b_d[i]));
    z_d[i] = _malloc(maxcount*sizeof(*z_d[i]));
  }
  for(i=0;i<2;i++) {
    x_d[i] = _malloc(maxcount*sizeof(*x_d[i]));
    y0_d[i] = _malloc(maxcount*sizeof(*y0_d[i]));
    y1_d[i] = _malloc(maxcount*sizeof(*y1_d[i]));
    t_d[i] = _malloc(maxcount*sizeof(*t_d[i]));
    rbit_d[i] = _malloc(maxcount*sizeof(*rbit_d[i]));
    z2_d[i] = _malloc(maxcount*sizeof(*z2_d[i]));
  }
  if(!b_d[0] || !b_d[1] ||
#if OUTER_BZ_DIGITS == 3
     !b_d[2] ||
#endif
     !x_d[0] || !x_d[1] ||
     !z_d[0] || !z_d[1] ||
#if OUTER_BZ_DIGITS == 3
     !z_d[2] ||
#endif
     !y0_d[0] || !y0_d[1] ||
     !y1_d[0] || !y1_d[1] ||
     !t_d[0] || !t_d[1] ||
     !rbit_d[0] || !rbit_d[1] ||
     !z2_d[0] || !z2_d[1]) {
    fprintf(stderr,"ERROR: outer witness allocation failure\n");
    goto fail;
  }

#if OUTER_BZ_DIGITS == 2
  decompose_centered_two_digits(b_d[0],b_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_B_BLOCK]*N],
                                ctx->b_len*N,OUTER_BZ_DECOMP_BASE);
#else
  decompose_centered_three_digits(b_d[0],b_d[1],b_d[2],
                                  &inst->z1coeffs[ctx->off[LNP22_B_BLOCK]*N],
                                  ctx->b_len*N,OUTER_BZ_DECOMP_BASE);
#endif
  decompose_centered_two_digits(x_d[0],x_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_X_BLOCK]*N],
                                ctx->x_len*N,OUTER_XY_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 2
  decompose_centered_two_digits(z_d[0],z_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_Z_BLOCK]*N],
                                ctx->z_len*N,OUTER_BZ_DECOMP_BASE);
#else
  decompose_centered_three_digits(z_d[0],z_d[1],z_d[2],
                                  &inst->z1coeffs[ctx->off[LNP22_Z_BLOCK]*N],
                                  ctx->z_len*N,OUTER_BZ_DECOMP_BASE);
#endif
  decompose_centered_two_digits(y0_d[0],y0_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_Y0_BLOCK]*N],
                                ctx->y_len*N,OUTER_XY_DECOMP_BASE);
  decompose_centered_two_digits(y1_d[0],y1_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_Y1_BLOCK]*N],
                                ctx->y_len*N,OUTER_XY_DECOMP_BASE);
  decompose_centered_two_digits(t_d[0],t_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_RHAT_BLOCK]*N],
                                leopard_outer_tail_len(ctx)*N,OUTER_XY_DECOMP_BASE);
  decompose_centered_two_digits(rbit_d[0],rbit_d[1],
                                &inst->z1coeffs[ctx->off[LNP22_RHAT_BITS_BLOCK]*N],
                                ctx->rhat_bits_len*N,OUTER_XY_DECOMP_BASE);
  decompose_centered_two_digits(z2_d[0],z2_d[1],inst->z2coeffs,ctx->m2*N,OUTER_Z2_DECOMP_BASE);

  if(set_witness_vector_raw(wt,OUT_B0,ctx->b_len,1,b_d[0]) ||
     set_witness_vector_raw(wt,OUT_B1,ctx->b_len,1,b_d[1]) ||
#if OUTER_BZ_DIGITS == 3
     set_witness_vector_raw(wt,OUT_B2,ctx->b_len,1,b_d[2]) ||
#endif
     set_witness_vector_raw(wt,OUT_X0,ctx->x_len,1,x_d[0]) ||
     set_witness_vector_raw(wt,OUT_X1,ctx->x_len,1,x_d[1]) ||
     set_witness_vector_raw(wt,OUT_Z0,ctx->z_len,1,z_d[0]) ||
     set_witness_vector_raw(wt,OUT_Z1,ctx->z_len,1,z_d[1]) ||
#if OUTER_BZ_DIGITS == 3
     set_witness_vector_raw(wt,OUT_Z2,ctx->z_len,1,z_d[2]) ||
#endif
     set_witness_vector_raw(wt,OUT_Y00,ctx->y_len,1,y0_d[0]) ||
     set_witness_vector_raw(wt,OUT_Y01,ctx->y_len,1,y0_d[1]) ||
     set_witness_vector_raw(wt,OUT_Y10,ctx->y_len,1,y1_d[0]) ||
     set_witness_vector_raw(wt,OUT_Y11,ctx->y_len,1,y1_d[1]) ||
     set_witness_vector_raw(wt,OUT_T0,leopard_outer_tail_len(ctx),1,t_d[0]) ||
     set_witness_vector_raw(wt,OUT_T1,leopard_outer_tail_len(ctx),1,t_d[1]) ||
     set_witness_vector_raw(wt,OUT_RBIT0,ctx->rhat_bits_len,1,rbit_d[0]) ||
     set_witness_vector_raw(wt,OUT_RBIT1,ctx->rhat_bits_len,1,rbit_d[1]) ||
     set_witness_vector_raw(wt,OUT_Z20,ctx->m2,1,z2_d[0]) ||
     set_witness_vector_raw(wt,OUT_Z21,ctx->m2,1,z2_d[1])) {
    fprintf(stderr,"ERROR: failed to set decomposed outer witness vectors\n");
    goto fail;
  }

  for(i=0;i<wt->r;i++)
    betasq += wt->normsq[i];
  betasq += 16*N;
  *betasq_out = betasq;

  for(i=0;i<3;i++) {
    free(b_d[i]);
    free(z_d[i]);
  }
  for(i=0;i<2;i++) {
    free(x_d[i]);
    free(y0_d[i]);
    free(y1_d[i]);
    free(t_d[i]);
    free(rbit_d[i]);
    free(z2_d[i]);
  }
  return 0;

fail:
  for(i=0;i<3;i++) {
    free(b_d[i]);
    free(z_d[i]);
  }
  for(i=0;i<2;i++) {
    free(x_d[i]);
    free(y0_d[i]);
    free(y1_d[i]);
    free(t_d[i]);
    free(rbit_d[i]);
    free(z2_d[i]);
  }
  free_witness(wt);
  memset(wt,0,sizeof(*wt));
  return 1;
}

static int leopard_prepare_outer_statement(prncplstmnt *st,
                                           const selector_lnp22_ctx *ctx,
                                           const polx *tA, const polx *w,
                                           const polx *c, const polx *vlin,
                                           const polx *tq, const polx *vq,
                                           uint64_t betasq,
                                           int prover_ctx,
                                           leopard_outer_statement_profile *prof) {
  size_t n[OUT_R];
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;
  polx *openfold = NULL, *open1 = NULL, *open2 = NULL, *openrhs = NULL;
  polx *quadz = NULL, *quadx = NULL, *quadbits = NULL, *quadz2 = NULL;
  polx *a3 = NULL;
  polx ufold_public[1], linrhs[1], quadrhs[1];
  (void)prover_ctx;
#if OUTER_BZ_DIGITS == 3
  const int64_t bzbase2 = OUTER_BZ_DECOMP_BASE * OUTER_BZ_DECOMP_BASE;
#endif
  const size_t tail_len = leopard_outer_tail_len(ctx);
  const size_t rbit_len = ctx->rhat_bits_len;
  size_t i;
  struct timespec t0, t1;

  if(prof != NULL)
    memset(prof,0,sizeof(*prof));
  leopard_outer_witness_shape(n,ctx);
  if(init_prncplstmnt_raw(st,OUT_R,n,betasq,3,1) != 0)
    return 1;
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
  a3 = alloc_polxvec(ctx->m2);
  if(!openfold || !open1 || !open2 || !openrhs ||
     !quadz || !quadx || !quadbits || !quadz2 || !a3)
    goto fail;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  derive_outer_openfold(openfold,ctx,tA,w,vlin,tq,vq,c);
  compute_openfold_coeffs(open1,open2,openrhs,ctx,tA,w,c,openfold);
  compute_public_fold_rhs(ufold_public,(selector_lnp22_ctx *)ctx);
  polx_mul(linrhs,c,ufold_public);
  polx_add(linrhs,linrhs,vlin);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->fold_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  compute_a3(a3,(selector_lnp22_ctx *)ctx);
  compute_quad_linear_coeffs(quadz,quadx,quadbits,quadz2,ctx,a3,c);
  polx_mul(quadrhs,c,tq);
  polx_add(quadrhs,quadrhs,vq);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->quad_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  cnst = &st->cnst[0];
  buf = init_sparsecnst_half(cnst,OUT_R,OUT_R,
                             OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
                             OUTER_BZ_DIGITS*ctx->z_len +
                             2*ctx->y_len + 2*ctx->y_len +
                             2*tail_len + 2*rbit_len + 2*ctx->m2,
                             1,0,0);
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
#if OUTER_BZ_DIGITS == 2
  leopard_fill_two_digit_views(cnst->phi[OUT_B0],cnst->phi[OUT_B1],
                               &open1[ctx->off[LNP22_B_BLOCK]],ctx->b_len,
                               OUTER_BZ_DECOMP_BASE);
#else
  leopard_fill_three_digit_views(cnst->phi[0],cnst->phi[1],cnst->phi[2],
                                 &open1[ctx->off[LNP22_B_BLOCK]],ctx->b_len,
                                 OUTER_BZ_DECOMP_BASE,bzbase2);
#endif
  leopard_fill_two_digit_views(cnst->phi[OUT_X0],cnst->phi[OUT_X1],
                               &open1[ctx->off[LNP22_X_BLOCK]],ctx->x_len,
                               OUTER_XY_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 2
  leopard_fill_two_digit_views(cnst->phi[OUT_Z0],cnst->phi[OUT_Z1],
                               &open1[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,
                               OUTER_BZ_DECOMP_BASE);
#else
  leopard_fill_three_digit_views(cnst->phi[OUT_Z0],cnst->phi[OUT_Z1],cnst->phi[OUT_Z2],
                                 &open1[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,
                                 OUTER_BZ_DECOMP_BASE,bzbase2);
#endif
  leopard_fill_two_digit_views(cnst->phi[OUT_Y00],cnst->phi[OUT_Y01],
                               &open1[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_Y10],cnst->phi[OUT_Y11],
                               &open1[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_T0],cnst->phi[OUT_T1],
                               &open1[ctx->off[LNP22_RHAT_BLOCK]],tail_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_RBIT0],cnst->phi[OUT_RBIT1],
                               &open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_Z20],cnst->phi[OUT_Z21],open2,ctx->m2,
                               OUTER_Z2_DECOMP_BASE);
  for(i=0;i<OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
             OUTER_BZ_DIGITS*ctx->z_len + 2*ctx->y_len + 2*ctx->y_len +
             2*tail_len + 2*rbit_len + 2*ctx->m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
  polxvec_copy(cnst->b,openrhs,1);
  OUTER_STMT_ABSORB_POLX(&shakectx,cnst->b);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->cnst1_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  cnst = &st->cnst[1];
  buf = init_sparsecnst_half(cnst,OUT_R,OUT_R,
                             OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
                             OUTER_BZ_DIGITS*ctx->z_len +
                             2*ctx->y_len + 2*ctx->y_len +
                             2*tail_len + 2*rbit_len + 2*ctx->m2,
                             1,0,0);
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
#if OUTER_BZ_DIGITS == 2
  leopard_fill_two_digit_views(cnst->phi[OUT_B0],cnst->phi[OUT_B1],
                               &ctx->lincoeff[ctx->off[LNP22_B_BLOCK]],ctx->b_len,
                               OUTER_BZ_DECOMP_BASE);
#else
  leopard_fill_three_digit_views(cnst->phi[OUT_B0],cnst->phi[OUT_B1],cnst->phi[OUT_B2],
                                 &ctx->lincoeff[ctx->off[LNP22_B_BLOCK]],ctx->b_len,
                                 OUTER_BZ_DECOMP_BASE,bzbase2);
#endif
  leopard_fill_two_digit_views(cnst->phi[OUT_X0],cnst->phi[OUT_X1],
                               &ctx->lincoeff[ctx->off[LNP22_X_BLOCK]],ctx->x_len,
                               OUTER_XY_DECOMP_BASE);
#if OUTER_BZ_DIGITS == 2
  leopard_fill_two_digit_views(cnst->phi[OUT_Z0],cnst->phi[OUT_Z1],
                               &ctx->lincoeff[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,
                               OUTER_BZ_DECOMP_BASE);
#else
  leopard_fill_three_digit_views(cnst->phi[OUT_Z0],cnst->phi[OUT_Z1],cnst->phi[OUT_Z2],
                                 &ctx->lincoeff[ctx->off[LNP22_Z_BLOCK]],ctx->z_len,
                                 OUTER_BZ_DECOMP_BASE,bzbase2);
#endif
  leopard_fill_two_digit_views(cnst->phi[OUT_Y00],cnst->phi[OUT_Y01],
                               &ctx->lincoeff[ctx->off[LNP22_Y0_BLOCK]],ctx->y_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_Y10],cnst->phi[OUT_Y11],
                               &ctx->lincoeff[ctx->off[LNP22_Y1_BLOCK]],ctx->y_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_T0],cnst->phi[OUT_T1],
                               &ctx->lincoeff[ctx->off[LNP22_RHAT_BLOCK]],tail_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_RBIT0],cnst->phi[OUT_RBIT1],
                               &ctx->lincoeff[ctx->off[LNP22_RHAT_BITS_BLOCK]],rbit_len,
                               OUTER_XY_DECOMP_BASE);
  leopard_fill_two_digit_views(cnst->phi[OUT_Z20],cnst->phi[OUT_Z21],
                               ctx->lincoeff2,ctx->m2,OUTER_Z2_DECOMP_BASE);
  for(i=0;i<OUTER_BZ_DIGITS*ctx->b_len + 2*ctx->x_len +
             OUTER_BZ_DIGITS*ctx->z_len + 2*ctx->y_len + 2*ctx->y_len +
             2*tail_len + 2*rbit_len + 2*ctx->m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
  polxvec_copy(cnst->b,linrhs,1);
  OUTER_STMT_ABSORB_POLX(&shakectx,cnst->b);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->cnst2_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  cnst = &st->cnst[2];
  buf = init_sparsecnst_half(cnst,OUT_R,OUTER_BZ_DIGITS + 6,
                             2*ctx->x_len + OUTER_BZ_DIGITS*ctx->z_len +
                             2*rbit_len + 2*ctx->m2,1,1,0);
  {
    size_t p = 0;
    const size_t q_x0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_X0,ctx->x_len);
    const size_t q_x1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_X1,ctx->x_len);
    const size_t q_z0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z0,ctx->z_len);
    const size_t q_z1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z1,ctx->z_len);
#if OUTER_BZ_DIGITS == 3
    const size_t q_z2 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z2,ctx->z_len);
#endif
    const size_t q_rbit0 = p; buf = append_cnst_block(cnst,&p,buf,OUT_RBIT0,rbit_len);
    const size_t q_rbit1 = p; buf = append_cnst_block(cnst,&p,buf,OUT_RBIT1,rbit_len);
    const size_t q_z20 = p; buf = append_cnst_block(cnst,&p,buf,OUT_Z20,ctx->m2);
    const size_t q_z21 = p; (void)append_cnst_block(cnst,&p,buf,OUT_Z21,ctx->m2);
    for(i=0;i<p;i++)
      polxvec_setzero(cnst->phi[i],cnst->len[i]);
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
  for(i=0;i<2*ctx->x_len + OUTER_BZ_DIGITS*ctx->z_len +
             2*rbit_len + 2*ctx->m2;i++)
    OUTER_STMT_ABSORB_POLX(&shakectx,&cnst->phi[0][i]);
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
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->cnst3_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  if(prof != NULL)
    prof->hash_sec = leopard_elapsed_sec(&t0,&t1);

  free(openfold);
  free(open1);
  free(open2);
  free(openrhs);
  free(quadz);
  free(quadx);
  free(quadbits);
  free(quadz2);
  free(a3);
  return 0;

fail:
  free_prncplstmnt(st);
  free(openfold);
  free(open1);
  free(open2);
  free(openrhs);
  free(quadz);
  free(quadx);
  free(quadbits);
  free(quadz2);
  free(a3);
  return 1;
}
#else
#ifndef OUTER_PACK_LINEAR_TAIL
#define OUTER_PACK_LINEAR_TAIL 1
#endif

/* Concatenate blocks that never occur in quadratic terms. This is a
 * permutation of coefficients, so the joint norm budget is unchanged. */
static void outer_linear_tail_layout(const size_t old_n[OUT_R],
                                     size_t map[OUT_R], size_t off[OUT_R],
                                     size_t new_n[OUT_R]) {
  size_t next=0, total=0;
  for(size_t i=0;i<OUT_R;i++) {
    off[i]=0;
    if(i==OUT_T1 || i==OUT_Z20 || i==OUT_Z21) continue;
    map[i]=next;
    new_n[next++]=old_n[i];
  }
  const size_t tails[]={OUT_T0,OUT_T1,OUT_Z20,OUT_Z21};
  for(size_t j=0;j<4;j++) {
    const size_t i=tails[j];
    map[i]=map[OUT_T0];
    off[i]=total;
    total+=old_n[i];
  }
  new_n[map[OUT_T0]]=total;
}

static void outer_pack_linear_witness(witness *wt) {
  size_t map[OUT_R], off[OUT_R], n[OUT_R];
  witness packed={0};
  outer_linear_tail_layout(wt->n,map,off,n);
  init_witness_raw(&packed,OUT_R-3,n);
  for(size_t i=0;i<OUT_R;i++) {
    memcpy(packed.s[map[i]]+off[i],wt->s[i],wt->n[i]*sizeof(poly));
  }
  for(size_t i=0;i<packed.r;i++)
    packed.normsq[i]=polyvec_sprodz(packed.s[i],packed.s[i],packed.n[i]);
  free_witness(wt);
  *wt=packed;
}

static int outer_pack_linear_statement(prncplstmnt *st) {
  size_t map[OUT_R], off[OUT_R], n[OUT_R];
  if(st->r!=OUT_R || st->build_postcommit_cnst) return 1;
  outer_linear_tail_layout(st->n,map,off,n);
  /* Refuse a future relation that uses any concatenated block quadratically. */
  for(size_t k=0;k<st->k;k++) {
    sparsemat *a=st->cnst[k].a;
    for(size_t j=0;j<a->len;j++) {
      const size_t x=a->rows[j], y=a->cols[j];
      if(map[x]==map[OUT_T0] || map[y]==map[OUT_T0]) return 1;
    }
  }
  for(size_t k=0;k<st->k;k++) {
    sparsecnst *c=&st->cnst[k];
    for(size_t j=0;j<c->nz;j++) {
      c->off[j]+=off[c->idx[j]];
      c->idx[j]=map[c->idx[j]];
    }
    for(size_t j=0;j<c->a->len;j++) {
      c->a->rows[j]=map[c->a->rows[j]];
      c->a->cols[j]=map[c->a->cols[j]];
    }
  }
  st->r=OUT_R-3;
  memcpy(st->n,n,st->r*sizeof(*n));
  shake128incctx hash;
  shake128_inc_init(&hash);
  absorb_bytes(&hash,"outer-linear-tail-concat-v1",sizeof("outer-linear-tail-concat-v1"));
  absorb_bytes(&hash,st->h,sizeof(st->h));
  for(size_t i=0;i<st->r;i++) absorb_u64(&hash,st->n[i]);
  shake128_inc_finalize(&hash);
  shake128_inc_squeeze(st->h,sizeof(st->h),&hash);
  return 0;
}

static int leopard_prepare_outer_witness(
    witness *wt, const outer_instance *inst,
    uint64_t *betasq_out) {
  const size_t tail_len =
      inst->ctx.rhat_len + inst->ctx.umc_eua_len +
      inst->ctx.umc_eub_len + inst->ctx.umc_r_len +
      inst->ctx.umc_va_len + inst->ctx.umc_vb_len;
  const size_t binary_len = inst->ctx.rhat_bits_len;
  const size_t p_len = inst->ctx.x_len;
  size_t n[OUT_R];
  size_t i;
  uint64_t betasq = 0;
  polx *p = NULL;
  __int128 *pcoeffs = NULL;
  poly *pdst[OUTER_DIRECT_P_DIGITS];
  polx *pbit = NULL;
  int64_t *pbitcoeffs = NULL;

  n[OUT_B0] = inst->ctx.b_len;
  n[OUT_B1] = inst->ctx.b_len;
  n[OUT_X0] = inst->ctx.x_len;
  n[OUT_X1] = inst->ctx.x_len;
  n[OUT_P0] = p_len;
  n[OUT_P1] = p_len;
  n[OUT_P2] = p_len;
  n[OUT_P3] = p_len;
  n[OUT_P4] = p_len;
  n[OUT_T0] = tail_len;
  n[OUT_T1] = tail_len;
  n[OUT_RBIT0] = binary_len;
  n[OUT_RBIT1] = binary_len;
  n[OUT_Z20] = inst->ctx.m2;
  n[OUT_Z21] = inst->ctx.m2;
  init_witness_raw(wt,OUT_R,n);

  p = alloc_polxvec(p_len);
  pcoeffs = _malloc(p_len*N*sizeof(*pcoeffs));
  compute_direct_outer_p(p,&inst->ctx,inst->z1);
  polxvec_centered_to_i128(pcoeffs,p,p_len);
  for(i=0;i<OUTER_DIRECT_P_DIGITS;i++)
    pdst[i] = wt->s[OUT_P0+i];
#if OUTER_RBIT_PRODUCT_DIGITS > 0
  for(i=0;i<binary_len;i++)
    polx_mul(&pbit[i],&inst->ctx.rbit_weight[i],
             &inst->z1[inst->ctx.off[LNP22_RHAT_BITS_BLOCK]+i]);
  polxvec_refresh(pbit,binary_len);
  polxvec_centered_to_i64(pbitcoeffs,pbit,binary_len);
  for(i=0;i<OUTER_RBIT_PRODUCT_DIGITS;i++)
    pbitdst[i] = wt->s[OUT_PBIT0+i];
#endif

  if(decompose_segment_two_digits(
         wt->s[OUT_B0],wt->s[OUT_B1],
         &inst->z1coeffs[inst->ctx.off[LNP22_B_BLOCK]*N],
         inst->ctx.b_len,OUTER_BZ_DECOMP_BASE,"B") ||
     decompose_direct_x_two_digits(
         wt->s[OUT_X0],wt->s[OUT_X1],
         &inst->z1coeffs[inst->ctx.off[LNP22_X_BLOCK]*N],
         &inst->ctx,OUTER_XY_DECOMP_BASE) ||
     decompose_segment_five_digits(
         pdst,pcoeffs,p_len,OUTER_DIRECT_P_BASE,"P") ||
     decompose_segment_two_digits(
         wt->s[OUT_T0],wt->s[OUT_T1],
         &inst->z1coeffs[inst->ctx.off[LNP22_RHAT_BLOCK]*N],
         tail_len,OUTER_XY_DECOMP_BASE,"tail") ||
     decompose_segment_two_digits(
         wt->s[OUT_RBIT0],wt->s[OUT_RBIT1],
         &inst->z1coeffs[
             inst->ctx.off[LNP22_RHAT_BITS_BLOCK]*N],
         binary_len,OUTER_XY_DECOMP_BASE,"RBIT") ||
#if OUTER_RBIT_PRODUCT_DIGITS > 0
     decompose_rbit_product_digits(
         pbitdst,pbitcoeffs,binary_len,"PBIT") ||
#endif
     decompose_segment_two_digits(
         wt->s[OUT_Z20],wt->s[OUT_Z21],inst->z2coeffs,
         inst->ctx.m2,OUTER_Z2_DECOMP_BASE,"z2")) {
    free(p);
    free(pcoeffs);
    free(pbit);
    free(pbitcoeffs);
    free_witness(wt);
    return 1;
  }

  for(i=0;i<OUT_R;i++) {
    wt->normsq[i] =
        polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    betasq += wt->normsq[i];
  }
  betasq += 16*N;
  *betasq_out = betasq;
  free(p);
  free(pcoeffs);
  free(pbit);
  free(pbitcoeffs);
  if(betasq > outer_admissible_betasq(&inst->ctx)) {
    free_witness(wt);
    return 1;
  }
  *betasq_out = outer_admissible_betasq(&inst->ctx);
  if(OUTER_PACK_LINEAR_TAIL) outer_pack_linear_witness(wt);
  return 0;
}

static int leopard_prepare_outer_statement(
    prncplstmnt *st,
    selector_lnp22_ctx *ctx,
    const polx *tA, const polx *w, const polx *c,
    const polx *vlin, const polx *tq, const polx *vq,
    uint64_t betasq,
    int prover_ctx,
    leopard_outer_statement_profile *prof) {
  polx *a3 = alloc_polxvec(ctx->m2);
  polx ufold_public[1];

  memset(prof,0,sizeof(*prof));
  compute_a3(a3,ctx);
  compute_public_fold_rhs(ufold_public,ctx);
  if(prover_ctx)
    build_direct_outer_statement(
        st,ctx,tA,w,vlin,tq,vq,c,a3,ufold_public,betasq,1);
  else
    build_selector_lnp22_outer_statement_public(
        st,ctx,tA,w,vlin,tq,vq,c,a3,ufold_public,betasq);
  free(a3);
  return OUTER_PACK_LINEAR_TAIL ? outer_pack_linear_statement(st) : 0;
}
#endif

static leopard_coeff *dup_coeffs(const leopard_coeff *src, size_t rows) {
  leopard_coeff *dst;

  dst = checked_malloc(rows * LEOPARD_RING_DEG * sizeof(*dst),"coeff copy");
  memcpy(dst,src,rows * LEOPARD_RING_DEG * sizeof(*dst));
  return dst;
}

static leopard_coeff *dup_slot_coeffs(const leopard_coeff *src, size_t slot_count, size_t rows) {
  leopard_coeff *dst;
  const size_t count = slot_count * rows * LEOPARD_RING_DEG;

  dst = checked_malloc(count * sizeof(*dst),"slot coeff copy");
  memcpy(dst,src,count * sizeof(*dst));
  return dst;
}

static int export_preproc_client_state_data(leopard_preq_data *out_preq,
                                            leopard_pre_state_data *out_st,
                                            const selector_lnp22_ctx *ctx) {
  const size_t slot_count = ctx->preproc_count;
  size_t j;

  memset(out_preq,0,sizeof(*out_preq));
  memset(out_st,0,sizeof(*out_st));
  out_preq->slot_count = slot_count;
  out_preq->cr_rows = ctx->cx_rows;
  out_preq->ar_seeds = checked_malloc(16 * slot_count,"preq ar seeds");
  out_preq->cr_coeffs = checked_malloc(slot_count * ctx->cx_rows * LEOPARD_RING_DEG *
                                       sizeof(*out_preq->cr_coeffs),"preq cr coeffs");

  out_st->slot_count = slot_count;
  out_st->next_j = 0;
  out_st->cr_rows = ctx->cx_rows;
  out_st->rhat_rows = ctx->rhat_len;
  out_st->eua_rows = ctx->umc_eua_len;
  out_st->eub_rows = ctx->umc_eub_len;
  out_st->r_rows = ctx->umc_r_len;
  out_st->va_rows = ctx->umc_va_len;
  out_st->vb_rows = ctx->umc_vb_len;
  out_st->ar_seeds = checked_malloc(16 * slot_count,"pre state ar seeds");
  out_st->cr_coeffs = checked_malloc(slot_count * ctx->cx_rows * LEOPARD_RING_DEG *
                                     sizeof(*out_st->cr_coeffs),"pre state cr coeffs");
  out_st->rhat_coeffs = checked_malloc(slot_count * ctx->rhat_len * LEOPARD_RING_DEG *
                                       sizeof(*out_st->rhat_coeffs),"pre state rhat coeffs");
  out_st->eua_coeffs = checked_malloc(slot_count * ctx->umc_eua_len * LEOPARD_RING_DEG *
                                      sizeof(*out_st->eua_coeffs),"pre state eua coeffs");
  out_st->eub_coeffs = checked_malloc(slot_count * ctx->umc_eub_len * LEOPARD_RING_DEG *
                                      sizeof(*out_st->eub_coeffs),"pre state eub coeffs");
  out_st->r_coeffs = checked_malloc(slot_count * ctx->umc_r_len * LEOPARD_RING_DEG *
                                    sizeof(*out_st->r_coeffs),"pre state r coeffs");
  out_st->va_coeffs = checked_malloc(slot_count * ctx->umc_va_len * LEOPARD_RING_DEG *
                                     sizeof(*out_st->va_coeffs),"pre state va coeffs");
  out_st->vb_coeffs = checked_malloc(slot_count * ctx->umc_vb_len * LEOPARD_RING_DEG *
                                     sizeof(*out_st->vb_coeffs),"pre state vb coeffs");

  memcpy(out_preq->ar_seeds,ctx->pre_ar_seed,16 * slot_count);
  memcpy(out_st->ar_seeds,ctx->pre_ar_seed,16 * slot_count);
  for(j=0;j<slot_count;j++) {
    export_polx_coeffs(&out_preq->cr_coeffs[j * ctx->cx_rows * LEOPARD_RING_DEG],
                       &ctx->pre_CX[j * ctx->cx_rows],ctx->cx_rows);
    export_polx_coeffs(&out_st->cr_coeffs[j * ctx->cx_rows * LEOPARD_RING_DEG],
                       &ctx->pre_CX[j * ctx->cx_rows],ctx->cx_rows);
    export_poly_coeffs(&out_st->rhat_coeffs[j * ctx->rhat_len * LEOPARD_RING_DEG],
                       &ctx->pre_rhat[j * ctx->rhat_len],ctx->rhat_len);
    export_poly_coeffs(&out_st->eua_coeffs[j * ctx->umc_eua_len * LEOPARD_RING_DEG],
                       &ctx->pre_eua[j * ctx->umc_eua_len],ctx->umc_eua_len);
    export_poly_coeffs(&out_st->eub_coeffs[j * ctx->umc_eub_len * LEOPARD_RING_DEG],
                       &ctx->pre_eub[j * ctx->umc_eub_len],ctx->umc_eub_len);
    export_poly_coeffs(&out_st->r_coeffs[j * ctx->umc_r_len * LEOPARD_RING_DEG],
                       &ctx->pre_r[j * ctx->umc_r_len],ctx->umc_r_len);
    export_poly_coeffs(&out_st->va_coeffs[j * ctx->umc_va_len * LEOPARD_RING_DEG],
                       &ctx->pre_va[j * ctx->umc_va_len],ctx->umc_va_len);
    export_poly_coeffs(&out_st->vb_coeffs[j * ctx->umc_vb_len * LEOPARD_RING_DEG],
                       &ctx->pre_vb[j * ctx->umc_vb_len],ctx->umc_vb_len);
  }
  return 0;
}

static int export_request_data(leopard_request_data *out,
                               const leopard_nizk1_request *req) {
  memset(out,0,sizeof(*out));
  out->slot_j = req->selected_j;
  out->cr_rows = req->cr_len;
  out->cx_rows = 1;
  out->cx_cols = req->cx_len;
  out->kmsis = req->kmsis;
  out->rbit_const_len = req->rbit_const_len;
  out->proof_len = req->proof_len;
  out->outer_betasq = req->outer_betasq;
  out->cr_coeffs = dup_coeffs(req->cr_coeffs,out->cr_rows);
  out->cx_coeffs = dup_coeffs(req->cx_coeffs,out->cx_rows * out->cx_cols);
  out->tA_coeffs = dup_coeffs(req->tA_coeffs,out->kmsis);
  out->w_coeffs = dup_coeffs(req->w_coeffs,out->kmsis);
  out->rbit_const_tg_coeffs = dup_coeffs(req->rbit_const_tg_coeffs,out->rbit_const_len);
  out->rbit_const_h_coeffs = dup_coeffs(req->rbit_const_h_coeffs,out->rbit_const_len);
  out->vlin_coeffs = dup_coeffs(req->vlin_coeffs,1);
  out->tq_coeffs = dup_coeffs(req->tq_coeffs,1);
  out->vq_coeffs = dup_coeffs(req->vq_coeffs,1);
  memcpy(out->ar_seed,req->ar_seed,sizeof(out->ar_seed));
  out->proof = checked_malloc(out->proof_len,"request outer proof");
  memcpy(out->proof,req->proof,out->proof_len);
  return 0;
}

static int export_state_data(leopard_state_data *out,
                             const leopard_nizk1_state *st) {
  memset(out,0,sizeof(*out));
  out->slot_j = st->selected_j;
  out->rhat_rows = st->rhat_len;
  out->rhat_coeffs = dup_coeffs(st->rhat_coeffs,out->rhat_rows);
  memcpy(out->tag_hash,st->tag_hash,sizeof(out->tag_hash));
  memcpy(out->xbits,st->xbits,sizeof(out->xbits));
  return 0;
}

static void request_data_to_internal(leopard_nizk1_request *dst,
                                     const leopard_request_data *src) {
  dst->selected_j = src->slot_j;
  import_polx_coeffs(dst->cr,src->cr_coeffs,src->cr_rows);
  import_polx_coeffs(dst->cx,src->cx_coeffs,src->cx_rows * src->cx_cols);
  import_polx_coeffs(dst->tA,src->tA_coeffs,src->kmsis);
  import_polx_coeffs(dst->w,src->w_coeffs,src->kmsis);
  import_polx_coeffs(dst->rbit_const_tg,src->rbit_const_tg_coeffs,src->rbit_const_len);
  import_polx_coeffs(dst->rbit_const_h,src->rbit_const_h_coeffs,src->rbit_const_len);
  import_polx_coeffs(dst->vlin,src->vlin_coeffs,1);
  import_polx_coeffs(dst->tq,src->tq_coeffs,1);
  import_polx_coeffs(dst->vq,src->vq_coeffs,1);
  memcpy(dst->cr_coeffs,src->cr_coeffs,src->cr_rows * LEOPARD_RING_DEG * sizeof(*dst->cr_coeffs));
  memcpy(dst->cx_coeffs,src->cx_coeffs,
         src->cx_rows * src->cx_cols * LEOPARD_RING_DEG * sizeof(*dst->cx_coeffs));
  memcpy(dst->tA_coeffs,src->tA_coeffs,src->kmsis * LEOPARD_RING_DEG * sizeof(*dst->tA_coeffs));
  memcpy(dst->w_coeffs,src->w_coeffs,src->kmsis * LEOPARD_RING_DEG * sizeof(*dst->w_coeffs));
  memcpy(dst->rbit_const_tg_coeffs,src->rbit_const_tg_coeffs,
         src->rbit_const_len * LEOPARD_RING_DEG * sizeof(*dst->rbit_const_tg_coeffs));
  memcpy(dst->rbit_const_h_coeffs,src->rbit_const_h_coeffs,
         src->rbit_const_len * LEOPARD_RING_DEG * sizeof(*dst->rbit_const_h_coeffs));
  memcpy(dst->vlin_coeffs,src->vlin_coeffs,LEOPARD_RING_DEG * sizeof(*dst->vlin_coeffs));
  memcpy(dst->tq_coeffs,src->tq_coeffs,LEOPARD_RING_DEG * sizeof(*dst->tq_coeffs));
  memcpy(dst->vq_coeffs,src->vq_coeffs,LEOPARD_RING_DEG * sizeof(*dst->vq_coeffs));
  memcpy(dst->ar_seed,src->ar_seed,sizeof(dst->ar_seed));
  dst->proof_len = src->proof_len;
  dst->outer_betasq = src->outer_betasq;
  dst->proof = checked_malloc(dst->proof_len,"request proof copy");
  memcpy(dst->proof,src->proof,dst->proof_len);
}

static void state_data_to_internal(leopard_nizk1_state *dst,
                                   const leopard_state_data *src) {
  dst->selected_j = src->slot_j;
  import_polx_coeffs(dst->rhat,src->rhat_coeffs,src->rhat_rows);
  memcpy(dst->tag_hash,src->tag_hash,sizeof(dst->tag_hash));
  memcpy(dst->xbits,src->xbits,sizeof(dst->xbits));
}

static int leopard_nizk1_request_prove(leopard_nizk1_request **out_req,
                                       leopard_nizk1_state **out_st,
                                       const leopard_pp *pp,
                                       const uint8_t *tag, size_t tag_len,
                                       const uint8_t xbytes[LEOPARD_X_BYTES]) {
  outer_instance inst = {};
  prncplstmnt outer_st = {};
  witness outer_wt = {};
  composite outer_pf = {};
  leopard_outer_statement_profile outer_stmt_prof = {};
  leopard_nizk1_request *req = NULL;
  leopard_nizk1_state *st = NULL;
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];
  uint64_t outer_betasq = 0;
  struct timespec t0, t1;

  *out_req = NULL;
  *out_st = NULL;
  memset(&g_last_request_profile,0,sizeof(g_last_request_profile));
  leopard_tag_hash(tag_hash,tag,tag_len);
  leopard_expand_input_bits(xbits,xbytes);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_generate_request_instance(&inst,pp,tag_hash,xbits))
    return 1;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_instance_sec = leopard_elapsed_sec(&t0,&t1);
#if !LEOPARD_REQUEST_SELF_CHECK
  g_last_request_profile.inner_selfcheck_sec = 0.0;
#endif

  req = leopard_request_alloc(inst.ctx.cx_rows,inst.ctx.m,inst.ctx.kmsis,
                              inst.ctx.rbit_const_lambda);
  st = leopard_state_alloc(inst.ctx.rhat_len);
  if(req == NULL || st == NULL)
    goto fail;
  polxvec_copy(st->rhat,&inst.ctx.s1[inst.ctx.off[LNP22_RHAT_BLOCK]],inst.ctx.rhat_len);
  export_polx_coeffs(st->rhat_coeffs,st->rhat,st->rhat_len);
  memcpy(st->tag_hash,tag_hash,sizeof(st->tag_hash));
  memcpy(st->xbits,xbits,sizeof(st->xbits));
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_outer_witness(&outer_wt,&inst,&outer_betasq) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_witness_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_outer_statement(&outer_st,&inst.ctx,inst.tA,inst.w,inst.c,
                                     inst.vlin,inst.tq,inst.vq,outer_betasq,
                                     1,&outer_stmt_prof) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_statement_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.outer_stmt_fold_sec = outer_stmt_prof.fold_sec;
  g_last_request_profile.outer_stmt_quad_sec = outer_stmt_prof.quad_sec;
  g_last_request_profile.outer_stmt_cnst1_sec = outer_stmt_prof.cnst1_sec;
  g_last_request_profile.outer_stmt_cnst2_sec = outer_stmt_prof.cnst2_sec;
  g_last_request_profile.outer_stmt_cnst3_sec = outer_stmt_prof.cnst3_sec;
  g_last_request_profile.outer_stmt_hash_sec = outer_stmt_prof.hash_sec;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(getenv("LEOPARD_PROOF_DIAGNOSTICS") && principle_verify(&outer_st,&outer_wt)) {
    fprintf(stderr,"request raw outer relation failed\n");
    goto fail;
  }
  release_inner_secret_buffers(&inst);
  if(composite_prove_principle(&outer_pf,&outer_st,&outer_wt) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_prove_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_serialize_composite(&req->proof,&req->proof_len,&outer_pf) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.serialize_sec = leopard_elapsed_sec(&t0,&t1);
  req->outer_betasq = outer_betasq;
  if(getenv("LEOPARD_PROOF_DIAGNOSTICS")) {
    composite copy={0};
    fprintf(stderr,"request original composite verify=%d\n",composite_verify_principle(&outer_pf,&outer_st));
    if(leopard_deserialize_composite(&copy,req->proof,req->proof_len)) goto fail;
    for(size_t level=0;level<copy.l;level++) {
      proof *a=outer_pf.pi[level], *b=copy.pi[level];
      const size_t total=a->cpp->u1len+a->cpp->u2len+LIFTS;
      size_t differences=0;
      for(size_t polyidx=0;polyidx<total;polyidx++) {
        polz diff;
        polz_sub(&diff,&a->u1[polyidx],&b->u1[polyidx]);
        polz_reduce(&diff); polz_center(&diff);
        if(!polz_iszero(&diff)) {
          differences++;
          fprintf(stderr,"roundtrip level=%zu poly=%zu top=%d/%d\n",level,polyidx,a->u1[polyidx].limbs[L-1].c[1],b->u1[polyidx].limbs[L-1].c[1]);
        }
      }
      fprintf(stderr,"roundtrip level=%zu poly diffs=%zu projection=%d\n",level,differences,memcmp(a->p,b->p,sizeof(a->p)));
    }
    fprintf(stderr,"request serialized composite verify=%d\n",composite_verify_principle(&copy,&outer_st));
    free_composite(&copy);
  }

  clock_gettime(CLOCK_MONOTONIC,&t0);
  polxvec_copy(req->cr,inst.ctx.CX,inst.ctx.cx_rows);
  polxvec_copy(req->cx,inst.ctx.Cout,inst.ctx.m);
  polxvec_copy(req->tA,inst.tA,inst.ctx.kmsis);
  polxvec_copy(req->w,inst.w,inst.ctx.kmsis);
  polxvec_copy(req->rbit_const_tg,inst.ctx.rbit_const_tg,inst.ctx.rbit_const_lambda);
  polxvec_copy(req->rbit_const_h,inst.ctx.rbit_const_h,inst.ctx.rbit_const_lambda);
  polxvec_copy(req->vlin,inst.vlin,1);
  polxvec_copy(req->tq,inst.tq,1);
  polxvec_copy(req->vq,inst.vq,1);
  export_polx_coeffs(req->cr_coeffs,req->cr,req->cr_len);
  export_polx_coeffs(req->cx_coeffs,req->cx,req->cx_len);
  export_polx_coeffs(req->tA_coeffs,req->tA,req->kmsis);
  export_polx_coeffs(req->w_coeffs,req->w,req->kmsis);
  export_polx_coeffs(req->rbit_const_tg_coeffs,req->rbit_const_tg,
                     req->rbit_const_len);
  export_polx_coeffs(req->rbit_const_h_coeffs,req->rbit_const_h,
                     req->rbit_const_len);
  export_polx_coeffs(req->vlin_coeffs,req->vlin,1);
  export_polx_coeffs(req->tq_coeffs,req->tq,1);
  export_polx_coeffs(req->vq_coeffs,req->vq,1);
  leopard_derive_ar_seed_from_cr(req->ar_seed,&inst.ctx,tag_hash);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.export_sec = leopard_elapsed_sec(&t0,&t1);

  free_prncplstmnt(&outer_st);
  free_witness(&outer_wt);
  free_composite(&outer_pf);
  free_outer_instance(&inst);
  *out_req = req;
  *out_st = st;
  return 0;

fail:
  free_prncplstmnt(&outer_st);
  free_witness(&outer_wt);
  free_composite(&outer_pf);
  free_outer_instance(&inst);
  leopard_nizk1_request_free(req);
  leopard_nizk1_state_free(st);
  return 1;
}

static int leopard_nizk1_request_prove_preprocessed(
    leopard_nizk1_request **out_req,
    leopard_nizk1_state **out_st,
    leopard_pre_state_data *io_pre_st,
    const leopard_pp *pp,
    const uint8_t *tag, size_t tag_len,
    const uint8_t xbytes[LEOPARD_X_BYTES]) {
  outer_instance inst = {};
  prncplstmnt outer_st = {};
  witness outer_wt = {};
  composite outer_pf = {};
  leopard_outer_statement_profile outer_stmt_prof = {};
  leopard_nizk1_request *req = NULL;
  leopard_nizk1_state *st = NULL;
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];
  uint64_t outer_betasq = 0;
  size_t slot_j;
  struct timespec t0, t1;

  if(io_pre_st == NULL || io_pre_st->slot_count == 0)
    return 1;

  *out_req = NULL;
  *out_st = NULL;
  memset(&g_last_request_profile,0,sizeof(g_last_request_profile));
  clock_gettime(CLOCK_MONOTONIC,&t0);
  leopard_tag_hash(tag_hash,tag,tag_len);
  leopard_expand_input_bits(xbits,xbytes);
  slot_j = io_pre_st->next_j % io_pre_st->slot_count;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_input_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_generate_request_instance_preprocessed(&inst,pp,io_pre_st,slot_j,
                                                    tag_hash,xbits) != 0)
    return 1;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.inner_instance_sec = leopard_elapsed_sec(&t0,&t1);
#if !LEOPARD_REQUEST_SELF_CHECK
  g_last_request_profile.inner_selfcheck_sec = 0.0;
#endif

  clock_gettime(CLOCK_MONOTONIC,&t0);
  req = leopard_request_alloc(inst.ctx.cx_rows,inst.ctx.m,inst.ctx.kmsis,
                              inst.ctx.rbit_const_lambda);
  st = leopard_state_alloc(inst.ctx.rhat_len);
  if(req == NULL || st == NULL)
    goto fail;
  req->selected_j = slot_j;
  st->selected_j = slot_j;
  polxvec_copy(st->rhat,&inst.ctx.s1[inst.ctx.off[LNP22_RHAT_BLOCK]],inst.ctx.rhat_len);
  export_polx_coeffs(st->rhat_coeffs,st->rhat,st->rhat_len);
  memcpy(st->tag_hash,tag_hash,sizeof(st->tag_hash));
  memcpy(st->xbits,xbits,sizeof(st->xbits));
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.protocol_state_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_outer_witness(&outer_wt,&inst,&outer_betasq) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_witness_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_outer_statement(&outer_st,&inst.ctx,inst.tA,inst.w,inst.c,
                                     inst.vlin,inst.tq,inst.vq,outer_betasq,
                                     1,&outer_stmt_prof) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_statement_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.outer_stmt_fold_sec = outer_stmt_prof.fold_sec;
  g_last_request_profile.outer_stmt_quad_sec = outer_stmt_prof.quad_sec;
  g_last_request_profile.outer_stmt_cnst1_sec = outer_stmt_prof.cnst1_sec;
  g_last_request_profile.outer_stmt_cnst2_sec = outer_stmt_prof.cnst2_sec;
  g_last_request_profile.outer_stmt_cnst3_sec = outer_stmt_prof.cnst3_sec;
  g_last_request_profile.outer_stmt_hash_sec = outer_stmt_prof.hash_sec;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(getenv("LEOPARD_PROOF_DIAGNOSTICS") && principle_verify(&outer_st,&outer_wt)) {
    fprintf(stderr,"request raw outer relation failed\n");
    goto fail;
  }
  release_inner_secret_buffers(&inst);
  if(composite_prove_principle(&outer_pf,&outer_st,&outer_wt) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.outer_prove_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_serialize_composite(&req->proof,&req->proof_len,&outer_pf) != 0)
    goto fail;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.serialize_sec = leopard_elapsed_sec(&t0,&t1);
  req->outer_betasq = outer_betasq;

  if(getenv("LEOPARD_PROOF_DIAGNOSTICS")) {
    composite copy={0};
    fprintf(stderr,"request original composite verify=%d\n",composite_verify_principle(&outer_pf,&outer_st));
    if(leopard_deserialize_composite(&copy,req->proof,req->proof_len)) goto fail;
    for(size_t level=0;level<copy.l;level++) {
      proof *a=outer_pf.pi[level], *b=copy.pi[level];
      const size_t total=a->cpp->u1len+a->cpp->u2len+LIFTS;
      size_t differences=0;
      for(size_t polyidx=0;polyidx<total;polyidx++) {
        polz diff;
        polz_sub(&diff,&a->u1[polyidx],&b->u1[polyidx]);
        polz_reduce(&diff); polz_center(&diff);
        if(!polz_iszero(&diff)) {
          differences++;
          fprintf(stderr,"roundtrip level=%zu poly=%zu top=%d/%d\n",level,polyidx,a->u1[polyidx].limbs[L-1].c[1],b->u1[polyidx].limbs[L-1].c[1]);
        }
      }
      fprintf(stderr,"roundtrip level=%zu poly diffs=%zu projection=%d\n",level,differences,memcmp(a->p,b->p,sizeof(a->p)));
    }
    fprintf(stderr,"request serialized composite verify=%d\n",composite_verify_principle(&copy,&outer_st));
    free_composite(&copy);
  }

  clock_gettime(CLOCK_MONOTONIC,&t0);
  polxvec_copy(req->cr,inst.ctx.CX,inst.ctx.cx_rows);
  polxvec_copy(req->cx,inst.ctx.Cout,inst.ctx.m);
  polxvec_copy(req->tA,inst.tA,inst.ctx.kmsis);
  polxvec_copy(req->w,inst.w,inst.ctx.kmsis);
  polxvec_copy(req->rbit_const_tg,inst.ctx.rbit_const_tg,inst.ctx.rbit_const_lambda);
  polxvec_copy(req->rbit_const_h,inst.ctx.rbit_const_h,inst.ctx.rbit_const_lambda);
  polxvec_copy(req->vlin,inst.vlin,1);
  polxvec_copy(req->tq,inst.tq,1);
  polxvec_copy(req->vq,inst.vq,1);
  export_polx_coeffs(req->cr_coeffs,req->cr,req->cr_len);
  export_polx_coeffs(req->cx_coeffs,req->cx,req->cx_len);
  export_polx_coeffs(req->tA_coeffs,req->tA,req->kmsis);
  export_polx_coeffs(req->w_coeffs,req->w,req->kmsis);
  export_polx_coeffs(req->rbit_const_tg_coeffs,req->rbit_const_tg,
                     req->rbit_const_len);
  export_polx_coeffs(req->rbit_const_h_coeffs,req->rbit_const_h,
                     req->rbit_const_len);
  export_polx_coeffs(req->vlin_coeffs,req->vlin,1);
  export_polx_coeffs(req->tq_coeffs,req->tq,1);
  export_polx_coeffs(req->vq_coeffs,req->vq,1);
  memcpy(req->ar_seed,&io_pre_st->ar_seeds[16 * slot_j],sizeof(req->ar_seed));
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_request_profile.export_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_request_profile.protocol_public_export_sec = g_last_request_profile.export_sec;

  io_pre_st->next_j = (slot_j + 1) % io_pre_st->slot_count;

  free_prncplstmnt(&outer_st);
  free_witness(&outer_wt);
  free_composite(&outer_pf);
  free_outer_instance(&inst);
  *out_req = req;
  *out_st = st;
  return 0;

fail:
  free_prncplstmnt(&outer_st);
  free_witness(&outer_wt);
  free_composite(&outer_pf);
  free_outer_instance(&inst);
  leopard_nizk1_request_free(req);
  leopard_nizk1_state_free(st);
  return 1;
}

static int leopard_nizk1_request_verify(const leopard_pp *pp,
                                        const uint8_t *tag, size_t tag_len,
                                        const leopard_nizk1_request *req) {
  selector_lnp22_ctx ctx;
  prncplstmnt outer_st = {};
  composite outer_pf = {};
  leopard_outer_statement_profile outer_stmt_prof = {};
  polx c_check[1];
  uint8_t ar_seed[16];
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  int ret = 1;
  struct timespec t0, t1;

  memset(&g_last_verify_profile,0,sizeof(g_last_verify_profile));
  if(req->cr_len != UMC_N+UMC_ELL+DEMO_M || req->cx_len != DEMO_M ||
     req->kmsis != LNP22_KMSIS || req->rbit_const_len != LNP22_RBIT_CONSTCHECK_LAMBDA)
    return 1;
  for(size_t i=0;i<req->cr_len*N;i++)
    if(req->cr_coeffs[i] < -(UMC_MODP/2) || req->cr_coeffs[i] > UMC_MODP/2)
      return 1;
  leopard_tag_hash(tag_hash,tag,tag_len);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_selector_public_ctx_snow(
         &ctx,pp,tag_hash,req->selected_j,
         req->cr,req->cx,
         req->rbit_const_tg,req->rbit_const_h,
         req->rbit_const_len,req->ar_seed,req->tA)) {
    selector_ctx_clear(&ctx);
    return 1;
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.public_ctx_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_verify_profile.a1open_cache = ctx.cache_a1open;
  if(req->cr_len != ctx.cx_rows || req->cx_len != ctx.m ||
     req->kmsis != ctx.kmsis || req->rbit_const_len != ctx.rbit_const_lambda) {
    selector_ctx_clear(&ctx);
    return 1;
  }

  if(req->outer_betasq != outer_admissible_betasq(&ctx))
    goto out;

  clock_gettime(CLOCK_MONOTONIC,&t0);
  fs_challenge(c_check,&ctx,req->tA,req->w,req->vlin,req->tq,req->vq);
  derive_ar_seed_from_commitment(ar_seed,&ctx,req->cr);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.fs_seed_sec = leopard_elapsed_sec(&t0,&t1);
  if(memcmp(ar_seed,req->ar_seed,sizeof(ar_seed)) != 0)
    goto out;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_deserialize_composite(&outer_pf,req->proof,req->proof_len) != 0)
    goto out;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.deserialize_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_prepare_outer_statement(&outer_st,&ctx,req->tA,req->w,c_check,
                                     req->vlin,req->tq,req->vq,req->outer_betasq,
                                     0,&outer_stmt_prof) != 0)
    goto out;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.outer_statement_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_verify_profile.outer_stmt_fold_sec = outer_stmt_prof.fold_sec;
  g_last_verify_profile.outer_stmt_quad_sec = outer_stmt_prof.quad_sec;
  g_last_verify_profile.outer_stmt_cnst1_sec = outer_stmt_prof.cnst1_sec;
  g_last_verify_profile.outer_stmt_cnst2_sec = outer_stmt_prof.cnst2_sec;
  g_last_verify_profile.outer_stmt_cnst3_sec = outer_stmt_prof.cnst3_sec;
  g_last_verify_profile.outer_stmt_hash_sec = outer_stmt_prof.hash_sec;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(composite_verify_principle(&outer_pf,&outer_st) != 0)
    goto out;
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.outer_verify_sec = leopard_elapsed_sec(&t0,&t1);
  ret = 0;

out:
  free_prncplstmnt(&outer_st);
  free_composite(&outer_pf);
  selector_ctx_clear(&ctx);
  return ret;
}

static size_t leopard_nizk1_request_bytes(const leopard_nizk1_request *req) {
  const size_t i16_bytes = 2;
  const size_t q_bytes = LEOPARD_Q_BYTES;

  return sizeof(uint32_t) +
         7 * sizeof(uint64_t) +
         sizeof(req->ar_seed) +
         req->cr_len * LEOPARD_RING_DEG * i16_bytes +
         req->cx_len * LEOPARD_RING_DEG * q_bytes +
         2 * req->kmsis * LEOPARD_RING_DEG * q_bytes +
         2 * req->rbit_const_len * LEOPARD_RING_DEG * q_bytes +
         3 * LEOPARD_RING_DEG * q_bytes +
         req->proof_len;
}

static size_t leopard_nizk1_cr_rows(const leopard_nizk1_request *req) {
  return req->cr_len;
}

static size_t leopard_nizk1_cx_rows(const leopard_nizk1_request *req) {
  return req->cx_len;
}

static size_t leopard_nizk1_rhat_rows(const leopard_nizk1_state *st) {
  return st->rhat_len;
}

static const leopard_coeff *leopard_nizk1_cr_coeffs(const leopard_nizk1_request *req) {
  return req->cr_coeffs;
}

static const leopard_coeff *leopard_nizk1_cx_coeffs(const leopard_nizk1_request *req) {
  return req->cx_coeffs;
}

static const leopard_coeff *leopard_nizk1_rhat_coeffs(const leopard_nizk1_state *st) {
  return st->rhat_coeffs;
}

static const uint8_t *leopard_nizk1_tag_hash(const leopard_nizk1_state *st) {
  return st->tag_hash;
}

static const uint8_t *leopard_nizk1_xbits(const leopard_nizk1_state *st) {
  return st->xbits;
}

static int leopard_nizk1_finalize_output(
    uint8_t y[LEOPARD_OUTPUT_BYTES],
    const leopard_nizk1_state *st,
    const leopard_coeff *ux_coeffs, size_t ux_rows,
    const leopard_coeff *vk_coeffs, size_t vk_rows) {
  polx *ux;
  polx *vk;
  polx rvk[1], z[1];

  if(ux_rows != 1 || vk_rows != st->rhat_len)
    return 1;

  ux = alloc_polxvec(ux_rows);
  vk = alloc_polxvec(vk_rows);
  import_polx_coeffs(ux,ux_coeffs,ux_rows);
  import_polx_coeffs(vk,vk_coeffs,vk_rows);
  polxvec_sprod(rvk,st->rhat,vk,st->rhat_len);
  polx_sub(z,ux,rvk);
  leopard_hash_output(y,st->tag_hash,st->xbits,z);
  free(ux);
  free(vk);
  return 0;
}

int leopard_nizk1_request_prove_data(leopard_request_data *out_req,
                                     leopard_state_data *out_st,
                                     const leopard_pp *pp,
                                     const uint8_t *tag, size_t tag_len,
                                     const uint8_t xbytes[LEOPARD_X_BYTES]) {
  leopard_nizk1_request *req = NULL;
  leopard_nizk1_state *st = NULL;
  int ret;

  ret = leopard_nizk1_request_prove(&req,&st,pp,tag,tag_len,xbytes);
  if(ret != 0)
    return ret;
  if(export_request_data(out_req,req) != 0) {
    leopard_nizk1_request_free(req);
    leopard_nizk1_state_free(st);
    return 1;
  }
  export_state_data(out_st,st);
  leopard_nizk1_request_free(req);
  leopard_nizk1_state_free(st);
  return 0;
}

int leopard_nizk1_preproc_client_data(leopard_preq_data *out_preq,
                                      leopard_pre_state_data *out_st,
                                      const leopard_pp *pp,
                                      size_t slot_count) {
  selector_lnp22_ctx ctx;
  size_t max_slots;

  selector_ctx_init(&ctx,pp->nizk1_seed);
  max_slots = ctx.preproc_count;
  if(slot_count == 0 || slot_count > max_slots) {
    selector_ctx_clear(&ctx);
    return 1;
  }
  ctx.preproc_count = slot_count;
  precompute_umc_commitments(&ctx);
  if(export_preproc_client_state_data(out_preq,out_st,&ctx) != 0) {
    selector_ctx_clear(&ctx);
    return 1;
  }
  selector_ctx_clear(&ctx);
  return 0;
}

int leopard_nizk1_request_prove_preprocessed_data(leopard_request_data *out_req,
                                                  leopard_state_data *out_st,
                                                  leopard_pre_state_data *io_pre_st,
                                                  const leopard_pp *pp,
                                                  const uint8_t *tag, size_t tag_len,
                                                  const uint8_t xbytes[LEOPARD_X_BYTES]) {
  leopard_nizk1_request *req = NULL;
  leopard_nizk1_state *st = NULL;
  int ret;

  ret = leopard_nizk1_request_prove_preprocessed(&req,&st,io_pre_st,pp,tag,
                                                 tag_len,xbytes);
  if(ret != 0)
    return ret;
  if(export_request_data(out_req,req) != 0) {
    leopard_nizk1_request_free(req);
    leopard_nizk1_state_free(st);
    return 1;
  }
  export_state_data(out_st,st);
  leopard_nizk1_request_free(req);
  leopard_nizk1_state_free(st);
  return 0;
}

int leopard_nizk1_request_fast_public_preprocessed_data(
    leopard_nizk1_fast_public_profile *out_prof,
    const leopard_pp *pp,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t *tag, size_t tag_len,
    const uint8_t xbytes[LEOPARD_X_BYTES],
    const leopard_request_data *expected_req) {
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];

  (void)pp;
  if(out_prof == NULL || pre == NULL || slot_j >= pre->slot_count)
    return 1;
  leopard_tag_hash(tag_hash,tag,tag_len);
  leopard_expand_input_bits(xbits,xbytes);
  return leopard_fast_public_compute_preprocessed(out_prof,pre,slot_j,tag_hash,
                                                  xbits,expected_req,NULL);
}

int leopard_nizk1_request_verify_data(const leopard_pp *pp,
                                      const uint8_t *tag, size_t tag_len,
                                      const leopard_request_data *req_data) {
  leopard_nizk1_request *req;
  int ret;
  if(req_data == NULL || req_data->cr_rows != UMC_N+UMC_ELL+DEMO_M ||
     req_data->cx_rows != 1 || req_data->cx_cols != DEMO_M ||
     req_data->kmsis != LNP22_KMSIS ||
     req_data->rbit_const_len != LNP22_RBIT_CONSTCHECK_LAMBDA ||
     req_data->cr_coeffs == NULL)
    return 1;
  for(size_t i=0;i<req_data->cr_rows*N;i++)
    if(req_data->cr_coeffs[i] < -(UMC_MODP/2) || req_data->cr_coeffs[i] > UMC_MODP/2)
      return 1;
  size_t cx_len = req_data->cx_rows * req_data->cx_cols;

  req = leopard_request_alloc(req_data->cr_rows,cx_len,req_data->kmsis,
                              req_data->rbit_const_len);
  request_data_to_internal(req,req_data);
  ret = leopard_nizk1_request_verify(pp,tag,tag_len,req);
  leopard_nizk1_request_free(req);
  return ret;
}

size_t leopard_nizk1_request_data_bytes(const leopard_request_data *req) {
  leopard_nizk1_request fake = {0};

  fake.cr_len = req->cr_rows;
  fake.cx_len = req->cx_rows * req->cx_cols;
  fake.kmsis = req->kmsis;
  fake.rbit_const_len = req->rbit_const_len;
  fake.proof_len = req->proof_len;
  fake.outer_betasq = req->outer_betasq;
  memcpy(fake.ar_seed,req->ar_seed,sizeof(fake.ar_seed));
  return leopard_nizk1_request_bytes(&fake);
}

int leopard_nizk1_finalize_output_data(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                       const leopard_state_data *st_data,
                                       const leopard_reply_data *rep) {
  (void)y;
  (void)st_data;
  (void)rep;
  return 1;
}

int leopard_nizk1_finalize_raw_output_data(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                           const leopard_state_data *st_data,
                                           const leopard_coeff *ux_coeffs,
                                           size_t ux_rows,
                                           const leopard_coeff *vk_coeffs,
                                           size_t vk_rows) {
  return leopard_finalize_output_binary_rhat(y,st_data,ux_coeffs,ux_rows,
                                             vk_coeffs,vk_rows);
}

int leopard_nizk1_eval(uint8_t y[LEOPARD_OUTPUT_BYTES],
                       const leopard_pp *pp,
                       const uint8_t *tag, size_t tag_len,
                       const uint8_t xbytes[LEOPARD_X_BYTES],
                       const leopard_coeff *k_coeffs, size_t k_rows) {
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];
  polx *k;
  polx *bxt;
  polx z[1];

  (void)pp;
  leopard_tag_hash(tag_hash,tag,tag_len);
  leopard_expand_input_bits(xbits,xbytes);
  if(k_rows != DEMO_M)
    return 1;
  k = alloc_polxvec(k_rows);
  bxt = alloc_polxvec(k_rows);
  if(leopard_compute_bxt_fast(bxt,tag_hash,xbits) != 0) {
    free(k);
    free(bxt);
    return 1;
  }
  import_polx_coeffs(k,k_coeffs,k_rows);
  polxvec_sprod(z,bxt,k,k_rows);
  leopard_hash_output(y,tag_hash,xbits,z);
  free(k);
  free(bxt);
  return 0;
}

const leopard_nizk1_request_profile *leopard_nizk1_last_request_profile(void) {
  return &g_last_request_profile;
}

const leopard_nizk1_verify_profile *leopard_nizk1_last_verify_profile(void) {
  return &g_last_verify_profile;
}
