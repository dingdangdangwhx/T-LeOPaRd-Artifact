#include "lazer.h"

#ifndef NIZK2_PARAMS_HEADER
#define NIZK2_PARAMS_HEADER "params.h"
#endif
#include NIZK2_PARAMS_HEADER

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "commitment.h"

#ifndef NIZK2_M
#define NIZK2_M 27
#endif
#ifndef NIZK2_N
#define NIZK2_N 24
#endif
#ifndef NIZK2_H
#define NIZK2_H 1
#endif
#ifndef NIZK2_D
#define NIZK2_D 64
#endif
#ifndef NIZK2_MOD_Q
#define NIZK2_MOD_Q 281474976710597LL
#endif
#ifndef NIZK2_ES_PRIME_MAX
#define NIZK2_ES_PRIME_MAX 7781
#endif
#ifndef NIZK2_DECOMP_BASE
#define NIZK2_DECOMP_BASE 100
#endif
#ifndef NIZK2_DECOMP_LEN
#define NIZK2_DECOMP_LEN 2
#endif
#ifndef NIZK2_AGGREGATE_MAX
#define NIZK2_AGGREGATE_MAX 1
#endif
#ifndef NIZK2_PROFILE_TIMERS
#define NIZK2_PROFILE_TIMERS 0
#endif
#ifndef NIZK2_LEGACY_I64_INSTANCE
#define NIZK2_LEGACY_I64_INSTANCE 0
#endif
#ifndef NIZK2_CHECK_LINEAR_RELATION
#define NIZK2_CHECK_LINEAR_RELATION 1
#endif

enum
{
  M = NIZK2_M,
  N = NIZK2_N,
  H = NIZK2_H,
  D = NIZK2_D,
  ES_PRIME_MAX = NIZK2_ES_PRIME_MAX,
  DECOMP_BASE = NIZK2_DECOMP_BASE,
  DECOMP_LEN = NIZK2_DECOMP_LEN,
  COM_RANK = NIZK2_SMALL_COMMITMENT ? NIZK2_COM_RANK : NIZK2_M,
  COM_MESSAGE_ROWS = NIZK2_N * (NIZK2_DOUBLE_COMMITMENT ? 2 : 1),
  COM_ROWS = COM_RANK + COM_MESSAGE_ROWS,
  WIT_KEY_OFF = 2 * COM_RANK + COM_MESSAGE_ROWS,
  WIT_ES_OFF = WIT_KEY_OFF + NIZK2_N,
  WIT_DIGIT_OFF = WIT_ES_OFF + NIZK2_M + NIZK2_N,
  REL_ROWS = COM_ROWS + NIZK2_M + NIZK2_N + NIZK2_H,
  WIT_ROWS = WIT_DIGIT_OFF + NIZK2_H * NIZK2_DECOMP_LEN,
};

#include "commitment_proof.h"

static double get_sec (void);

static void
sample_polyvec_bnd (polyvec_t v, int64_t lo64, int64_t hi64,
                    const uint8_t seed[32], uint32_t dom)
{
  int_srcptr q = polyring_get_mod (polyvec_get_ring (v));
  INT_T (lo, q->nlimbs);
  INT_T (hi, q->nlimbs);

  int_set_i64 (lo, lo64);
  int_set_i64 (hi, hi64);
  polyvec_urandom_bnd (v, lo, hi, seed, dom);
}

#if NIZK2_LEGACY_I64_INSTANCE || NIZK2_SMALL_COMMITMENT
static void
sample_polymat_bnd (polymat_t mat, int64_t lo64, int64_t hi64,
                    const uint8_t seed[32], uint32_t dom)
{
  unsigned int i;
  polyvec_t row;

  for (i = 0; i < polymat_get_nrows (mat); i++)
    {
      polymat_get_row (row, mat, i);
      sample_polyvec_bnd (row, lo64, hi64, seed, dom + i);
    }
}
#endif

static void
set_const_poly (poly_t poly, int64_t c)
{
  int_srcptr q = polyring_get_mod (poly->ring);
  INT_T (val, q->nlimbs);

  int_set_i64 (val, c);
  int_set (poly_get_coeff (poly, 0), val);
}

#if !NIZK2_LEGACY_I64_INSTANCE
static void
sample_polyvec_centered_mod (polyvec_t v, const uint8_t seed[32],
                             uint32_t dom)
{
  int_srcptr q = polyring_get_mod (polyvec_get_ring (v));
  INT_T (lo, q->nlimbs);
  INT_T (hi, q->nlimbs);

  int_set (hi, q);
  int_rshift (hi, hi, 1);
  int_neg (lo, hi);
  polyvec_urandom_bnd (v, lo, hi, seed, dom);
}

static void
sample_polymat_centered_mod (polymat_t mat, const uint8_t seed[32],
                             uint32_t dom)
{
  unsigned int i;
  polyvec_t row;

  for (i = 0; i < polymat_get_nrows (mat); i++)
    {
      polymat_get_row (row, mat, i);
      sample_polyvec_centered_mod (row, seed, dom + i);
    }
}

static void
set_const_poly_int (poly_t poly, const int_t c)
{
  int_set (poly_get_coeff (poly, 0), c);
}

static void
set_const_diag_int (polymat_t mat, const int_t c)
{
  unsigned int i;
  const unsigned int nrows = polymat_get_nrows (mat);
  const unsigned int ncols = polymat_get_ncols (mat);
  const unsigned int n = nrows < ncols ? nrows : ncols;

  polymat_set_zero (mat);
  for (i = 0; i < n; i++)
    set_const_poly_int (polymat_get_elem (mat, i, i), c);
}

static void
set_const_diag_q_div_pow2 (polymat_t mat, unsigned int logdiv)
{
  int_srcptr q = polyring_get_mod (mat->ring);
  INT_T (c, q->nlimbs);

  int_set (c, q);
  int_rshift (c, c, logdiv);
  set_const_diag_int (mat, c);
}

#endif

#if NIZK2_LEGACY_I64_INSTANCE || NIZK2_SMALL_COMMITMENT
static void
set_const_diag (polymat_t mat, int64_t c)
{
  unsigned int i;
  const unsigned int nrows = polymat_get_nrows (mat);
  const unsigned int ncols = polymat_get_ncols (mat);
  const unsigned int n = nrows < ncols ? nrows : ncols;

  polymat_set_zero (mat);
  for (i = 0; i < n; i++)
    set_const_poly (polymat_get_elem (mat, i, i), c);
}
#endif

static void
set_gi (polymat_t gi)
{
  unsigned int i, j;
  int64_t factor = 1;

  polymat_set_zero (gi);
  for (j = 0; j < DECOMP_LEN; j++)
    {
      for (i = 0; i < H; i++)
        set_const_poly (polymat_get_elem (gi, i, j * H + i), factor);
      factor *= DECOMP_BASE;
    }
}

static void
decompose_es_prime (polyvec_t es_prime_dec, polyvec_t es_prime)
{
  int64_t raw[D * H];
  int64_t dec[D * H * DECOMP_LEN];
  unsigned int idx;

  for (idx = 0; idx < D * H * DECOMP_LEN; idx++)
    dec[idx] = 0;

  polyvec_get_coeffvec_i64 (raw, es_prime);
#ifdef NIZK2_GAUSSIAN_ERRORS
  _Static_assert(DECOMP_LEN == 2,"Gaussian error encoding needs two digits");
  for(idx=0;idx<D*H;idx++) {
    const int64_t x=raw[idx];
    const int64_t high=x>=0 ? (x+DECOMP_BASE/2)/DECOMP_BASE : -((-x+DECOMP_BASE/2)/DECOMP_BASE);
    dec[idx]=x-DECOMP_BASE*high;
    dec[D*H+idx]=high;
  }
#else
  for (idx = 0; idx < D * H; idx++)
    {
      uint64_t v = (uint64_t)raw[idx];

      for (unsigned int j = 0; j < DECOMP_LEN; j++)
        {
          dec[j * (D * H) + idx] = (int64_t)(v % DECOMP_BASE);
          v /= DECOMP_BASE;
        }
    }

#endif
  polyvec_set_coeffvec_i64 (es_prime_dec, dec);
}


static void
generate_instance (polymat_t A, polyvec_t s, polyvec_t t)
{
  const uint8_t instseed[32] = { 0 };
#if NIZK2_LEGACY_I64_INSTANCE
  const int64_t q = NIZK2_MOD_Q;
  const int64_t qhalf = q / 2;
  const int64_t q16 = q / 16;
#endif

  polymat_t A1, Id, B, q16I, Ar, Cx, GI;
  polyvec_t s_r, s_e, s_ei, s_k, s_es, s_esp, es_prime_raw;
  polyvec_t t_c, t_cu, t_vk, t_uk, tmp_n, tmp_h;

  polymat_set_zero (A);
  polyvec_set_zero (s);
  polyvec_set_zero (t);

  /* Witness layout: s = [r, e, e_i, k, e_s, es_prime_dec]. */
  polyvec_get_subvec (s_r, s, 0, M, 1);
  polyvec_get_subvec (s_e, s, M, M, 1);
  polyvec_get_subvec (s_ei, s, 2 * M, N, 1);
  polyvec_get_subvec (s_k, s, 2 * M + N, N, 1);
  polyvec_get_subvec (s_es, s, 2 * M + 2 * N, M + N, 1);
  polyvec_get_subvec (s_esp, s, (M + N) * 3, H * DECOMP_LEN, 1);

  /* Public statement block layout. */
  polymat_get_submat (A1, A, 0, 0, M, M, 1, 1);
  polymat_get_submat (Id, A, 0, M, M, M, 1, 1);
#if NIZK2_LEGACY_I64_INSTANCE
  sample_polymat_bnd (A1, -qhalf, qhalf, instseed, 0);
#else
  sample_polymat_centered_mod (A1, instseed, 0);
#endif
  polymat_set_one (Id);

  polymat_get_submat (B, A, M, 0, N, M, 1, 1);
  polymat_get_submat (q16I, A, M, 2 * M + N, N, N, 1, 1);
  polymat_get_submat (Id, A, M, 2 * M, N, N, 1, 1);
#if NIZK2_LEGACY_I64_INSTANCE
  sample_polymat_bnd (B, -qhalf, qhalf, instseed, 100);
  set_const_diag (q16I, q16);
#else
  sample_polymat_centered_mod (B, instseed, 100);
  set_const_diag_q_div_pow2 (q16I, 4);
#endif
  polymat_set_one (Id);

  polymat_get_submat (Ar, A, M + N, 2 * M + N, M + N, N, 1, 1);
  polymat_get_submat (Id, A, M + N, 2 * M + 2 * N, M + N, M + N, 1, 1);
#if NIZK2_LEGACY_I64_INSTANCE
  sample_polymat_bnd (Ar, -qhalf, qhalf, instseed, 200);
#else
  sample_polymat_centered_mod (Ar, instseed, 200);
#endif
  polymat_set_one (Id);

  polymat_get_submat (Cx, A, (M + N) * 2, 2 * M + N, H, N, 1, 1);
  polymat_get_submat (GI, A, (M + N) * 2, (M + N) * 3, H, H * DECOMP_LEN, 1,
                      1);
#if NIZK2_LEGACY_I64_INSTANCE
  sample_polymat_bnd (Cx, -qhalf, qhalf, instseed, 300);
#else
  sample_polymat_centered_mod (Cx, instseed, 300);
#endif
  set_gi (GI);

  /* Sample witness blocks following the Sage generator ranges. */
  sample_polyvec_bnd (s_r, 0, NIZK2_AGGREGATE_MAX, instseed, 1000);
  sample_polyvec_bnd (s_e, 0, NIZK2_AGGREGATE_MAX, instseed, 1001);
  sample_polyvec_bnd (s_ei, 0, NIZK2_AGGREGATE_MAX, instseed, 1002);
  sample_polyvec_bnd (s_k, 0, 15 * NIZK2_AGGREGATE_MAX, instseed, 1003);
  sample_polyvec_bnd (s_es, 0, 14 * NIZK2_AGGREGATE_MAX, instseed, 1004);

  polyvec_alloc (es_prime_raw, polyvec_get_ring (s), H);
  sample_polyvec_bnd (es_prime_raw, 0, ES_PRIME_MAX, instseed, 1005);
  decompose_es_prime (s_esp, es_prime_raw);
  polyvec_free (es_prime_raw);

  /* Compute t = -(c, cu, vk, uk). */
  polyvec_get_subvec (t_c, t, 0, M, 1);
  polyvec_get_subvec (t_cu, t, M, N, 1);
  polyvec_get_subvec (t_vk, t, M + N, M + N, 1);
  polyvec_get_subvec (t_uk, t, (M + N) * 2, H, 1);

  polyvec_mul (t_c, A1, s_r);
  polyvec_add (t_c, t_c, s_e, 0);

  polyvec_alloc (tmp_n, polyvec_get_ring (s), N);
  polyvec_mul (t_cu, B, s_r);
  polyvec_mul (tmp_n, q16I, s_k);
  polyvec_add (t_cu, t_cu, tmp_n, 0);
  polyvec_add (t_cu, t_cu, s_ei, 0);
  polyvec_free (tmp_n);

  polyvec_mul (t_vk, Ar, s_k);
  polyvec_add (t_vk, t_vk, s_es, 0);

  polyvec_alloc (tmp_h, polyvec_get_ring (s), H);
  polyvec_mul (t_uk, Cx, s_k);
  polyvec_mul (tmp_h, GI, s_esp);
  polyvec_add (t_uk, t_uk, tmp_h, 0);
  polyvec_free (tmp_h);

  polyvec_mod (t, t);
  polyvec_redc (t, t);
  polyvec_neg_self (t);
  polyvec_mod (t, t);
  polyvec_redc (t, t);
}

static size_t
prover (uint8_t *proof, polyvec_t s, polymat_t A, polyvec_t t,
        const uint8_t pp[32])
{
  return mixed_prover(proof,s,A,t,pp);
}

#if NIZK2_CHECK_LINEAR_RELATION
static int
check_linear_relation (polymat_t A, polyvec_t s, polyvec_t t)
{
  const unsigned int nrels = polyvec_get_nelems (t);
  const size_t ncoeffs = (size_t)nrels * D;
  polyvec_t lhs;
  int64_t *raw;
  size_t i;
  size_t nonzero = 0;
  int64_t maxabs = 0;

  raw = malloc (ncoeffs * sizeof (*raw));
  if (raw == NULL)
    return 1;
  polyvec_alloc (lhs, polyvec_get_ring (t), nrels);
  polyvec_mul (lhs, A, s);
  polyvec_add (lhs, lhs, t, 0);
  polyvec_mod (lhs, lhs);
  polyvec_redc (lhs, lhs);
  polyvec_get_coeffvec_i64 (raw, lhs);
  for (i = 0; i < ncoeffs; i++)
    {
      int64_t a = raw[i] < 0 ? -raw[i] : raw[i];
      if (raw[i] != 0)
        nonzero++;
      if (a > maxabs)
        maxabs = a;
    }
  if (nonzero != 0)
    printf ("linear relation check failed: nonzero=%zu maxabs=%lld\n",
            nonzero, (long long)maxabs);
  polyvec_free (lhs);
  free (raw);
  return nonzero == 0 ? 0 : 1;
}
#endif

static int
verifier (const uint8_t *proof, polymat_t A, polyvec_t t, const uint8_t pp[32])
{
  return mixed_verifier(proof,A,t,pp);
}

static double
get_sec (void)
{
  struct timeval tv;

  gettimeofday (&tv, NULL);
  return ((double)tv.tv_sec) + 1E-6 * ((double)tv.tv_usec);
}

