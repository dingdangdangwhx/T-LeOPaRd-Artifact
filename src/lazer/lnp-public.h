#ifndef LAZER_LNP_PUBLIC_H
#define LAZER_LNP_PUBLIC_H

#include "lazer.h"
#include "memory.h"
#include "brandom.h"

#ifndef LAZER_SZ_SPLIT_ACCUM
#define LAZER_SZ_SPLIT_ACCUM 1
#endif
#ifndef LAZER_Q72_FAST
#define LAZER_Q72_FAST 1
#endif
#ifndef LAZER_PUBLIC_SPARSE_PRODUCTS
#define LAZER_PUBLIC_SPARSE_PRODUCTS 1
#endif

#if LAZER_PUBLIC_SPARSE_PRODUCTS
/* Only call on public statement matrices, never on secret witness entries. */
static int
__lnp_public_poly_is_zero (poly_t p)
{
  intvec_srcptr coefficients = poly_get_coeffvec (p);
  for (unsigned int i = 0; i < coefficients->nelems; i++)
    {
      int_srcptr c = intvec_get_elem_src (coefficients, i);
      for (unsigned int j = 0; j < c->nlimbs; j++)
        if (c->limbs[j] != 0)
          return 0;
    }
  return 1;
}

static void
__lnp_public_addmul (polyvec_t out, polymat_t matrix, polyvec_t witness)
{
  poly_t product, sum;
  poly_alloc (product, matrix->ring);
  poly_alloc (sum, matrix->ring);
  for (unsigned int i = 0; i < matrix->nrows; i++)
    {
      poly_set_zero (sum);
      for (unsigned int j = 0; j < matrix->ncols; j++)
        {
          poly_ptr entry = polymat_get_elem (matrix, i, j);
          if (__lnp_public_poly_is_zero (entry))
            continue;
          poly_mul (product, entry, polyvec_get_elem (witness, j));
          poly_add (sum, sum, product, 1);
        }
      poly_ptr dst = polyvec_get_elem (out, i);
      poly_add (dst, dst, sum, 0);
    }
  poly_free (sum);
  poly_free (product);
}

static void
__lnp_public_leftmul (polymat_t out, polymat_t left, polymat_t matrix)
{
  poly_t product;
  poly_alloc (product, matrix->ring);
  polymat_set_zero (out);
  /* Check each public entry once, before multiplication changes its CRT cache. */
  for (unsigned int j = 0; j < matrix->ncols; j++)
    for (unsigned int i = 0; i < matrix->nrows; i++)
      {
        poly_ptr entry = polymat_get_elem (matrix, i, j);
        if (__lnp_public_poly_is_zero (entry))
          continue;
        for (unsigned int k = 0; k < left->nrows; k++)
          {
            poly_ptr dst = polymat_get_elem (out, k, j);
            poly_mul (product, polymat_get_elem (left, k, i), entry);
            poly_add (dst, dst, product, 1);
          }
      }
  poly_free (product);
}

/* Branches depend only on public matrix coefficients, not on the witness. */
static void
__lnp_public_leftvec (polyvec_t out, polyvec_t left, polymat_t matrix)
{
  poly_t product;
  poly_alloc (product, matrix->ring);
  polyvec_set_zero (out);
  for (unsigned int j = 0; j < matrix->ncols; j++)
    for (unsigned int i = 0; i < matrix->nrows; i++)
      {
        poly_ptr entry = polymat_get_elem (matrix, i, j);
        if (__lnp_public_poly_is_zero (entry))
          continue;
        poly_ptr dst = polyvec_get_elem (out, j);
        poly_mul (product, polyvec_get_elem (left, i), entry);
        poly_add (dst, dst, product, 1);
      }
  poly_free (product);
}

static void
__lnp_public_auto (polymat_t out, polymat_t matrix)
{
  for (unsigned int i = 0; i < matrix->nrows; i++)
    for (unsigned int j = 0; j < matrix->ncols; j++)
      {
        poly_ptr src = polymat_get_elem (matrix, i, j);
        poly_ptr dst = polymat_get_elem (out, i, j);
        if (__lnp_public_poly_is_zero (src))
          poly_set_zero (dst);
        else
          poly_auto (dst, src);
      }
}

static int
__lnp_public_constant_matrix (polymat_t matrix)
{
  for (unsigned int i = 0; i < matrix->nrows; i++)
    for (unsigned int j = 0; j < matrix->ncols; j++)
      {
        intvec_srcptr coefficients = poly_get_coeffvec (polymat_get_elem (matrix, i, j));
        for (unsigned int k = 1; k < coefficients->nelems; k++)
          {
            int_srcptr c = intvec_get_elem_src (coefficients, k);
            for (unsigned int limb = 0; limb < c->nlimbs; limb++)
              if (c->limbs[limb] != 0)
                return 0;
          }
      }
  return 1;
}

/* Encode exactly the dense public matrix, including every zero coefficient.
 * This changes neither the hashed statement nor its canonical byte encoding. */
static void
__lnp_encode_public_matrix (coder_state_t state, polymat_t matrix,
                            const int_t q, unsigned int log2q)
{
  const size_t bits = (size_t)matrix->ring->d * log2q;
  for (unsigned int i = 0; i < matrix->nrows; i++)
    for (unsigned int j = 0; j < matrix->ncols; j++)
      {
        poly_ptr p = polymat_get_elem (matrix, i, j);
        if (state->bit_off == 0 && (bits & 7) == 0
            && __lnp_public_poly_is_zero (p))
          {
            memset (state->out, 0, bits / 8);
            state->out += bits / 8;
            state->byte_off += bits / 8;
          }
        else
          {
            poly_redp (p, p);
            coder_enc_urandom2 (state, p, q, log2q);
          }
      }
}
#endif

#if LAZER_SZ_SPLIT_ACCUM
static void
__sz_ternary_product (intmat_t out, const intmat_t V, const int_t q,
                      unsigned int rows, unsigned int cols,
                      const uint8_t seed[32], unsigned int domain_offset)
{
  const size_t len = (size_t)rows * cols;
  int64_t *restrict low = _alloc (len * sizeof (*low));
  int64_t *restrict high = _alloc (len * sizeof (*high));
  int8_t *signs = _alloc (cols);
  const unsigned __int128 modulus = (unsigned __int128)q->limbs[0]
    | (q->nlimbs > 1 ? (unsigned __int128)q->limbs[1] << 64 : 0);

  memset (low, 0, len * sizeof (*low));
  memset (high, 0, len * sizeof (*high));
  for (unsigned int i = 0; i < 256; i++)
    {
      _brandom (signs, cols, 1, seed, domain_offset + i);
      for (unsigned int k = 0; k < rows; k++)
        {
          int_srcptr challenge = intmat_get_elem_src (V, k, i);
          const int64_t sign = 1 - 2 * (int64_t)challenge->neg;
          const int64_t lo = sign * (int64_t)(challenge->limbs[0] & UINT32_MAX);
          const uint64_t upper = (challenge->limbs[0] >> 32)
            | (challenge->nlimbs > 1 ? (uint64_t)challenge->limbs[1] << 32 : 0);
          const int64_t hi = sign * (int64_t)upper;
          int64_t *restrict alo = low + (size_t)k * cols;
          int64_t *restrict ahi = high + (size_t)k * cols;

          /* For q <= 2^72, the high-part sum is below 2^48, low below 2^40.
             Avoid int128 operations and ternary branches in the SIMD loop. */
          for (unsigned int j = 0; j < cols; j++)
            {
              alo[j] += lo * signs[j];
              ahi[j] += hi * signs[j];
            }
        }
    }

  for (unsigned int k = 0; k < rows; k++)
    for (unsigned int j = 0; j < cols; j++)
      {
        const size_t pos = (size_t)k * cols + j;
        const __int128 value = (__int128)high[pos] * ((__int128)1 << 32) + low[pos];
        __int128 rem = value % (__int128)modulus;
        int_ptr dst = intmat_get_elem (out, k, j);
        if (rem < 0)
          rem += (__int128)modulus;
        dst->neg = 0;
        dst->limbs[0] = (limb_t)rem;
        for (unsigned int limbidx = 1; limbidx < dst->nlimbs; limbidx++)
          dst->limbs[limbidx] = 0;
        if (dst->nlimbs > 1)
          dst->limbs[1] = (limb_t)((unsigned __int128)rem >> 64);
      }
  _free (signs, cols);
  _free (high, len * sizeof (*high));
  _free (low, len * sizeof (*low));
}
#endif

#endif

