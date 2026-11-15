#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../third_party/nizk1/src/core/aesctr.h"

#include "common.h"
#include "nizk2.h"
#include "wire.h"

#include "../../third_party/nizk2/server_proof.c"

static void *checked_malloc(size_t size, const char *what) {
  void *ptr = malloc(size);

  if(ptr == NULL) {
    fprintf(stderr,"allocation failed for %s (%zu bytes)\n",what,size);
    exit(1);
  }
  return ptr;
}

static leopard_coeff *alloc_coeffs(size_t rows) {
  leopard_coeff *coeffs;

  coeffs = checked_malloc(rows * LEOPARD_RING_DEG * sizeof(*coeffs),"coefficients");
  memset(coeffs,0,rows * LEOPARD_RING_DEG * sizeof(*coeffs));
  return coeffs;
}

static _Thread_local leopard_nizk2_blindeval_profile g_last_blindeval_profile;
static _Thread_local leopard_nizk2_verify_profile g_last_verify_profile;

static double leopard_elapsed_sec(const struct timespec *start,
                                  const struct timespec *end) {
  return (double)(end->tv_sec - start->tv_sec) +
         1e-9 * (double)(end->tv_nsec - start->tv_nsec);
}

static void leopard_shake_absorb_u64(shake128_state_t st, uint64_t x) {
  uint8_t buf[8];
  size_t i;

  for(i=0;i<sizeof(buf);i++)
    buf[i] = (uint8_t)(x >> (8 * i));
  shake128_absorb(st,buf,sizeof(buf));
}

static void leopard_shake_absorb_coeffs(shake128_state_t st,
                                        const leopard_coeff *coeffs, size_t rows) {
  if(coeffs != NULL && rows != 0)
    shake128_absorb(st,(const uint8_t *)coeffs,
                    rows * LEOPARD_RING_DEG * sizeof(*coeffs));
}

static void leopard_noise_seed(uint8_t seed[32], const leopard_pp *pp,
                               const leopard_sk_data *sk,
                               const uint8_t ar_seed[16],
                               size_t slot_j, uint64_t purpose) {
  static const uint8_t domain[] = "threshold-poprf-noise-v1";
  shake128_state_t st;

  if(NIZK2_DOUBLE_COMMITMENT) {
    /* Replicated holders share the opening, but not their error coins. */
    nizk2_private_seed(seed);
    return;
  }

  shake128_init(st);
  shake128_absorb(st,domain,sizeof(domain));
  shake128_absorb(st,pp->nizk2_pp,sizeof(pp->nizk2_pp));
  leopard_shake_absorb_u64(st,purpose);
  leopard_shake_absorb_u64(st,slot_j);
  if(ar_seed != NULL)
    shake128_absorb(st,ar_seed,16);
  if(sk != NULL) {
    leopard_shake_absorb_u64(st,sk->component_mask);
    leopard_shake_absorb_coeffs(st,sk->r_coeffs,sk->r_rows);
    leopard_shake_absorb_coeffs(st,sk->e_coeffs,sk->e_rows);
    leopard_shake_absorb_coeffs(st,sk->ei_coeffs,sk->ei_rows);
    leopard_shake_absorb_coeffs(st,sk->k_coeffs,sk->k_rows);
  }
  shake128_squeeze(st,seed,32);
  shake128_clear(st);
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

static leopard_coeff leopard_center_mod_q(leopard_ucoeff x) {
  const leopard_ucoeff q = LEOPARD_MOD_Q;

  if(x >= q)
    x -= q;
  if(x > q / 2)
    return -(leopard_coeff)(q - x);
  return (leopard_coeff)x;
}

static leopard_coeff leopard_q_div_pow2_i64(unsigned int logdiv) {
  const unsigned __int128 q = (unsigned __int128)LEOPARD_MOD_Q;
  return (leopard_coeff)(q >> logdiv);
}

typedef union {
  __m512i v[LEOPARD_RING_DEG / 32];
  int16_t c[LEOPARD_RING_DEG];
} leopard_ar_vecn;

typedef struct {
  leopard_ar_vecn limbs[LEOPARD_AR_LIMBS];
} leopard_ar_polz;

static void leopard_ar_polz_bitunpack(leopard_ar_polz *r,
                                      const uint8_t buf[LEOPARD_RING_DEG * LEOPARD_Q_BYTES]) {
  int i, j, k, bits;
  __m512i f, g, h, mask;

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<LEOPARD_AR_LIMBS;i++) {
    bits = (i + 1 < LEOPARD_AR_LIMBS) ? 14 : (LEOPARD_LOGQ - 14 * (LEOPARD_AR_LIMBS - 1));
    mask = _mm512_set1_epi16((int16_t)((1u << bits) - 1u));

    for(j=0;j<LEOPARD_RING_DEG / 32;j++) {
      if(k < bits) {
        f = _mm512_load_si512((const __m512i *)buf);
        g = _mm512_slli_epi16(f,k);
        h = _mm512_add_epi16(g,h);
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(f,bits-k);
        buf += 64;
      }
      else {
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(h,bits);
      }
      _mm512_store_si512(&r->limbs[i].v[j],g);
      k = (k - bits) & 0xF;
    }
  }
}

static int leopard_request_ar_coeffs(leopard_coeff *out, size_t rows, size_t cols,
                                     const uint8_t ar_seed[16]) {
  const size_t total_polys = rows * cols;
  const size_t poly_bytes = LEOPARD_RING_DEG * LEOPARD_Q_BYTES;
  uint64_t nonce = 1600;
  size_t poly_off = 0;

  while(poly_off < total_polys) {
    const size_t chunk_polys = (total_polys - poly_off > 32) ? 32 : (total_polys - poly_off);
    const size_t chunk_bytes = chunk_polys * poly_bytes;
    const size_t nblocks =
      (chunk_bytes + AES128CTR_BLOCKBYTES - 1) / AES128CTR_BLOCKBYTES;
    uint8_t *buf = checked_malloc(nblocks * AES128CTR_BLOCKBYTES,"request Ar bytes");
    aes128ctr_ctx aesctx;
    size_t local_idx, j;

    aes128ctr_init(&aesctx,ar_seed,nonce);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);

    for(local_idx=0;local_idx<chunk_polys;local_idx++) {
      const size_t src_poly = poly_off + local_idx;
      const size_t row = src_poly % rows;
      const size_t col = src_poly / rows;
      const size_t dst_poly = row * cols + col;
      const uint8_t *src = &buf[local_idx * poly_bytes];
      leopard_coeff *dst = &out[dst_poly * LEOPARD_RING_DEG];
      leopard_ar_polz unpacked;
      __attribute__((aligned(64)))
      uint8_t packed_poly[LEOPARD_RING_DEG * LEOPARD_Q_BYTES];

      memcpy(packed_poly,src,sizeof(packed_poly));
      leopard_ar_polz_bitunpack(&unpacked,packed_poly);
      for(j=0;j<LEOPARD_RING_DEG;j++) {
        leopard_ucoeff x = 0;
        size_t limb;

        for(limb=0;limb<LEOPARD_AR_LIMBS;limb++)
          x |= ((leopard_ucoeff)(uint16_t)unpacked.limbs[limb].c[j]) << (14 * limb);

        dst[j] = leopard_center_mod_q(x);
      }
    }

    free(buf);
    poly_off += chunk_polys;
    nonce += (uint64_t)1 << 32;
  }
  return 0;
}

static void raw_to_root_polyvec(polyvec_t dst, const leopard_coeff *src) {
#if LEOPARD_LOGQ > 64
  for(size_t i=0;i<polyvec_get_nelems(dst);i++) {
    poly_ptr p=polyvec_get_elem(dst,i);
    for(size_t j=0;j<LEOPARD_RING_DEG;j++) {
      int_ptr c=intvec_get_elem(p->coeffs,j);
      const leopard_coeff x=src[i*LEOPARD_RING_DEG+j];
      leopard_ucoeff mag=x<0 ? -x : x;
      int_set_zero(c);
      c->neg=x<0;
      c->limbs[0]=(uint64_t)mag;
      if(c->nlimbs<2) abort();
      c->limbs[1]=(uint64_t)(mag>>64);
    }
    p->crt=0;
  }
#else
  polyvec_set_coeffvec_i64(dst,src);
#endif
}

static void root_to_raw_polyvec(leopard_coeff *dst, polyvec_t src) {
#if LEOPARD_LOGQ > 64
  polyvec_redc(src,src);
  for(size_t i=0;i<polyvec_get_nelems(src);i++) {
    poly_ptr p=polyvec_get_elem(src,i);
    _fromcrt(p);
    for(size_t j=0;j<LEOPARD_RING_DEG;j++) {
      int_srcptr c=intvec_get_elem_src(p->coeffs,j);
      leopard_ucoeff mag=c->limbs[0];
      if(c->nlimbs>1) mag|=(leopard_ucoeff)c->limbs[1]<<64;
      for(size_t limb=2;limb<c->nlimbs;limb++) if(c->limbs[limb]) abort();
      dst[i*LEOPARD_RING_DEG+j]=c->neg ? -(leopard_coeff)mag : (leopard_coeff)mag;
    }
  }
#else
  polyvec_get_coeffvec_i64(dst,src);
#endif
}

static void raw_rows_to_root_polymat(polymat_t dst, const leopard_coeff *src,
                                     size_t rows, size_t cols) {
  size_t row;

  for(row=0;row<rows;row++) {
    polyvec_t rowv;
    polymat_get_row(rowv,dst,(unsigned int)row);
    raw_to_root_polyvec(rowv,&src[row * cols * LEOPARD_RING_DEG]);
  }
}

static void leopard_server_matrix_a1(polymat_t A1, const leopard_pp *pp) {
#if NIZK2_SMALL_COMMITMENT
  sample_polymat_bnd(A1,-NIZK2_COM_Q/2,NIZK2_COM_Q/2,pp->nizk2_pp,200000);
#elif NIZK2_LEGACY_I64_INSTANCE
  const leopard_coeff qhalf = leopard_q_div_pow2_i64(1);

  sample_polymat_bnd(A1,-qhalf,qhalf,pp->nizk2_pp,0);
#else
  sample_polymat_centered_mod(A1,pp->nizk2_pp,0);
#endif
}

static void leopard_server_matrix_b(polymat_t B, const leopard_pp *pp) {
#if NIZK2_SMALL_COMMITMENT
  sample_polymat_bnd(B,-NIZK2_COM_Q/2,NIZK2_COM_Q/2,pp->nizk2_pp,210000);
#elif NIZK2_LEGACY_I64_INSTANCE
  const leopard_coeff qhalf = leopard_q_div_pow2_i64(1);

  sample_polymat_bnd(B,-qhalf,qhalf,pp->nizk2_pp,100);
#else
  sample_polymat_centered_mod(B,pp->nizk2_pp,100);
#endif
}

static uint64_t leopard_double_scale(void) {
  const unsigned __int128 q=(unsigned __int128)LEOPARD_MOD_Q;
  uint64_t lo=0, hi=(uint64_t)1<<((LEOPARD_LOGQ+1)/2);
  while(lo+1<hi) {
    uint64_t mid=lo+(hi-lo)/2;
    if((unsigned __int128)mid*mid<=q) lo=mid; else hi=mid;
  }
  return lo;
}

static void leopard_commitment_encoding(polymat_t encoding) {
#if NIZK2_DOUBLE_COMMITMENT
  polymat_set_zero(encoding);
  for(unsigned int i=0;i<N;i++) {
    set_const_poly(polymat_get_elem(encoding,i,i),1);
    set_const_poly(polymat_get_elem(encoding,N+i,i),(leopard_coeff)leopard_double_scale());
  }
#elif NIZK2_SMALL_COMMITMENT
  set_const_diag(encoding,NIZK2_COM_SCALE);
#elif NIZK2_LEGACY_I64_INSTANCE
  set_const_diag(encoding,leopard_q_div_pow2_i64(4));
#else
  set_const_diag_q_div_pow2(encoding,4);
#endif
}

static void leopard_uniform_key(polyvec_t k, const uint8_t seed[32], uint32_t dom) {
  int_srcptr q=polyring_get_mod(polyvec_get_ring(k));
  INT_T(lo,q->nlimbs);
  INT_T(hi,q->nlimbs);
  int_rshift(hi,q,1);
  int_neg(lo,hi);
  polyvec_urandom_bnd(k,lo,hi,seed,dom);
}

static int leopard_keygen_with_nonce(leopard_pk_data *out_pk,
                                     leopard_sk_data *out_sk,
                                     const leopard_pp *pp,
                                     const leopard_coeff *given_k_coeffs,
                                     uint64_t component_mask,
                                     uint64_t nonce_base) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  polymat_t A1, B, q16I;
  polyvec_t r, e, ei, k, c, cu, tmp;
#if NIZK2_SMALL_COMMITMENT || NIZK2_DOUBLE_COMMITMENT
  uint8_t private_seed[32];
  nizk2_private_seed(private_seed);
  const uint8_t *sampling_seed = private_seed;
#else
  const uint8_t *sampling_seed = pp->nizk2_pp;
#endif

  memset(out_pk,0,sizeof(*out_pk));
  memset(out_sk,0,sizeof(*out_sk));
  out_pk->component_mask = component_mask;
  out_pk->ck_rows = COM_ROWS;
  out_pk->ck_coeffs = alloc_coeffs(out_pk->ck_rows);
  out_sk->component_mask = component_mask;
  out_sk->r_rows = COM_RANK;
  out_sk->e_rows = COM_RANK;
  out_sk->ei_rows = COM_MESSAGE_ROWS;
  out_sk->k_rows = N;
  out_sk->r_coeffs = alloc_coeffs(out_sk->r_rows);
  out_sk->e_coeffs = alloc_coeffs(out_sk->e_rows);
  out_sk->ei_coeffs = alloc_coeffs(out_sk->ei_rows);
  out_sk->k_coeffs = alloc_coeffs(out_sk->k_rows);

  polymat_alloc(A1,Rq,COM_RANK,COM_RANK);
  polymat_alloc(B,Rq,COM_MESSAGE_ROWS,COM_RANK);
  polymat_alloc(q16I,Rq,COM_MESSAGE_ROWS,N);
  polyvec_alloc(r,Rq,COM_RANK);
  polyvec_alloc(e,Rq,COM_RANK);
  polyvec_alloc(ei,Rq,COM_MESSAGE_ROWS);
  polyvec_alloc(k,Rq,N);
  polyvec_alloc(c,Rq,COM_RANK);
  polyvec_alloc(cu,Rq,COM_MESSAGE_ROWS);
  polyvec_alloc(tmp,Rq,COM_MESSAGE_ROWS);

  leopard_server_matrix_a1(A1,pp);
  leopard_server_matrix_b(B,pp);
  leopard_commitment_encoding(q16I);

  sample_polyvec_bnd(r,0,1,sampling_seed,nonce_base + 0);
  sample_polyvec_bnd(e,0,1,sampling_seed,nonce_base + 1);
  sample_polyvec_bnd(ei,0,1,sampling_seed,nonce_base + 2);
  if(given_k_coeffs != NULL)
    raw_to_root_polyvec(k,given_k_coeffs);
  else if(NIZK2_DOUBLE_COMMITMENT)
    leopard_uniform_key(k,sampling_seed,(uint32_t)(nonce_base + 3));
  else
    sample_polyvec_bnd(k,0,15,sampling_seed,nonce_base + 3);

  polyvec_mul(c,A1,r);
  polyvec_add(c,c,e,0);
  polyvec_mul(cu,B,r);
  polyvec_mul(tmp,q16I,k);
  polyvec_add(cu,cu,tmp,0);
  polyvec_add(cu,cu,ei,0);

  root_to_raw_polyvec(out_sk->r_coeffs,r);
  root_to_raw_polyvec(out_sk->e_coeffs,e);
  root_to_raw_polyvec(out_sk->ei_coeffs,ei);
  root_to_raw_polyvec(out_sk->k_coeffs,k);
  root_to_raw_polyvec(out_pk->ck_coeffs,c);
  root_to_raw_polyvec(&out_pk->ck_coeffs[COM_RANK * LEOPARD_RING_DEG],cu);
#if NIZK2_SMALL_COMMITMENT
  for(size_t i=0;i<out_pk->ck_rows*LEOPARD_RING_DEG;i++) {
    leopard_coeff v=out_pk->ck_coeffs[i] % NIZK2_COM_Q;
    if(v>NIZK2_COM_Q/2) v-=NIZK2_COM_Q;
    if(v< -NIZK2_COM_Q/2) v+=NIZK2_COM_Q;
    out_pk->ck_coeffs[i]=v;
  }
#endif

  polyvec_free(r);
  polyvec_free(e);
  polyvec_free(ei);
  polyvec_free(k);
  polyvec_free(c);
  polyvec_free(cu);
  polyvec_free(tmp);
  polymat_free(A1);
  polymat_free(B);
  polymat_free(q16I);
  return 0;
}

static int leopard_compute_vk_from_ar_seed(leopard_coeff *out_vk_coeffs,
                                           leopard_coeff *out_es_coeffs,
                                           const leopard_pp *pp,
                                           const leopard_sk_data *sk,
                                           const uint8_t ar_seed[16],
                                           size_t slot_j) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  polymat_t Ar;
  polyvec_t k, es, vk;
  leopard_coeff *ar_coeffs = NULL;
  uint8_t noise_seed[32];

  if(sk->k_rows != (size_t)N)
    return 1;

  ar_coeffs = alloc_coeffs((M + N) * N);
  leopard_request_ar_coeffs(ar_coeffs,M + N,N,ar_seed);

  polymat_alloc(Ar,Rq,M + N,N);
  polyvec_alloc(k,Rq,N);
  polyvec_alloc(es,Rq,M + N);
  polyvec_alloc(vk,Rq,M + N);
  raw_rows_to_root_polymat(Ar,ar_coeffs,M + N,N);
  raw_to_root_polyvec(k,sk->k_coeffs);
  leopard_noise_seed(noise_seed,pp,sk,ar_seed,slot_j,1);
#ifdef NIZK2_GAUSSIAN_ERRORS
  polyvec_grandom(es,NIZK2_ES_LOG2STDEV,noise_seed,0);
#else
  sample_polyvec_bnd(es,0,14,noise_seed,0);
#endif
  polyvec_mul(vk,Ar,k);
  polyvec_add(vk,vk,es,0);
  root_to_raw_polyvec(out_vk_coeffs,vk);
  if(out_es_coeffs != NULL)
    root_to_raw_polyvec(out_es_coeffs,es);

  free(ar_coeffs);
  polyvec_free(vk);
  polyvec_free(es);
  polyvec_free(k);
  polymat_free(Ar);
  return 0;
}

int leopard_nizk2_keygen_data(leopard_pk_data *out_pk,
                              leopard_sk_data *out_sk,
                              const leopard_pp *pp) {
  return leopard_keygen_with_nonce(out_pk,out_sk,pp,NULL,0,1000);
}

int leopard_nizk2_keygen_given_key_data(leopard_pk_data *out_pk,
                                        leopard_sk_data *out_sk,
                                        const leopard_pp *pp,
                                        const leopard_coeff *k_coeffs,
                                        uint64_t component_mask,
                                        uint64_t nonce_base) {
  return leopard_keygen_with_nonce(out_pk,out_sk,pp,k_coeffs,component_mask,
                                   nonce_base);
}

size_t leopard_nizk2_pk_data_bytes(const leopard_pk_data *pk) {
#if NIZK2_SMALL_COMMITMENT
  return 28 + (pk->ck_rows * LEOPARD_RING_DEG * NIZK2_COM_BITS + 7)/8;
#else
  return pk->ck_rows * LEOPARD_RING_DEG * LEOPARD_Q_BYTES;
#endif
}

static int leopard_build_reply_instance(polymat_t A, polyvec_t s, polyvec_t t,
                                        polyvec_t vk_out, polyvec_t ux_out,
                                        const leopard_pp *pp,
                                        const leopard_pk_data *pk,
                                        const leopard_sk_data *sk,
                                        const leopard_request_data *req,
                                        const leopard_coeff *pre_vk_coeffs,
                                        const leopard_coeff *pre_es_coeffs,
                                        double *protocol_eval_sec) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  polymat_t A1, Id, B, q16I, Ar, Cx, GI;
  polyvec_t s_r, s_e, s_ei, s_k, s_es, s_esp, es_prime_raw, tmp_h;
  polyvec_t t_c, t_cu, t_vk, t_uk;
  leopard_coeff *ar_coeffs = NULL;
  uint8_t noise_seed[32];
  struct timespec t0, t1;

  if(req->cx_rows != (size_t)H || req->cx_cols != (size_t)N ||
     pk->ck_rows != (size_t)COM_ROWS ||
     (sk != NULL && (sk->r_rows != (size_t)COM_RANK || sk->e_rows != (size_t)COM_RANK ||
                     sk->ei_rows != (size_t)COM_MESSAGE_ROWS || sk->k_rows != (size_t)N))) {
    return 1;
  }

  polymat_set_zero(A);
  polyvec_set_zero(s);
  polyvec_set_zero(t);

  polyvec_get_subvec(s_r,s,0,COM_RANK,1);
  polyvec_get_subvec(s_e,s,COM_RANK,COM_RANK,1);
  polyvec_get_subvec(s_ei,s,2 * COM_RANK,COM_MESSAGE_ROWS,1);
  polyvec_get_subvec(s_k,s,WIT_KEY_OFF,N,1);
  polyvec_get_subvec(s_es,s,WIT_ES_OFF,M + N,1);
  polyvec_get_subvec(s_esp,s,WIT_DIGIT_OFF,H * DECOMP_LEN,1);

  polymat_get_submat(A1,A,0,0,COM_RANK,COM_RANK,1,1);
  polymat_get_submat(Id,A,0,COM_RANK,COM_RANK,COM_RANK,1,1);
  leopard_server_matrix_a1(A1,pp);
  polymat_set_one(Id);

  polymat_get_submat(B,A,COM_RANK,0,COM_MESSAGE_ROWS,COM_RANK,1,1);
  polymat_get_submat(q16I,A,COM_RANK,WIT_KEY_OFF,COM_MESSAGE_ROWS,N,1,1);
  polymat_get_submat(Id,A,COM_RANK,2 * COM_RANK,COM_MESSAGE_ROWS,COM_MESSAGE_ROWS,1,1);
  leopard_server_matrix_b(B,pp);
  leopard_commitment_encoding(q16I);
  polymat_set_one(Id);

  polymat_get_submat(Ar,A,COM_ROWS,WIT_KEY_OFF,M + N,N,1,1);
  polymat_get_submat(Id,A,COM_ROWS,WIT_ES_OFF,M + N,M + N,1,1);
  ar_coeffs = alloc_coeffs((M + N) * N);
  leopard_request_ar_coeffs(ar_coeffs,M + N,N,req->ar_seed);
  raw_rows_to_root_polymat(Ar,ar_coeffs,M + N,N);
  polymat_set_one(Id);

  polymat_get_submat(Cx,A,COM_ROWS + M + N,WIT_KEY_OFF,H,N,1,1);
  polymat_get_submat(GI,A,COM_ROWS + M + N,WIT_DIGIT_OFF,H,H * DECOMP_LEN,1,1);
  raw_rows_to_root_polymat(Cx,req->cx_coeffs,req->cx_rows,req->cx_cols);
  set_gi(GI);

  if(protocol_eval_sec != NULL)
    *protocol_eval_sec = 0.0;
  if(sk != NULL) {
    raw_to_root_polyvec(s_r,sk->r_coeffs);
    raw_to_root_polyvec(s_e,sk->e_coeffs);
    raw_to_root_polyvec(s_ei,sk->ei_coeffs);
    raw_to_root_polyvec(s_k,sk->k_coeffs);
    if(pre_es_coeffs != NULL)
      raw_to_root_polyvec(s_es,pre_es_coeffs);
    else {
      leopard_noise_seed(noise_seed,pp,sk,req->ar_seed,req->slot_j,1);
#ifdef NIZK2_GAUSSIAN_ERRORS
      polyvec_grandom(s_es,NIZK2_ES_LOG2STDEV,noise_seed,0);
#else
      sample_polyvec_bnd(s_es,0,14,noise_seed,0);
#endif
    }
    if(pre_vk_coeffs != NULL)
      raw_to_root_polyvec(vk_out,pre_vk_coeffs);
    else {
      polyvec_mul(vk_out,Ar,s_k);
      polyvec_add(vk_out,vk_out,s_es,0);
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    polyvec_alloc(es_prime_raw,Rq,H);
    polyvec_alloc(tmp_h,Rq,H);
    leopard_noise_seed(noise_seed,pp,sk,req->ar_seed,req->slot_j,2);
#ifdef NIZK2_GAUSSIAN_ERRORS
    polyvec_grandom(es_prime_raw,NIZK2_ESP_LOG2STDEV,noise_seed,0);
#else
    sample_polyvec_bnd(es_prime_raw,0,ES_PRIME_MAX,noise_seed,0);
#endif
    decompose_es_prime(s_esp,es_prime_raw);
    polyvec_mul(ux_out,Cx,s_k);
    polyvec_mul(tmp_h,GI,s_esp);
    polyvec_add(ux_out,ux_out,tmp_h,0);
    polyvec_free(tmp_h);
    polyvec_free(es_prime_raw);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    if(protocol_eval_sec != NULL)
      *protocol_eval_sec = leopard_elapsed_sec(&t0,&t1);
  }

  polyvec_get_subvec(t_c,t,0,COM_RANK,1);
  polyvec_get_subvec(t_cu,t,COM_RANK,COM_MESSAGE_ROWS,1);
  polyvec_get_subvec(t_vk,t,COM_ROWS,M + N,1);
  polyvec_get_subvec(t_uk,t,COM_ROWS + M + N,H,1);
  raw_to_root_polyvec(t_c,pk->ck_coeffs);
  raw_to_root_polyvec(t_cu,&pk->ck_coeffs[COM_RANK * LEOPARD_RING_DEG]);
  polyvec_set(t_vk,vk_out);
  polyvec_set(t_uk,ux_out);
  polyvec_neg_self(t);
  polyvec_mod(t,t);
  polyvec_redc(t,t);
  free(ar_coeffs);
  (void)Rq;
  return 0;
}

int leopard_nizk2_blindeval_data(leopard_reply_data *out_rep,
                                 const leopard_pp *pp,
                                 const leopard_pk_data *pk,
                                 const leopard_sk_data *sk,
                                 const uint8_t *tag, size_t tag_len,
                                 const leopard_request_data *req) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  const unsigned int nrels = REL_ROWS;
  const unsigned int nwit = WIT_ROWS;
  polymat_t A;
  polyvec_t s, t, vk, ux;
  struct timespec t0, t1;

  memset(out_rep,0,sizeof(*out_rep));
  memset(&g_last_blindeval_profile,0,sizeof(g_last_blindeval_profile));
  (void)tag;
  (void)tag_len;
  out_rep->component_mask = pk->component_mask;
  out_rep->ux_rows = H;
  out_rep->proof_len = lin_params_get_prooflen(params);
  out_rep->proof = checked_malloc(out_rep->proof_len,"reply proof");
  out_rep->ux_coeffs = alloc_coeffs(out_rep->ux_rows);

  polymat_alloc(A,Rq,nrels,nwit);
  polyvec_alloc(s,Rq,nwit);
  polyvec_alloc(t,Rq,nrels);
  polyvec_alloc(vk,Rq,M + N);
  polyvec_alloc(ux,Rq,H);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_build_reply_instance(A,s,t,vk,ux,pp,pk,sk,req,NULL,NULL,NULL) != 0) {
    polymat_free(A);
    polyvec_free(s);
    polyvec_free(t);
    polyvec_free(vk);
    polyvec_free(ux);
    return 1;
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.build_instance_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  out_rep->proof_len = prover(out_rep->proof,s,A,t,pp->nizk2_pp);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.prover_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  root_to_raw_polyvec(out_rep->ux_coeffs,ux);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.export_sec = leopard_elapsed_sec(&t0,&t1);

  polymat_free(A);
  polyvec_free(s);
  polyvec_free(t);
  polyvec_free(vk);
  polyvec_free(ux);
  return 0;
}

int leopard_nizk2_preproc_server_data(leopard_prep_data *out_prep,
                                      leopard_server_pre_state_data *out_st,
                                      const leopard_pp *pp,
                                      const leopard_sk_data *sk,
                                      const leopard_preq_data *preq) {
  size_t j;

  if(preq == NULL || preq->slot_count == 0)
    return 1;

  memset(out_prep,0,sizeof(*out_prep));
  memset(out_st,0,sizeof(*out_st));
  out_prep->slot_count = preq->slot_count;
  out_prep->component_mask = sk->component_mask;
  out_prep->vk_rows = M + N;
  out_prep->vk_coeffs = alloc_coeffs(out_prep->slot_count * out_prep->vk_rows);

  out_st->slot_count = preq->slot_count;
  out_st->component_mask = sk->component_mask;
  out_st->vk_rows = M + N;
  out_st->es_rows = M + N;
  out_st->vk_coeffs = alloc_coeffs(out_st->slot_count * out_st->vk_rows);
  out_st->es_coeffs = alloc_coeffs(out_st->slot_count * out_st->es_rows);

  for(j=0;j<preq->slot_count;j++) {
    leopard_coeff *vk_slot = &out_st->vk_coeffs[j * out_st->vk_rows * LEOPARD_RING_DEG];
    leopard_coeff *es_slot = &out_st->es_coeffs[j * out_st->es_rows * LEOPARD_RING_DEG];

    if(leopard_compute_vk_from_ar_seed(vk_slot,es_slot,pp,sk,
                                       &preq->ar_seeds[16 * j],j) != 0) {
      leopard_prep_data_clear(out_prep);
      leopard_server_pre_state_data_clear(out_st);
      return 1;
    }
    memcpy(&out_prep->vk_coeffs[j * out_prep->vk_rows * LEOPARD_RING_DEG],
           vk_slot,
           out_prep->vk_rows * LEOPARD_RING_DEG * sizeof(*vk_slot));
  }
  return 0;
}

int leopard_nizk2_blindeval_preprocessed_data(
    leopard_reply_data *out_rep,
    const leopard_pp *pp,
    const leopard_pk_data *pk,
    const leopard_sk_data *sk,
    const leopard_server_pre_state_data *pre_st,
    const uint8_t *tag, size_t tag_len,
    const leopard_request_data *req) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  const unsigned int nrels = REL_ROWS;
  const unsigned int nwit = WIT_ROWS;
  polymat_t A;
  polyvec_t s, t, vk, ux;
  const leopard_coeff *pre_vk_coeffs;
  const leopard_coeff *pre_es_coeffs;
  struct timespec t0, t1;

  if(pre_st == NULL || req->slot_j >= pre_st->slot_count ||
     pre_st->vk_rows != (size_t)(M + N) ||
     pre_st->es_rows != (size_t)(M + N)) {
    return 1;
  }

  memset(out_rep,0,sizeof(*out_rep));
  memset(&g_last_blindeval_profile,0,sizeof(g_last_blindeval_profile));
  (void)tag;
  (void)tag_len;
  out_rep->component_mask = pk->component_mask;
  out_rep->ux_rows = H;
  out_rep->proof_len = lin_params_get_prooflen(params);
  out_rep->proof = checked_malloc(out_rep->proof_len,"reply proof");
  out_rep->ux_coeffs = alloc_coeffs(out_rep->ux_rows);
  pre_vk_coeffs = &pre_st->vk_coeffs[req->slot_j * pre_st->vk_rows * LEOPARD_RING_DEG];
  pre_es_coeffs = &pre_st->es_coeffs[req->slot_j * pre_st->es_rows * LEOPARD_RING_DEG];

  polymat_alloc(A,Rq,nrels,nwit);
  polyvec_alloc(s,Rq,nwit);
  polyvec_alloc(t,Rq,nrels);
  polyvec_alloc(vk,Rq,M + N);
  polyvec_alloc(ux,Rq,H);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_build_reply_instance(A,s,t,vk,ux,pp,pk,sk,req,pre_vk_coeffs,pre_es_coeffs,
                                  &g_last_blindeval_profile.protocol_eval_sec) != 0) {
    polymat_free(A);
    polyvec_free(s);
    polyvec_free(t);
    polyvec_free(vk);
    polyvec_free(ux);
    return 1;
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.build_instance_sec = leopard_elapsed_sec(&t0,&t1);
  g_last_blindeval_profile.nizk_build_sec =
      g_last_blindeval_profile.build_instance_sec - g_last_blindeval_profile.protocol_eval_sec;
  if(g_last_blindeval_profile.nizk_build_sec < 0.0)
    g_last_blindeval_profile.nizk_build_sec = 0.0;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  out_rep->proof_len = prover(out_rep->proof,s,A,t,pp->nizk2_pp);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.prover_sec = leopard_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
  root_to_raw_polyvec(out_rep->ux_coeffs,ux);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_blindeval_profile.export_sec = leopard_elapsed_sec(&t0,&t1);

  polymat_free(A);
  polyvec_free(s);
  polyvec_free(t);
  polyvec_free(vk);
  polyvec_free(ux);
  return 0;
}

int leopard_nizk2_reply_verify_data(const leopard_pp *pp,
                                    const leopard_pk_data *pk,
                                    const leopard_prep_data *prep,
                                    const uint8_t *tag, size_t tag_len,
                                    const leopard_request_data *req,
                                    const leopard_reply_data *rep) {
  polyring_srcptr Rq = params->tbox_params->tbox->ring;
  const unsigned int nrels = REL_ROWS;
  const unsigned int nwit = WIT_ROWS;
  polymat_t A;
  polyvec_t t, vk, ux, dummy_s;
  struct timespec total0, total1, t0, t1;
  int accept = 0;
  int ret = 1;

  memset(&g_last_verify_profile,0,sizeof(g_last_verify_profile));
  (void)tag;
  (void)tag_len;
  if(rep->ux_rows != (size_t)H || prep == NULL ||
     prep->slot_count <= req->slot_j || prep->vk_rows != (size_t)(M + N))
    return 1;

  clock_gettime(CLOCK_MONOTONIC,&total0);
  polymat_alloc(A,Rq,nrels,nwit);
  polyvec_alloc(dummy_s,Rq,nwit);
  polyvec_alloc(t,Rq,nrels);
  polyvec_alloc(vk,Rq,M + N);
  polyvec_alloc(ux,Rq,H);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  raw_to_root_polyvec(vk,&prep->vk_coeffs[req->slot_j * prep->vk_rows * LEOPARD_RING_DEG]);
  raw_to_root_polyvec(ux,rep->ux_coeffs);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.import_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  if(leopard_build_reply_instance(A,dummy_s,t,vk,ux,pp,pk,NULL,req,NULL,NULL,NULL) != 0) {
    goto end;
  }
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.build_instance_sec = leopard_elapsed_sec(&t0,&t1);

  clock_gettime(CLOCK_MONOTONIC,&t0);
  accept = verifier(rep->proof,A,t,pp->nizk2_pp);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  g_last_verify_profile.verifier_sec = leopard_elapsed_sec(&t0,&t1);
  ret = accept ? 0 : 1;

end:
  polymat_free(A);
  polyvec_free(dummy_s);
  polyvec_free(t);
  polyvec_free(vk);
  polyvec_free(ux);
  clock_gettime(CLOCK_MONOTONIC,&total1);
  g_last_verify_profile.total_sec = leopard_elapsed_sec(&total0,&total1);
  return ret;
}

size_t leopard_nizk2_reply_data_bytes(const leopard_reply_data *rep) {
  return sizeof(uint64_t) +
         rep->ux_rows * LEOPARD_RING_DEG * LEOPARD_Q_BYTES + rep->proof_len;
}

const leopard_nizk2_blindeval_profile *leopard_nizk2_last_blindeval_profile(void) {
  return &g_last_blindeval_profile;
}

const leopard_nizk2_verify_profile *leopard_nizk2_last_verify_profile(void) {
  return &g_last_verify_profile;
}
