#include <stdint.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "data.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"

#ifndef POLX_FROMPOLY_PARALLEL_MIN_LEN
#define POLX_FROMPOLY_PARALLEL_MIN_LEN 2048
#endif
#ifndef POLX_FROMPOLY_MAX_THREADS
#define POLX_FROMPOLY_MAX_THREADS 8
#endif
#ifndef POLX_EXTENSION_PARALLEL_MIN_LEN
#define POLX_EXTENSION_PARALLEL_MIN_LEN 1024
#endif
#ifndef POLX_FROMINT64_PARALLEL_MIN_LEN
#define POLX_FROMINT64_PARALLEL_MIN_LEN 1024
#endif
#ifndef POLX_FROMINT64_MAX_THREADS
#define POLX_FROMINT64_MAX_THREADS 8
#endif
#ifndef POLX_EXTENSION_MAX_THREADS
#define POLX_EXTENSION_MAX_THREADS 8
#endif
#ifndef POLX_POINTWISE_BLOCK
#define POLX_POINTWISE_BLOCK 64
#endif
#ifndef POLX_FROMPOLY_BLOCK
#define POLX_FROMPOLY_BLOCK 64
#endif
#ifndef POLX_SPROD_BLOCK
#define POLX_SPROD_BLOCK 64
#endif
#if POLX_SPROD_BLOCK < 0 || POLX_SPROD_BLOCK % 2
#error "POLX_SPROD_BLOCK must be zero or a positive even number"
#endif

static labrador_scalar cmodq(labrador_scalar a) {
  const __int128 q = ((__int128)1 << LOGQ) - QOFF;
  const __int128 half = q/2;
  __int128 x = a;

  x %= q;
  if(x > half)
    x -= q;
  if(x < -half)
    x += q;
  return (labrador_scalar)x;
}

static int16_t scalar_to_mont_prime(labrador_scalar s, const pdata *prime) {
  __int128 t = (__int128)cmodq(s);

  t %= prime->p;
  t *= ((__int128)1 << 16) % prime->p;
  t %= prime->p;
  return (int16_t)t;
}

static int16_t scalar_to_prime_times(labrador_scalar s, int16_t factor, const pdata *prime) {
  __int128 t = (__int128)cmodq(s);

  t %= prime->p;
  t *= factor;
  t %= prime->p;
  return (int16_t)t;
}

void polx_print(const polx *a) {
  polz t[1];
  polz_frompolx(t,a);
  polz_center(t);
  polz_print(t);
}

void polx_getcoeff(zz *r, const polx *a, int k) {
  polz z;

  polz_frompolx(&z,a);
  polz_getcoeff(r,&z,k);
}

void polxvec_setzero(polx *r, size_t len) {
  polyvec_setzero(&r->vec[0],len*K);
}

void polxvec_copy(polx *r, const polx *a, size_t len) {
  polyvec_copy(&r->vec[0],&a->vec[0],len*K);
}

void polxvec_fromint64vec(polx *r, size_t len, size_t deg, const int64_t a[len*deg*N]) {
  size_t i;

#if defined(_OPENMP)
  int threads = omp_get_max_threads();
  if(threads > POLX_FROMINT64_MAX_THREADS)
    threads = POLX_FROMINT64_MAX_THREADS;

  if(len >= POLX_FROMINT64_PARALLEL_MIN_LEN && threads > 1) {
#pragma omp parallel for schedule(static) num_threads(threads)
    for(i=0;i<len;i++) {
      polz t[deg];
      polzvec_fromint64vec(t,1,deg,&a[i*deg*N]);
      polzvec_topolxvec(&r[i*deg],t,deg);
    }
    return;
  }
#endif

  for(i=0;i<len;i++) {
    polz t[deg];
    polzvec_fromint64vec(t,1,deg,&a[i*deg*N]);
    polzvec_topolxvec(&r[i*deg],t,deg);
  }
}

int polx_iszero(const polx *a) {
  polz t;

  polz_frompolx(&t,a);
  polz_center(&t);
  return polz_iszero(&t);
}

int polxvec_iszero(const polx *a, size_t len) {
  polz t[len];

  polzvec_frompolxvec(t,a,len);
  polzvec_center(t,len);
  return polzvec_iszero(t,len);
}

int polx_iszero_constcoeff(const polx *a) {
  polz t;

  polz_frompolx(&t,a);
  polz_center(&t);
  return polz_iszero_constcoeff(&t);
}

void polx_monomial(polx *r, labrador_scalar v, int k) {
  size_t i;

  v = cmodq(v);
  for(i=0;i<K;i++)
    poly_monomial_ntt(&r->vec[i],v % primes[i].p,k,&primes[i]);
}

void polxvec_ternary(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[32];

  while(len >= 32) {
    polyvec_ternary(t,32,seed,nonce);
    polxvec_frompolyvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polyvec_ternary(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_quarternary(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[32];

  while(len >= 32) {
    polyvec_quarternary(t,32,seed,nonce);
    polxvec_frompolyvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polyvec_quarternary(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_challenge(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[10];

  while(len >= 10) {
    polyvec_challenge(t,10,seed,nonce);
    polxvec_frompolyvec(r,t,10);
    nonce += (uint64_t)1 << 32;
    r += 10;
    len -= 10;
  }

  polyvec_challenge(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_almostuniform(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  polz t[32];

  while(len >= 32) {
    polzvec_almostuniform(t,32,seed,nonce);
    polzvec_topolxvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polzvec_almostuniform(t,len,seed,nonce);
  polzvec_topolxvec(r,t,len);
}

void polx_frompoly(polx *r, const poly *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],a,primes[i].s,&primes[i]);

  polx_ntt(r,r);
}

void polxvec_frompolyvec(polx *r, const poly *a, size_t len) {
  size_t i;

#if POLX_FROMPOLY_BLOCK > 0
  if(len > POLX_FROMPOLY_BLOCK) {
#if defined(_OPENMP)
    const int threads = MIN(omp_get_max_threads(),POLX_FROMPOLY_MAX_THREADS);
#pragma omp parallel for schedule(static) num_threads(threads) if(len >= POLX_FROMPOLY_PARALLEL_MIN_LEN && threads > 1)
#endif
    for(size_t off=0;off<len;off+=POLX_FROMPOLY_BLOCK) {
      const size_t count = MIN(len-off,POLX_FROMPOLY_BLOCK);
      for(size_t limb=0;limb<K;limb++) {
        polyvec_scale_widening(&r[off].vec[limb],&a[off],count,K,primes[limb].s,&primes[limb]);
        polyvec_ntt(&r[off].vec[limb],&r[off].vec[limb],count,K,&primes[limb]);
      }
    }
    return;
  }
#endif
#if defined(_OPENMP)
  int threads = omp_get_max_threads();
  if(threads > POLX_FROMPOLY_MAX_THREADS)
    threads = POLX_FROMPOLY_MAX_THREADS;
  if(threads > K)
    threads = K;

  if(len >= POLX_FROMPOLY_PARALLEL_MIN_LEN && threads > 1) {
#pragma omp parallel num_threads(threads)
    {
#pragma omp for schedule(static)
      for(i=0;i<K;i++)
        polyvec_scale_widening(&r->vec[i],a,len,K,primes[i].s,&primes[i]);

#pragma omp for schedule(static)
      for(i=0;i<K;i++)
        polyvec_ntt(&r->vec[i],&r->vec[i],len,K,&primes[i]);
    }
    return;
  }
#endif

  for(i=0;i<K;i++)
    polyvec_scale_widening(&r->vec[i],a,len,K,primes[i].s,&primes[i]);

  polxvec_ntt(r,r,len);
}

void polx_refresh(polx *r) {
  polz a;

  polz_frompolx(&a,r);
  polz_center(&a);
  polz_topolx(r,&a);
}

void polxvec_refresh(polx *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polx_refresh(&r[i]);
}

void polx_reduce(polx *r) {
  size_t i;

  for(i=0;i<K;i++)
    poly_reduce(&r->vec[i],&primes[i]);
}

void polxvec_reduce(polx *r, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_reduce(&r->vec[i],len,K,&primes[i]);
}

void polx_neg(polx *r, const polx *a) {
  polyvec_neg(&r->vec[0],&a->vec[0],K);
}

void polxvec_neg(polx *r, const polx *a, size_t len) {
  polyvec_neg(&r->vec[0],&a->vec[0],K*len);
}

void polx_add(polx *r, const polx *a, const polx *b) {
  polyvec_add(&r->vec[0],&a->vec[0],&b->vec[0],K);
  polx_reduce(r);
}

void polxvec_add(polx *r, const polx *a, const polx *b, size_t len) {
  polyvec_add(&r->vec[0],&a->vec[0],&b->vec[0],K*len);
  polxvec_reduce(r,len);
}

void polx_sub(polx *r,const polx *a,const polx *b) {
  polyvec_sub(&r->vec[0],&a->vec[0],&b->vec[0],K);
  polx_reduce(r);
}

void polxvec_sub(polx *r, const polx *a, const polx *b, size_t len) {
  polyvec_sub(&r->vec[0],&a->vec[0],&b->vec[0],K*len);
  polxvec_reduce(r,len);
}

void polx_ntt(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_ntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_ntt(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_ntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

void polx_invntt(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_invntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_invntt(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_invntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

void polx_mul(polx *r, const polx *a, const polx *b) {
  size_t i;

  for(i=0;i<K;i++)
    poly_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],&primes[i]);
}

void polx_poly_mul(polx *r, const polx *a, const poly *b) {
  polx t;

  polx_frompoly(&t,b);
  polx_mul(r,a,&t);
}

void polxvec_mul(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polxvec_polx_mul(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

#if POLX_POINTWISE_BLOCK > 0
  if(len > POLX_POINTWISE_BLOCK) {
    const polx scalar = *a;
    /* Visit all CRT limbs while this coefficient block is still in cache. */
    for(size_t off=0;off<len;off+=POLX_POINTWISE_BLOCK) {
      const size_t count = MIN(len-off,POLX_POINTWISE_BLOCK);
      for(i=0;i<K;i++)
        polyvec_poly_pointwise(&r[off].vec[i],&scalar.vec[i],&b[off].vec[i],count,K,&primes[i]);
    }
    return;
  }
#endif
  for(i=0;i<K;i++)
    polyvec_poly_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polx_mul_add(polx *r, const polx *a, const polx *b) {
  size_t i;

  for(i=0;i<K;i++)
    poly_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],&primes[i]);

  polx_reduce(r);
}

void polxvec_mul_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);

  polxvec_reduce(r,len);
}

void polxvec_polx_mul_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

#if POLX_POINTWISE_BLOCK > 0
  if(len > POLX_POINTWISE_BLOCK) {
    const polx scalar = *a;
    for(size_t off=0;off<len;off+=POLX_POINTWISE_BLOCK) {
      const size_t count = MIN(len-off,POLX_POINTWISE_BLOCK);
      for(i=0;i<K;i++)
        polyvec_poly_pointwise_add(&r[off].vec[i],&scalar.vec[i],&b[off].vec[i],count,K,&primes[i]);
      polxvec_reduce(&r[off],count);
    }
    return;
  }
#endif
  for(i=0;i<K;i++)
    polyvec_poly_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);

  polxvec_reduce(r,len);
}

void polxvec_sprod(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  if(len == 0) {
    polxvec_setzero(r,1);
    return;
  }

#if POLX_SPROD_BLOCK > 0
  if(len > POLX_SPROD_BLOCK) {
    polx_mul(r,a,b);
    /* Even blocks preserve the original pairwise reduction schedule. */
    for(size_t off=1;off<len;off+=POLX_SPROD_BLOCK) {
      const size_t count = MIN(len-off,POLX_SPROD_BLOCK);
      for(i=0;i<K;i++)
        polyvec_sprod_pointwise_add(&r->vec[i],&a[off].vec[i],&b[off].vec[i],count,K,&primes[i]);
    }
    return;
  }
#endif
  for(i=0;i<K;i++)
    polyvec_sprod_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polxvec_sprod_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

#if POLX_SPROD_BLOCK > 0
  if(len > POLX_SPROD_BLOCK) {
    for(size_t off=0;off<len;off+=POLX_SPROD_BLOCK) {
      const size_t count = MIN(len-off,POLX_SPROD_BLOCK);
      for(i=0;i<K;i++)
        polyvec_sprod_pointwise_add(&r->vec[i],&a[off].vec[i],&b[off].vec[i],count,K,&primes[i]);
    }
    return;
  }
#endif
  for(i=0;i<K;i++)
    polyvec_sprod_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polx_scale(polx *r, const polx *a, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],&a->vec[i],scalar_to_mont_prime(s,&primes[i]),&primes[i]);
}

void polx_scale_frompoly(polx *r, const poly *a, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],a,scalar_to_prime_times(s,primes[i].s,&primes[i]),&primes[i]);

  polx_ntt(r,r);
}

void polxvec_scale(polx *r, const polx *a, size_t len, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_scale(&r->vec[i],&a->vec[i],len,K,scalar_to_mont_prime(s,&primes[i]),&primes[i]);
}

void polxvec_scale_frompolyvec(polx *r, const poly *a, size_t len, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_scale_widening(&r->vec[i],a,len,K,scalar_to_prime_times(s,primes[i].s,&primes[i]),&primes[i]);

  polxvec_ntt(r,r,len);
}

void polx_scale_add(polx *r, const polx *a, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    poly_scale_add(&r->vec[i],&a->vec[i],scalar_to_mont_prime(s,&primes[i]),&primes[i]);

  polx_reduce(r);
}

void polxvec_scale_add(polx *r, const polx *a, size_t len, labrador_scalar s) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_scale_add(&r->vec[i],&a->vec[i],len,K,scalar_to_mont_prime(s,&primes[i]),&primes[i]);

  polxvec_reduce(r,len);
}

size_t polxvec_mul_extension(polx *c, const polx *a, const polx *b, size_t len, size_t deg, size_t mult) {
  size_t i,j,k = 0;

#if defined(_OPENMP)
  int threads = omp_get_max_threads();
  const size_t tasks = mult*K;

  if(threads > POLX_EXTENSION_MAX_THREADS)
    threads = POLX_EXTENSION_MAX_THREADS;
  if((size_t)threads > tasks)
    threads = (int)tasks;

  if(len >= POLX_EXTENSION_PARALLEL_MIN_LEN && threads > 1) {
#pragma omp parallel for schedule(static) num_threads(threads) reduction(max:k)
    for(size_t task=0;task<tasks;task++) {
      const size_t ii = task / K;
      const size_t jj = task % K;
      const size_t kk = polyvec_pointwise_extension(&c[ii].vec[jj],&a->vec[jj],&b[ii].vec[jj],
                                                    len,K*mult,deg,&primes[jj]);
      if(kk > k)
        k = kk;
    }
    return k;
  }
#endif

  for(i=0;i<mult;i++)
    for(j=0;j<K;j++)
      k = polyvec_pointwise_extension(&c[i].vec[j],&a->vec[j],&b[i].vec[j],len,K*mult,deg,&primes[j]);

  return k;
}

size_t polxvec_collaps_add_extension(polx *c, const polx *a, const polx *b, size_t len, size_t deg, size_t mult) {
  size_t i,j,k = 0;

#if defined(_OPENMP)
  int threads = omp_get_max_threads();
  const size_t tasks = mult*K;

  if(threads > POLX_EXTENSION_MAX_THREADS)
    threads = POLX_EXTENSION_MAX_THREADS;
  if((size_t)threads > tasks)
    threads = (int)tasks;

  if(len >= POLX_EXTENSION_PARALLEL_MIN_LEN && threads > 1) {
#pragma omp parallel for schedule(static) num_threads(threads) reduction(max:k)
    for(size_t task=0;task<tasks;task++) {
      const size_t ii = task / K;
      const size_t jj = task % K;
      const size_t kk = polyvec_collaps_add_extension(&c[ii].vec[jj],&a[ii].vec[jj],&b->vec[jj],
                                                      len,K*mult,deg,&primes[jj]);
      if(kk > k)
        k = kk;
    }
    return k;
  }
#endif

  for(i=0;i<mult;i++)
    for(j=0;j<K;j++)
      k = polyvec_collaps_add_extension(&c[i].vec[j],&a[i].vec[j],&b->vec[j],len,K*mult,deg,&primes[j]);

  return k;
}

void polxvec_decompose(poly *r, const polx *a, size_t len, size_t t, size_t d) {
  size_t i;
  const size_t stride = len;
  polz b[16];

  while(len >= 16) {
    polzvec_frompolxvec(b,a,16);
    polzvec_center(b,16);
    for(i=0;i<16;i++)
      polz_decompose(&r[i],&b[i],stride,t,d);
    r += 16;
    a += 16;
    len -= 16;
  }

  polzvec_frompolxvec(b,a,len);
  polzvec_center(b,len);
  for(i=0;i<len;i++)
    polz_decompose(&r[i],&b[i],stride,t,d);
}

void polxvec_reconstruct(polx *r, const poly *a, size_t len, size_t t, size_t d) {
  size_t i;
  const size_t stride = len;
  polz b[16];

  while(len >= 16) {
    for(i=0;i<16;i++)
      polz_reconstruct(&b[i],&a[i],stride,t,d);
    polzvec_topolxvec(r,b,16);
    r += 16;
    a += 16;
    len -= 16;
  }

  for(i=0;i<len;i++)
    polz_reconstruct(&b[i],&a[i],stride,t,d);
  polzvec_topolxvec(r,b,len);
}

void polx_sigmam1(polx *r, const polx *a) {
  polyvec_sigmam1_ntt(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigmam1(polx *r, const polx *a, size_t len) {
  polyvec_sigmam1_ntt(&r->vec[0],&a->vec[0],K*len);
}

/* FIXME: NTT rep! */
void polx_sigma5(polx *r, const polx *a) {
  polyvec_sigma5(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigma5(polx *r, const polx *a, size_t len) {
  polyvec_sigma5(&r->vec[0],&a->vec[0],K*len);
}

void polx_sigma5inv(polx *r, const polx *a) {
  polyvec_sigma5inv(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigma5inv(polx *r, const polx *a, size_t len) {
  polyvec_sigma5inv(&r->vec[0],&a->vec[0],K*len);
}

void polx_flip(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_flip_ntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_flip(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_flip_ntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

