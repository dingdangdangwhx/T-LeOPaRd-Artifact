#ifndef LEOPARD_PUBLIC_FASTPATH_H
#define LEOPARD_PUBLIC_FASTPATH_H

/* Private helpers, included after the selector arithmetic implementation. */
typedef struct {
  double forward_sec;
  double product_sec;
  double decompose_sec;
} leopard_public_chain_profile;

static int leopard_public_threads(void) {
  int threads = 1;
#ifdef _OPENMP
  const char *value = getenv("LEOPARD_PUBLIC_THREADS");
  threads = omp_get_max_threads();
  if(threads > 8)
    threads = 8;
  if(value != NULL) {
    char *end;
    long requested = strtol(value,&end,10);
    if(end != value && *end == '\0' && requested > 0 && requested <= 1024)
      threads = requested < omp_get_max_threads() ? (int)requested : omp_get_max_threads();
  }
#endif
  return threads;
}

static void leopard_public_almostuniform(polx *out, size_t len,
    const uint8_t seed[16], uint64_t nonce, int threads) {
  if(threads <= 1 || len < 64) {
    polxvec_almostuniform(out,len,seed,nonce);
    return;
  }
  /* Match polxvec_almostuniform's 32-polynomial chunks and high nonce word. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for(size_t chunk=0;chunk<(len + 31)/32;chunk++) {
    polz buffer[32];
    const size_t offset = 32 * chunk;
    const size_t count = len - offset < 32 ? len - offset : 32;
    polzvec_almostuniform(buffer,count,seed,nonce + ((uint64_t)chunk << 32));
    polzvec_topolxvec(out + offset,buffer,count);
  }
}

static void leopard_public_decompose(poly *out, const polx *in, size_t digits) {
  polz z[1];
  const unsigned __int128 beta4 = (unsigned __int128)GADGET_BASE *
      GADGET_BASE * GADGET_BASE * GADGET_BASE;
  const unsigned __int128 corr =
      (beta4 * beta4 * (digits == 9 ? GADGET_BASE : 1) - 1) / 2;

  if((digits != 8 && digits != 9) || GADGET_BASE % 2 == 0) {
    gadget_decompose_balanced3(out,in,digits);
    return;
  }
  polz_frompolx(z,in);
  polz_caddq(z);
  for(size_t j=0;j<N;j++) {
    zz coeff;
    size_t i = 0;
    unsigned __int128 wide;
    uint64_t narrow;

    polz_getcoeff(&coeff,z,j);
    wide = (unsigned __int128)zz_to_u64(&coeff) + corr;
    /* Preserve the canonical digits; only the arithmetic width changes. */
    while(wide > UINT64_MAX && i < digits) {
      out[i++].vec->c[j] = (int16_t)(wide % GADGET_BASE) - GADGET_CENTER;
      wide /= GADGET_BASE;
    }
    narrow = (uint64_t)wide;
    for(;i<digits;i++) {
      out[i].vec->c[j] = (int16_t)(narrow % GADGET_BASE) - GADGET_CENTER;
      narrow /= GADGET_BASE;
    }
  }
}

static poly *leopard_public_chain(
    polx *bnextx, poly *buffer0, poly *buffer1,
    const poly *g0, const poly *g1, const polx *A0, const polx *A1,
    const uint8_t *xbits, int threads, int parallel_conversion,
    leopard_public_chain_profile *profile) {
  const size_t m = DEMO_M;
  const size_t block = gadget_len_base3();
  const size_t len = m * block;

  polyvec_copy(buffer1,xbits[DEMO_L - 1] ? g1 : g0,len);
  if(profile != NULL)
    memset(profile,0,sizeof(*profile));
  if(threads <= 1 || profile != NULL) {
    poly *next = buffer1, *cur = buffer0;
    struct timespec t0, t1;

    for(size_t t=DEMO_L - 1;t-- > 0;) {
      const polx *matrix = xbits[t] ? A1 : A0;
      if(profile != NULL)
        clock_gettime(CLOCK_MONOTONIC,&t0);
      polxvec_frompolyvec(bnextx,next,len);
      if(profile != NULL) {
        clock_gettime(CLOCK_MONOTONIC,&t1);
        profile->forward_sec += leopard_elapsed_sec(&t0,&t1);
      }
      for(size_t row=0;row<m;row++) {
        polx rhs[1];
        if(profile != NULL)
          clock_gettime(CLOCK_MONOTONIC,&t0);
        polxvec_sprod(rhs,&matrix[row * len],bnextx,len);
        if(profile != NULL) {
          clock_gettime(CLOCK_MONOTONIC,&t1);
          profile->product_sec += leopard_elapsed_sec(&t0,&t1);
          t0 = t1;
        }
        leopard_public_decompose(&cur[row * block],rhs,block);
        if(profile != NULL) {
          clock_gettime(CLOCK_MONOTONIC,&t1);
          profile->decompose_sec += leopard_elapsed_sec(&t0,&t1);
        }
      }
      poly *tmp = next;
      next = cur;
      cur = tmp;
    }
    return next;
  }

  /* One team for the entire dependent chain, with independent work per layer.
     The pointers are thread-local; each worksharing barrier publishes a layer. */
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    poly *next = buffer1, *cur = buffer0;
    for(size_t t=DEMO_L - 1;t-- > 0;) {
      const polx *matrix = xbits[t] ? A1 : A0;
      if(parallel_conversion) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for(size_t prime=0;prime<K;prime++) {
          polyvec_scale_widening(&bnextx->vec[prime],next,len,K,
                                primes[prime].s,&primes[prime]);
          polyvec_ntt(&bnextx->vec[prime],&bnextx->vec[prime],len,K,&primes[prime]);
        }
      } else {
#ifdef _OPENMP
#pragma omp single
#endif
        polxvec_frompolyvec(bnextx,next,len);
      }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(size_t row=0;row<m;row++) {
        polx rhs[1];
        polxvec_sprod(rhs,&matrix[row * len],bnextx,len);
        leopard_public_decompose(&cur[row * block],rhs,block);
      }
      poly *tmp = next;
      next = cur;
      cur = tmp;
    }
  }
  return (DEMO_L - 1) % 2 ? buffer0 : buffer1;
}

#endif
