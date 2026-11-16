#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "polx.h"
#include "labrador.h"
#include "chihuahua.h"
#include "fips202.h"
#include "jlproj.h"
#include "polz.h"
#include "fips202.h"

static size_t checks;

static void equal(const polx *a, const polx *b, size_t len,
                  const char *op) {
  if(memcmp(a,b,len*sizeof(*a))) {
    fprintf(stderr,"FAIL: %s, len=%zu, LOGQ=%d\n",op,len,LOGQ);
    exit(1);
  }
  checks++;
}

static void reference_mul(polx *r, const polx *a, const polx *b,
                          size_t len, int add) {
  for(size_t i=0;i<K;i++) {
    if(add)
      polyvec_poly_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
    else
      polyvec_poly_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
  }
  if(add) polxvec_reduce(r,len);
}

static void reference_dot(polx *r, const polx *a, const polx *b,
                          size_t len, int add) {
  for(size_t i=0;i<K;i++) {
    if(add)
      polyvec_sprod_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
    else
      polyvec_sprod_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
  }
}

static void check_projection(void) {
  size_t lengths[5] = {17,2060,23,3001,0};
  size_t parts[5] = {0,3,0,3,0};
  const size_t width = 1024, rows = 6, bytes = rows*width*JL_BYTES;
  uint8_t *reference = aligned_alloc(64,bytes);
  uint8_t *candidate = aligned_alloc(64,bytes);
  witness wt[1];
  init_witness_raw(wt,5,lengths);
  size_t max_nonce = 0;
  for(unsigned trial=0;trial<12;trial++) {
    uint8_t seed[16] = {(uint8_t)trial};
    for(size_t i=0;i<5;i++) {
      polyvec_ternary(wt->s[i],lengths[i],seed,i);
      wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],lengths[i]);
    }
    statement st = {.n=width}, other;
    memcpy(st.h,seed,16);
    const statement initial = st;
    proof pi = {.r=5,.n=lengths,.nu=parts}, pi_other;
    int32_t hostile[256];
    for(size_t u=0;u<256;u++) hostile[u]=INT32_MIN;
    if(jlproj_normsq_wide(hostile)!=((unsigned __int128)1<<70) ||
       jlproj_normsq(hostile)!=UINT64_MAX) exit(1);
    checks++;
    omp_set_num_threads(1);
    if(project(&st,&pi,(void *)reference,wt)) exit(1);
    if(pi.jlnonce > max_nonce) max_nonce = pi.jlnonce;
    const int counts[] = {2,3,8};
    for(size_t i=0;i<sizeof(counts)/sizeof(counts[0]);i++) {
      other = initial;
      pi_other = pi;
      memset(candidate,0xa5,bytes);
      omp_set_num_threads(counts[i]);
      if(project(&other,&pi_other,(void *)candidate,wt) ||
         pi.jlnonce != pi_other.jlnonce ||
         memcmp(pi.p,pi_other.p,sizeof(pi.p)) ||
         memcmp(st.h,other.h,16) || memcmp(reference,candidate,bytes)) {
        fprintf(stderr,"FAIL: JL stream/projection/padding/hash, trial=%u threads=%d\n",
                trial,counts[i]);
        exit(1);
      }
      checks++;
      other=initial;
      memset(candidate,0xa5,bytes);
      if(reduce_project(&other,(void *)candidate,&pi,5,JLMAXNORMSQ) ||
         memcmp(st.h,other.h,16) || memcmp(reference,candidate,bytes)) {
        fprintf(stderr,"FAIL: verifier JL stream/padding/hash, trial=%u threads=%d\n",
                trial,counts[i]);
        exit(1);
      }
      checks++;
    }
    if(trial==0) {
      pi_other=pi;
      for(size_t u=0;u<256;u++) pi_other.p[u]=(int32_t)1<<30;
      other=initial;
      if(reduce_project(&other,(void *)candidate,&pi_other,5,JLMAXNORMSQ)==0) {
        fprintf(stderr,"FAIL: overflowing hostile JL norm accepted\n");
        exit(1);
      }
      checks++;
    }
    constraint slow={0}, fast={0};
    init_constraint_raw(&slow,rows,width,1,0);
    init_constraint_raw(&fast,rows,width,1,0);
    statement s1={.r=rows,.n=width}, s2=s1;
    memcpy(s1.h,st.h,16);
    memcpy(s2.h,st.h,16);
    collaps_jlproj_raw(&slow,rows,width,s1.h,pi.p,(void *)reference);
    collaps_jlproj(&fast,&s2,&pi,(void *)reference);
    equal(slow.phi,fast.phi,rows*width,"cached padded Bin1 columns");
    equal(slow.b,fast.b,1,"cached padded Bin1 RHS");
    if(memcmp(s1.h,s2.h,16)) exit(1);
    free_constraint(&slow); free_constraint(&fast);
  }
  printf("JL serial/parallel equality PASS (12 seeds, 2/3/8 threads, max nonce=%zu)\n",max_nonce);
  free_witness(wt);
  free(reference);
  free(candidate);
}

static void check_projection_collapse(void) {
  const size_t sizes[] = {1, 15, 16, 17, 63, 4097};
  const size_t capacity = 4097;
  uint8_t *mat = aligned_alloc(64,capacity*JL_BYTES);
  poly *s = aligned_alloc(64,capacity*sizeof(*s));
  polx *x = aligned_alloc(64,capacity*sizeof(*x));
  polx *phi = aligned_alloc(64,capacity*sizeof(*phi));
  polx *reference = aligned_alloc(64,capacity*sizeof(*reference));
  uint8_t alpha[256*QBYTES], seed[16] = {71};
  int32_t projection[256];
  polx dot[1];
  polz z[1];
  if(!mat || !s || !x || !phi || !reference) abort();
  shake128(mat,capacity*JL_BYTES,seed,sizeof(seed));
  seed[0]++;
  shake128(alpha,sizeof(alpha),seed,sizeof(seed));
  polyvec_ternary(s,capacity,seed,0);
  polxvec_frompolyvec(x,s,capacity);
  for(size_t k=0;k<sizeof(sizes)/sizeof(sizes[0]);k++) {
    size_t n=sizes[k];
    memset(projection,0,sizeof(projection));
    polyvec_jlproj_add(projection,s,n,mat);
    polxvec_jlproj_collapsmat(phi,mat,n,alpha);
    polxvec_jlproj_collapsmat_reference(reference,mat,n,alpha);
    equal(phi,reference,n,"fused Bin1 versus two-plane CRT reference");
    polxvec_sprod(dot,phi,x,n);
    polz_frompolx(z,dot);
    zz value;
    for(size_t j=0;j<L;j++) value.limbs[j]=z->limbs[j].c[0];
    if(scalar_modq(zz_toint128(&value)) !=
       scalar_modq(jlproj_collapsproj(projection,alpha))) {
      fprintf(stderr,"FAIL: Bin1 projection/collapse identity n=%zu\n",n);
      exit(1);
    }
    checks++;
  }
  for(size_t pattern=0;pattern<4;pattern++) {
    if(pattern<3) memset(alpha,pattern==0 ? 0 : pattern==1 ? 255 : 128,sizeof(alpha));
    else shake128(alpha,sizeof(alpha),seed,sizeof(seed));
    for(size_t t=0;t<3;t++) {
      const int threads[]={1,3,8};
      omp_set_num_threads(threads[t]);
      polxvec_jlproj_collapsmat(phi,mat,capacity,alpha);
      polxvec_jlproj_collapsmat_reference(reference,mat,capacity,alpha);
      equal(phi,reference,capacity,"fused Bin1 full-width alpha and thread counts");
    }
  }
  free(mat); free(s); free(x); free(phi); free(reference);
  puts("Bin1 projection/collapse identity PASS");
}

#if LOGQ == 72
static void reference_full_scalar(polx *r, const polx *a, size_t len, labrador_scalar s) {
  const int64_t base=(int64_t)1<<22;
  int64_t digits[5];
  size_t ndigits=0;
  polx *tmp=aligned_alloc(64,(len+1)*sizeof(*tmp));
  do { digits[ndigits++]=(int64_t)(s%base); s/=base; } while(s);
  polxvec_setzero(tmp,len);
  for(size_t j=ndigits;j>0;j--) {
    if(j!=ndigits) {
      polxvec_scale(tmp,tmp,len,base);
      polxvec_refresh(tmp,len);
    }
    if(digits[j-1]) {
      polxvec_scale_add(tmp,a,len,digits[j-1]);
      polxvec_refresh(tmp,len);
    }
  }
  polxvec_add(r,r,tmp,len);
  polxvec_refresh(r,len);
  free(tmp);
}

static void check_full_scalar(void) {
  const size_t capacity=4098;
  const size_t lengths[]={0,1,63,64,65,255,256,257,4097};
  const labrador_scalar half=(((labrador_scalar)1<<72)-107)/2;
  const labrador_scalar scalars[]={0,1,-1,half,-half,half-1,-half+1,((labrador_scalar)1<<64)+123};
  polx *src=aligned_alloc(64,capacity*sizeof(*src));
  polx *want=aligned_alloc(64,capacity*sizeof(*want));
  polx *got=aligned_alloc(64,capacity*sizeof(*got));
  const uint8_t seed[16]={137,19};
  for(size_t trial=0;trial<4;trial++) {
    polxvec_almostuniform(src,capacity,seed,trial);
    if(trial==1) {
      polz z;
      memset(&z,0,sizeof(z));
      for(size_t u=0;u<N;u++) polz_setcoeff_fromint128(&z,u%2 ? half : -half,u);
      for(size_t i=0;i<capacity;i++) polz_topolx(src+i,&z);
    }
    if(trial>=2) polxvec_scale(src,src,capacity,(labrador_scalar)1<<30);
    omp_set_num_threads(trial%2 ? 8 : 1);
    for(size_t j=0;j<sizeof(lengths)/sizeof(*lengths);j++)
      for(size_t k=0;k<sizeof(scalars)/sizeof(*scalars);k++) {
        const size_t len=lengths[j];
        memcpy(want,src,(len+1)*sizeof(*src));
        memcpy(got,src,(len+1)*sizeof(*src));
        reference_full_scalar(want,trial==3 ? want : src,len,scalars[k]);
        polxvec_scale_add_modq(got,trial==3 ? got : src,len,scalars[k]);
        equal(want,got,len+1,"canonical full-width scalar versus Horner, guard and alias");
      }
  }
  free(src); free(want); free(got);
  puts("Full-width scalar/CRT headroom regression PASS");
}
#endif

static void check_aggregation(void) {
  size_t ranks[]={513,1101,100,333}, parts[]={0,1,0,1};
  proof pi={.r=4,.n=ranks,.nu=parts};
  for(unsigned trial=0;trial<3;trial++) {
    sparsecnst cnst[12]={0};
    const uint8_t seed[16]={(uint8_t)trial,41};
    statement got={.r=2,.n=4096}, want=got;
    init_constraint_raw(got.cnst,2,4096,1,0);
    init_constraint_raw(want.cnst,2,4096,1,0);
    polxvec_almostuniform(got.cnst->phi,8192,seed,0);
    memcpy(want.cnst->phi,got.cnst->phi,8192*sizeof(polx));
    polxvec_setzero(got.cnst->b,1);
    polxvec_setzero(want.cnst->b,1);
    memcpy(got.h,seed,16);
    uint8_t hash[32];
    shake128(hash,32,seed,16);
    polxvec_reduce(want.cnst->phi,8192);
    for(size_t i=0;i<12;i++) {
      size_t idx[]={0,1,2,3};
      size_t lens[]={300,890,67,230};
      init_sparsecnst_raw(&cnst[i],4,4,idx,lens,1,0,0);
      cnst[i].scale=malloc(4*sizeof(*cnst[i].scale));
      polxvec_almostuniform(cnst[i].b,1,seed,10+i);
      polx alpha;
      polxvec_quarternary(&alpha,1,hash+16,i);
      if(i!=5) polxvec_sprod_add(want.cnst->b,&alpha,cnst[i].b,1);
      for(size_t j=0;j<4;j++) {
        cnst[i].off[j]=13;
        cnst[i].scale[j]=j==0 ? 1 : j==1 ? -1 : j==2 ? 40960 : 40960L*40960;
        polxvec_almostuniform(cnst[i].phi[j],lens[j],seed,100+4*i+j);
        if(i==5) continue;
        size_t base=0,v=0;
        for(size_t u=0;u<j;u++) {
          base+=parts[u]*4096;
          v=parts[u] ? 0 : v+ranks[u];
        }
        polx factor;
        if(cnst[i].scale[j]==1) factor=alpha;
        else polx_scale(&factor,&alpha,cnst[i].scale[j]);
        polxvec_polx_mul_add(want.cnst->phi+base+v+13,&factor,cnst[i].phi[j],lens[j]);
      }
      if(i==5) cnst[i].deg=0;
    }
    const int threads[]={1,3,8};
    omp_set_num_threads(threads[trial]);
    aggregate_sparsecnst(&got,&pi,cnst,12);
    equal(want.cnst->phi,got.cnst->phi,8192,"bucketed aggregate / serial reference");
    equal(want.cnst->b,got.cnst->b,1,"aggregate RHS");
    if(memcmp(got.h,hash,16)) exit(1);
    for(size_t i=0;i<12;i++) free_sparsecnst(&cnst[i]);
    free_constraint(got.cnst); free_constraint(want.cnst);
  }
  puts("Bucketed sparse aggregation regression PASS");
}

int main(void) {
  check_aggregation();
#if LOGQ == 72
  check_full_scalar();
#endif
  const size_t lengths[] = {0,1,2,31,32,63,64,65,127,128,129,1024,4097,8193};
  const size_t capacity = 8194;
  polx *a = aligned_alloc(64,capacity*sizeof(polx));
  polx *b = aligned_alloc(64,capacity*sizeof(polx));
  polx *want = aligned_alloc(64,capacity*sizeof(polx));
  polx *got = aligned_alloc(64,capacity*sizeof(polx));
  poly *coeffs = aligned_alloc(64,capacity*sizeof(poly));
  uint8_t seed[16] = {0};
  if(!a || !b || !want || !got || !coeffs) return 2;

  for(unsigned trial=0;trial<16;trial++) {
    seed[0] = (uint8_t)trial;
    polxvec_almostuniform(a,capacity,seed,0);
    polxvec_almostuniform(b,capacity,seed,1);
    polyvec_ternary(coeffs,capacity,seed,2);
    if(trial > 1) {
      for(size_t i=0;i<capacity;i++)
        for(size_t u=0;u<N;u++) {
          uint16_t bits = trial == 2 ? 0x8000 : trial == 3 ? 0x7fff :
            (uint16_t)(i*257+u*997+trial*8191);
          memcpy(&coeffs[i].vec->c[u],&bits,sizeof(bits));
        }
    }
    if(trial == 0) polxvec_setzero(a,capacity);
    if(trial == 1) polxvec_setzero(b,capacity);
    for(size_t j=0;j<sizeof(lengths)/sizeof(lengths[0]);j++) {
      const size_t len = lengths[j];
      memcpy(want,b,(len+1)*sizeof(polx));
      memcpy(got,b,(len+1)*sizeof(polx));
      for(size_t limb=0;limb<K;limb++)
        polyvec_scale_widening(&want->vec[limb],coeffs,len,K,primes[limb].s,&primes[limb]);
      polxvec_ntt(want,want,len);
      polxvec_frompolyvec(got,coeffs,len);
      equal(want,got,len+1,"coefficient conversion/NTT (including guard)");
      memcpy(want,b,(len+1)*sizeof(polx));
      memcpy(got,b,(len+1)*sizeof(polx));
      for(size_t limb=0;limb<K;limb++)
        if(polyvec_collaps_add_extension(&want->vec[limb],&a->vec[limb],
                                        &b->vec[limb],len,K,1,&primes[limb])!=len)
          abort();
      if(polxvec_collaps_add_extension(got,a,b,len,1,1)!=len) abort();
      equal(want,got,len+1,"degree-one extension collapse (including guard)");
      memcpy(want,b,(len+1)*sizeof(polx));
      memcpy(got,b,(len+1)*sizeof(polx));
      for(size_t i=0;i<len;i++) polx_refresh(&want[i]);
      polxvec_refresh(got,len);
      equal(want,got,len+1,"parallel CRT refresh (including guard)");
      for(int add=0;add<2;add++) {
        for(int alias=0;alias<3;alias++) {
          memcpy(want,b,(len+1)*sizeof(polx));
          memcpy(got,b,(len+1)*sizeof(polx));
          const polx *wa = alias == 2 ? want : a;
          const polx *ga = alias == 2 ? got : a;
          const polx *wb = alias == 1 ? want : b;
          const polx *gb = alias == 1 ? got : b;
          reference_mul(want,wa,wb,len,add);
          if(add) polxvec_polx_mul_add(got,ga,gb,len);
          else polxvec_polx_mul(got,ga,gb,len);
          equal(want,got,len+1,"broadcast multiply/add (including guard)");
        }
        memcpy(want,a,(len+1)*sizeof(polx));
        memcpy(got,a,(len+1)*sizeof(polx));
        reference_dot(want,want,b,len,add);
        if(add) polxvec_sprod_add(got,got,b,len);
        else polxvec_sprod(got,got,b,len);
        equal(want,got,len+1,"inner product/add, output aliases first input");
      }
    }
  }
  check_projection();
  check_projection_collapse();
  printf("LOGQ=%d: %zu exact arithmetic/matrix/bound checks PASS\n",LOGQ,checks);
  free(a);
  free(b);
  free(want);
  free(got);
  free(coeffs);
  return 0;
}
