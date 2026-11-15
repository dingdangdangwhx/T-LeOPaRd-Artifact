#include <stdint.h>
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "data.h"
#include "polz.h"
#include "jlproj.h"

#ifndef JLPROJ_Q64_MAX_THREADS
#define JLPROJ_Q64_MAX_THREADS 8
#endif
#ifndef JLPROJ_Q64_USE_AVX512
#define JLPROJ_Q64_USE_AVX512 1
#endif

static int64_t jl_center_u64_modq(uint64_t x) {
  const unsigned __int128 q = ((unsigned __int128)1 << LOGQ) - QOFF;
  unsigned __int128 y = x;

  if(y >= q)
    y -= q;
  if(y > q/2)
    return (int64_t)((__int128)y - (__int128)q);
  return (int64_t)y;
}

/* matrix format (32x16):
 * 4 bits -> nibble: 8 row stride
 * 2 nibbles -> byte: 8 col stride
 * 4 bytes -> dword: 1 row stride
 * 2 dwords -> qword: 4 row stride
 * 8 qwords -> zvec: 1 col stride
 */

void poly_jlproj_add(int32_t out[256], const poly *in, const uint8_t mat[256*N/8]) {
  int i;
  __m512i a,b,c,d;
  __m512i e,f,g,h;
  __m512i k,l,m,n;
  __m512i o,p,q,r;
  __m512i s,t,u,v;
  __m512i w,x,y,z;
  const __m512i zero = _mm512_setzero_si512();
  const __m512i stencil1 = _mm512_set1_epi8(0x11);
  const __m512i stencil2 = _mm512_set1_epi8(0x22);
  const __m512i stencil4 = _mm512_set1_epi8(0x44);
  const __m512i stencil8 = _mm512_set1_epi8((char)0x88);
  const __m512i vpermi2bidx = _mm512_set_epi8(95,87,79,71,
                                              94,86,78,70,
                                              93,85,77,69,
                                              92,84,76,68,
                                              91,83,75,67,
                                              90,82,74,66,
                                              89,81,73,65,
                                              88,80,72,64,
                                              31,23,15, 7,
                                              30,22,14, 6,
                                              29,21,13, 5,
                                              28,20,12, 4,
                                              27,19,11, 3,
                                              26,18,10, 2,
                                              25,17, 9, 1,
                                              24,16, 8, 0);

  /* Prepare 16 vectors containing all 16 signed combinations of 4 polynomial coeffs that are 16 apart */
  a = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&in->vec->c[48]));
  b = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&in->vec->c[32]));

  e = _mm512_add_epi32(a,b);     // ++
  f = _mm512_sub_epi32(a,b);     // +-
  g = _mm512_sub_epi32(b,a);     // -+
  h = _mm512_sub_epi32(zero,e);  // --

  a = _mm512_shuffle_i64x2(e,f,0x44);  // ++,+-
  b = _mm512_shuffle_i64x2(e,f,0xEE);  // ++,+-
  c = _mm512_shuffle_i64x2(g,h,0x44);  // -+,--
  d = _mm512_shuffle_i64x2(g,h,0xEE);  // -+,--

  e = _mm512_shuffle_i64x2(a,c,0x88);  // ++,+-,-+,--
  f = _mm512_shuffle_i64x2(a,c,0xDD);  // ++,+-,-+,--
  g = _mm512_shuffle_i64x2(b,d,0x88);  // ++,+-,-+,--
  h = _mm512_shuffle_i64x2(b,d,0xDD);  // ++,+-,-+,--

  a = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&in->vec->c[16]));
  b = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&in->vec->c[ 0]));

  c = _mm512_add_epi32(a,b);  // ++
  d = _mm512_sub_epi32(a,b);  // +-

  k = _mm512_shuffle_i64x2(c,c,0x00);  // ++,++,++,++
  l = _mm512_shuffle_i64x2(c,c,0x55);  // ++,++,++,++
  m = _mm512_shuffle_i64x2(c,c,0xAA);  // ++,++,++,++
  n = _mm512_shuffle_i64x2(c,c,0xFF);  // ++,++,++,++
  o = _mm512_shuffle_i64x2(d,d,0x00);  // +-,+-,+-,+-
  p = _mm512_shuffle_i64x2(d,d,0x55);  // +-,+-,+-,+-
  q = _mm512_shuffle_i64x2(d,d,0xAA);  // +-,+-,+-,+-
  r = _mm512_shuffle_i64x2(d,d,0xFF);  // +-,+-,+-,+-

  s = _mm512_add_epi32(e,k);  // ++++,+-++,-+++,--++
  t = _mm512_add_epi32(e,o);  // +++-,+-+-,-++-,--+-
  u = _mm512_sub_epi32(e,o);  // ++-+,+--+,-+-+,---+
  v = _mm512_sub_epi32(e,k);  // ++--,+---,-+--,----
  w = _mm512_add_epi32(f,l);  // ++++,+-++,-+++,--++
  x = _mm512_add_epi32(f,p);  // +++-,+-+-,-++-,--+-
  y = _mm512_sub_epi32(f,p);  // ++-+,+--+,-+-+,---+
  z = _mm512_sub_epi32(f,l);  // ++--,+---,-+--,----
  a = _mm512_add_epi32(g,m);  // ++++,+-++,-+++,--++
  b = _mm512_add_epi32(g,q);  // +++-,+-+-,-++-,--+-
  c = _mm512_sub_epi32(g,q);  // ++-+,+--+,-+-+,---+
  d = _mm512_sub_epi32(g,m);  // ++--,+---,-+--,----
  e = _mm512_add_epi32(h,n);  // ++++,+-++,-+++,--++
  f = _mm512_add_epi32(h,r);  // +++-,+-+-,-++-,--+-
  g = _mm512_sub_epi32(h,r);  // ++-+,+--+,-+-+,---+
  h = _mm512_sub_epi32(h,n);  // ++--,+---,-+--,----

  k = _mm512_unpacklo_epi32(s,t);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  l = _mm512_unpacklo_epi32(u,v);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  m = _mm512_unpackhi_epi32(s,t);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  n = _mm512_unpackhi_epi32(u,v);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  o = _mm512_unpacklo_epi32(w,x);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  p = _mm512_unpacklo_epi32(y,z);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  q = _mm512_unpackhi_epi32(w,x);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  r = _mm512_unpackhi_epi32(y,z);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  s = _mm512_unpacklo_epi32(a,b);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  t = _mm512_unpacklo_epi32(c,d);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  u = _mm512_unpackhi_epi32(a,b);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  v = _mm512_unpackhi_epi32(c,d);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  w = _mm512_unpacklo_epi32(e,f);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  x = _mm512_unpacklo_epi32(g,h);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  y = _mm512_unpackhi_epi32(e,f);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  z = _mm512_unpackhi_epi32(g,h);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----

  a = _mm512_unpacklo_epi64(k,l);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  b = _mm512_unpackhi_epi64(k,l);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  c = _mm512_unpacklo_epi64(m,n);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  d = _mm512_unpackhi_epi64(m,n);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  e = _mm512_unpacklo_epi64(o,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  f = _mm512_unpackhi_epi64(o,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  g = _mm512_unpacklo_epi64(q,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  h = _mm512_unpackhi_epi64(q,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  k = _mm512_unpacklo_epi64(s,t);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  l = _mm512_unpackhi_epi64(s,t);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  m = _mm512_unpacklo_epi64(u,v);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  n = _mm512_unpackhi_epi64(u,v);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  o = _mm512_unpacklo_epi64(w,x);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  p = _mm512_unpackhi_epi64(w,x);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  q = _mm512_unpacklo_epi64(y,z);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  r = _mm512_unpackhi_epi64(y,z);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----

  for(i=0;i<256/32;i++) {
    s = _mm512_load_si512((__m512i*)&mat[32*16/8*i+256*16/8*0]);
    t = _mm512_load_si512((__m512i*)&mat[32*16/8*i+256*16/8*1]);
    u = _mm512_load_si512((__m512i*)&mat[32*16/8*i+256*16/8*2]);
    v = _mm512_load_si512((__m512i*)&mat[32*16/8*i+256*16/8*3]);

    x = _mm512_slli_epi32(t,1);
    y = _mm512_slli_epi32(u,2);
    z = _mm512_slli_epi32(v,3);
    w = _mm512_and_si512(s,stencil1);
    x = _mm512_and_si512(x,stencil2);
    y = _mm512_and_si512(y,stencil4);
    z = _mm512_and_si512(z,stencil8);
    w = _mm512_add_epi8(w,x);
    x = _mm512_add_epi8(y,z);
    w = _mm512_add_epi8(w,x);

    x = _mm512_srli_epi32(s,1);
    x = _mm512_and_si512(x,stencil1);
    y = _mm512_and_si512(t,stencil2);
    x = _mm512_add_epi8(x,y);
    y = _mm512_slli_epi32(u,1);
    z = _mm512_slli_epi32(v,2);
    y = _mm512_and_si512(y,stencil4);
    z = _mm512_and_si512(z,stencil8);
    y = _mm512_add_epi8(y,z);
    x = _mm512_add_epi8(x,y);

    y = _mm512_permutex2var_epi8(w,vpermi2bidx,x);
    z = _mm512_add_epi8(vpermi2bidx,_mm512_set1_epi8(32));
    z = _mm512_permutex2var_epi8(w,z,x);

#define BLOCK(T0,T1)                  \
  x = _mm512_permutexvar_epi32(y,T0); \
  w = _mm512_add_epi32(w,x);          \
  x = _mm512_permutexvar_epi32(z,T1); \
  y = _mm512_srli_epi32(y,4);         \
  z = _mm512_srli_epi32(z,4);         \
  w = _mm512_add_epi32(w,x)

    w = _mm512_loadu_si512((__m512i*)&out[32*i+ 0]);
    BLOCK(a,e);
    BLOCK(k,o);
    BLOCK(b,f);
    BLOCK(l,p);
    BLOCK(c,g);
    BLOCK(m,q);
    BLOCK(d,h);
    BLOCK(n,r);
    _mm512_storeu_si512((__m512i*)&out[32*i+ 0],w);

    w = _mm512_srli_epi32(s,2);
    x = _mm512_srli_epi32(t,1);
    z = _mm512_slli_epi32(v,1);
    w = _mm512_and_si512(w,stencil1);
    x = _mm512_and_si512(x,stencil2);
    y = _mm512_and_si512(u,stencil4);
    z = _mm512_and_si512(z,stencil8);
    w = _mm512_add_epi8(w,x);
    x = _mm512_add_epi8(y,z);
    w = _mm512_add_epi8(w,x);

    x = _mm512_srli_epi32(s,3);
    y = _mm512_srli_epi32(t,2);
    x = _mm512_and_si512(x,stencil1);
    y = _mm512_and_si512(y,stencil2);
    x = _mm512_add_epi8(x,y);
    y = _mm512_srli_epi32(u,1);
    y = _mm512_and_si512(y,stencil4);
    z = _mm512_and_si512(v,stencil8);
    y = _mm512_add_epi8(y,z);
    x = _mm512_add_epi8(x,y);

    y = _mm512_permutex2var_epi8(w,vpermi2bidx,x);
    z = _mm512_add_epi8(vpermi2bidx,_mm512_set1_epi8(32));
    z = _mm512_permutex2var_epi8(w,z,x);

    w = _mm512_loadu_si512((__m512i*)&out[32*i+16]);
    BLOCK(a,e);
    BLOCK(k,o);
    BLOCK(b,f);
    BLOCK(l,p);
    BLOCK(c,g);
    BLOCK(m,q);
    BLOCK(d,h);
    BLOCK(n,r);
    _mm512_storeu_si512((__m512i*)&out[32*i+16],w);

#undef BLOCK
  }
}

void polyvec_jlproj_add(int32_t r[256], const poly *p, size_t len, const uint8_t *mat) {
  size_t i;

  for(i=0;i<len;i++)
    poly_jlproj_add(r,&p[i],&mat[256*N/8*i]);
}

static void jlproj_collapsmat32(int64_t *row, const uint8_t *mat, size_t len, const int64_t alpha[32]) {
  size_t i;
  __m512i a,b,c,d;
  __m512i e,f,g,h;
  __m512i k,l,m,n;
  __m512i o,p,q,r;
  __m512i s,t,u,v;
  __m512i w,x,y,z;
  const __m512i zero = _mm512_setzero_si512();

  /* Prepare 8 pairs of vectors containing all 16 signed combinations of 4 challenge coeffs that are 8 apart */
  a = _mm512_load_si512((__m512i*)&alpha[24]);
  b = _mm512_load_si512((__m512i*)&alpha[16]);

  e = _mm512_add_epi64(a,b);     // ++
  f = _mm512_sub_epi64(a,b);     // +-
  g = _mm512_sub_epi64(b,a);     // -+
  h = _mm512_sub_epi64(zero,e);  // --

  s = _mm512_shuffle_i64x2(e,f,0x00);  // ++,++,+-,+-
  t = _mm512_shuffle_i64x2(e,f,0x55);  // ++,++,+-,+-
  u = _mm512_shuffle_i64x2(e,f,0xAA);  // ++,++,+-,+-
  v = _mm512_shuffle_i64x2(e,f,0xFF);  // ++,++,+-,+-
  w = _mm512_shuffle_i64x2(g,h,0x00);  // -+,-+,--,--
  x = _mm512_shuffle_i64x2(g,h,0x55);  // -+,-+,--,--
  y = _mm512_shuffle_i64x2(g,h,0xAA);  // -+,-+,--,--
  z = _mm512_shuffle_i64x2(g,h,0xFF);  // -+,-+,--,--

  a = _mm512_load_si512((__m512i*)&alpha[ 8]);
  b = _mm512_load_si512((__m512i*)&alpha[ 0]);

  c = _mm512_add_epi64(a,b);  // ++
  d = _mm512_sub_epi64(b,a);  // -+

  k = _mm512_shuffle_i64x2(c,d,0x44);  // ++,-+
  l = _mm512_shuffle_i64x2(c,d,0xEE);  // ++,-+
  a = _mm512_shuffle_i64x2(k,k,0x88);  // ++,-+,++,-+
  b = _mm512_shuffle_i64x2(k,k,0xDD);  // ++,-+,++,-+
  c = _mm512_shuffle_i64x2(l,l,0x88);  // ++,-+,++,-+
  d = _mm512_shuffle_i64x2(l,l,0xDD);  // ++,-+,++,-+
  e = _mm512_shuffle_i64x2(k,k,0x22);  // -+,++,-+,++
  f = _mm512_shuffle_i64x2(k,k,0x77);  // -+,++,-+,++
  g = _mm512_shuffle_i64x2(l,l,0x22);  // -+,++,-+,++
  h = _mm512_shuffle_i64x2(l,l,0x77);  // -+,++,-+,++

  k = _mm512_add_epi64(s,a);  // ++++,++-+,+-++,+--+
  l = _mm512_add_epi64(t,b);  // ++++,++-+,+-++,+--+
  m = _mm512_add_epi64(u,c);  // ++++,++-+,+-++,+--+
  n = _mm512_add_epi64(v,d);  // ++++,++-+,+-++,+--+
  o = _mm512_sub_epi64(s,e);  // +++-,++--,+-+-,+---
  p = _mm512_sub_epi64(t,f);  // +++-,++--,+-+-,+---
  q = _mm512_sub_epi64(u,g);  // +++-,++--,+-+-,+---
  r = _mm512_sub_epi64(v,h);  // +++-,++--,+-+-,+---
  s = _mm512_add_epi64(w,a);  // -+++,-+-+,--++,+--+
  t = _mm512_add_epi64(x,b);  // -+++,-+-+,--++,+--+
  u = _mm512_add_epi64(y,c);  // -+++,-+-+,--++,+--+
  v = _mm512_add_epi64(z,d);  // -+++,-+-+,--++,+--+
  w = _mm512_sub_epi64(w,e);  // -++-,-+--,--+-,----
  x = _mm512_sub_epi64(x,f);  // -++-,-+--,--+-,----
  y = _mm512_sub_epi64(y,g);  // -++-,-+--,--+-,----
  z = _mm512_sub_epi64(z,h);  // -++-,-+--,--+-,----

  a = _mm512_unpacklo_epi64(k,o);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  b = _mm512_unpackhi_epi64(k,o);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  c = _mm512_unpacklo_epi64(l,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  d = _mm512_unpackhi_epi64(l,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  e = _mm512_unpacklo_epi64(m,q);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  f = _mm512_unpackhi_epi64(m,q);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  g = _mm512_unpacklo_epi64(n,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  h = _mm512_unpackhi_epi64(n,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  k = _mm512_unpacklo_epi64(s,w);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  l = _mm512_unpackhi_epi64(s,w);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  m = _mm512_unpacklo_epi64(t,x);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  n = _mm512_unpackhi_epi64(t,x);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  o = _mm512_unpacklo_epi64(u,y);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  p = _mm512_unpackhi_epi64(u,y);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  q = _mm512_unpacklo_epi64(v,z);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  r = _mm512_unpackhi_epi64(v,z);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----

#if 0
  const __m125i revidx = _mm512_set_epi64(0,1,2,3,4,5,6,7);

  a = _mm512_loadu_si256((__m512i*)&alpha[ 0]));
  b = _mm512_loadu_si256((__m512i*)&alpha[ 8]));
  c = _mm512_loadu_si256((__m512i*)&alpha[16]));
  d = _mm512_loadu_si256((__m512i*)&alpha[24]));

  e = _mm512_add_epi64(a,b);  // ++ (a+b)
  f = _mm512_sub_epi64(a,b);  // +- (a-b)
  g = _mm512_add_epi64(c,d);  // ++ (c+d)
  h = _mm512_sub_epi64(c,d);  // +- (c-d)

  k = _mm512_add_epi64(e,g);  // ++++
  l = _mm512_add_epi64(e,h);  // +++-
  m = _mm512_sub_epi64(e,h);  // ++-+
  n = _mm512_sub_epi64(e,g);  // ++--
  o = _mm512_add_epi64(f,g);  // +-++
  p = _mm512_add_epi64(f,h);  // +-+-
  q = _mm512_sub_epi64(f,h);  // +--+
  r = _mm512_sub_epi64(f,g);  // +---

  a = _mm512_shuffle_i64x2(k,m,0x44);  // ++++,++-+
  b = _mm512_shuffle_i64x2(k,m,0xEE);  // ++++,++-+
  c = _mm512_shuffle_i64x2(l,n,0x44);  // +++-,++--
  d = _mm512_shuffle_i64x2(l,n,0xEE);  // +++-,++--
  e = _mm512_shuffle_i64x2(o,q,0x44);  // +-++,+--+
  f = _mm512_shuffle_i64x2(o,q,0xEE);  // +-++,+--+
  g = _mm512_shuffle_i64x2(p,r,0x44);  // +-+-,+---
  h = _mm512_shuffle_i64x2(p,r,0xEE);  // +-+-,+---

  k = _mm512_shuffle_i64x2(a,e,0x88);  // ++++,++-+,+-++,+--+
  l = _mm512_shuffle_i64x2(a,e,0xDD);  // ++++,++-+,+-++,+--+
  m = _mm512_shuffle_i64x2(b,f,0x88);  // ++++,++-+,+-++,+--+
  n = _mm512_shuffle_i64x2(b,f,0xDD);  // ++++,++-+,+-++,+--+
  o = _mm512_shuffle_i64x2(c,g,0x88);  // +++-,++--,+-+-,+---
  p = _mm512_shuffle_i64x2(c,g,0xDD);  // +++-,++--,+-+-,+---
  q = _mm512_shuffle_i64x2(d,h,0x88);  // +++-,++--,+-+-,+---
  r = _mm512_shuffle_i64x2(d,h,0xDD);  // +++-,++--,+-+-,+---

  a = _mm512_unpacklo_epi64(k,o);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  b = _mm512_unpackhi_epi64(k,o);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  c = _mm512_unpacklo_epi64(l,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  d = _mm512_unpackhi_epi64(l,p);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  e = _mm512_unpacklo_epi64(m,q);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  f = _mm512_unpackhi_epi64(m,q);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  g = _mm512_unpacklo_epi64(n,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---
  h = _mm512_unpackhi_epi64(n,r);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---

  k = _mm512_sub_epi64(zero,a);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  l = _mm512_sub_epi64(zero,b);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  m = _mm512_sub_epi64(zero,c);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  n = _mm512_sub_epi64(zero,d);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  o = _mm512_sub_epi64(zero,e);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  p = _mm512_sub_epi64(zero,f);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  q = _mm512_sub_epi64(zero,g);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++
  r = _mm512_sub_epi64(zero,h);  // ----,---+,--+-,--++,-+--,-+-+,-++-,-+++

  k = _mm512_permutexvar_epi64(revidx,k);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  l = _mm512_permutexvar_epi64(revidx,l);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  m = _mm512_permutexvar_epi64(revidx,m);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  n = _mm512_permutexvar_epi64(revidx,n);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  o = _mm512_permutexvar_epi64(revidx,o);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  p = _mm512_permutexvar_epi64(revidx,p);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  q = _mm512_permutexvar_epi64(revidx,q);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
  r = _mm512_permutexvar_epi64(revidx,r);  // -+++,-++-,-+-+,-+--,--++,--+-,---+,----
#endif

#define BLOCK(T0,T1)                      \
  u = _mm512_permutex2var_epi64(T0,s,T1); \
  v = _mm512_permutex2var_epi64(T0,t,T1); \
  s = _mm512_srli_epi64(s,4);             \
  t = _mm512_srli_epi64(t,4);             \
  w = _mm512_add_epi64(w,u);              \
  y = _mm512_add_epi64(y,v);              \
  u = _mm512_permutex2var_epi64(T0,s,T1); \
  v = _mm512_permutex2var_epi64(T0,t,T1); \
  s = _mm512_srli_epi64(s,4);             \
  t = _mm512_srli_epi64(t,4);             \
  x = _mm512_add_epi64(x,u);              \
  z = _mm512_add_epi64(z,v)

  for(i=0;i<len/32;i++) {
    s = _mm512_load_si512((__m512i*)&mat[256*32/8*i+       0]);
    t = _mm512_load_si512((__m512i*)&mat[256*32/8*i+256*16/8]);
    w = _mm512_load_si512((__m512i*)&row[32*i+ 0]);
    x = _mm512_load_si512((__m512i*)&row[32*i+ 8]);
    y = _mm512_load_si512((__m512i*)&row[32*i+16]);
    z = _mm512_load_si512((__m512i*)&row[32*i+24]);
    BLOCK(a,k);
    BLOCK(b,l);
    BLOCK(c,m);
    BLOCK(d,n);
    BLOCK(e,o);
    BLOCK(f,p);
    BLOCK(g,q);
    BLOCK(h,r);
    _mm512_store_si512((__m512i*)&row[32*i+ 0],w);
    _mm512_store_si512((__m512i*)&row[32*i+ 8],x);
    _mm512_store_si512((__m512i*)&row[32*i+16],y);
    _mm512_store_si512((__m512i*)&row[32*i+24],z);
  }

#undef BLOCK
}

#if LOGQ == 64
static int64_t jl_center_i128_modq(__int128 x) {
  const __int128 q = ((__int128)1 << LOGQ) - QOFF;

  x %= q;
  if(x > q/2)
    x -= q;
  if(x < -q/2)
    x += q;
  return (int64_t)x;
}

static __int128 jl_combo_q64(uint8_t nib, size_t rowlow, const int64_t alpha[32]) {
  size_t b;
  __int128 acc = 0;

  for(b=0;b<4;b++) {
    const __int128 a = alpha[rowlow + 8*b];
    if(nib & (uint8_t)(1u << b))
      acc -= a;
    else
      acc += a;
  }
  return acc;
}

static void jlproj_collapsmat32_q64(__int128 *row, const uint8_t *mat, size_t len, const int64_t alpha[32]) {
  size_t chunk, lane, rowlow, nib;
  __int128 combo[16][8];

  for(nib=0;nib<16;nib++)
    for(rowlow=0;rowlow<8;rowlow++)
      combo[nib][rowlow] = jl_combo_q64((uint8_t)nib,rowlow,alpha);

  for(chunk=0;chunk<len/32;chunk++) {
    const size_t base = (256*32/8)*chunk;
    for(lane=0;lane<8;lane++) {
      uint64_t s, t;
      memcpy(&s,&mat[base + 8*lane],8);
      memcpy(&t,&mat[base + 256*16/8 + 8*lane],8);
      for(rowlow=0;rowlow<8;rowlow++) {
        const uint8_t s0 = (uint8_t)((s >> (8*rowlow)) & 0x0f);
        const uint8_t s1 = (uint8_t)((s >> (8*rowlow + 4)) & 0x0f);
        const uint8_t t0 = (uint8_t)((t >> (8*rowlow)) & 0x0f);
        const uint8_t t1 = (uint8_t)((t >> (8*rowlow + 4)) & 0x0f);
        const size_t off = 32*chunk;

        row[off + lane] += combo[s0][rowlow];
        row[off + lane + 8] += combo[s1][rowlow];
        row[off + lane + 16] += combo[t0][rowlow];
        row[off + lane + 24] += combo[t1][rowlow];
      }
    }
  }
}

#if JLPROJ_Q64_USE_AVX512
static int64_t jl_q64_reduce_parts(uint64_t lo, int64_t hi) {
  const __int128 x = (__int128)(unsigned __int128)lo + (__int128)hi * QOFF;
  return jl_center_i128_modq(x);
}

typedef struct {
  __m512i lo;
  __m512i hi;
} jl128;

static inline jl128 jl128_load(const int64_t *a) {
  jl128 r;
  r.lo = _mm512_load_si512((__m512i*)a);
  r.hi = _mm512_srai_epi64(r.lo,63);
  return r;
}

static inline jl128 jl128_load_parts(const uint64_t *lo, const int64_t *hi) {
  jl128 r;
  r.lo = _mm512_load_si512((__m512i*)lo);
  r.hi = _mm512_load_si512((__m512i*)hi);
  return r;
}

static inline void jl128_store_parts(uint64_t *lo, int64_t *hi, jl128 a) {
  _mm512_store_si512((__m512i*)lo,a.lo);
  _mm512_store_si512((__m512i*)hi,a.hi);
}

static inline jl128 jl128_add(jl128 a, jl128 b) {
  const __m512i one = _mm512_set1_epi64(1);
  jl128 r;
  r.lo = _mm512_add_epi64(a.lo,b.lo);
  r.hi = _mm512_add_epi64(a.hi,b.hi);
  r.hi = _mm512_mask_add_epi64(r.hi,_mm512_cmplt_epu64_mask(r.lo,a.lo),r.hi,one);
  return r;
}

static inline jl128 jl128_sub(jl128 a, jl128 b) {
  const __m512i one = _mm512_set1_epi64(1);
  jl128 r;
  r.lo = _mm512_sub_epi64(a.lo,b.lo);
  r.hi = _mm512_sub_epi64(a.hi,b.hi);
  r.hi = _mm512_mask_sub_epi64(r.hi,_mm512_cmplt_epu64_mask(a.lo,b.lo),r.hi,one);
  return r;
}

#define JL128_SHUFFLE(A,B,IMM) ((jl128){_mm512_shuffle_i64x2((A).lo,(B).lo,(IMM)), \
                                        _mm512_shuffle_i64x2((A).hi,(B).hi,(IMM))})
#define JL128_UNPACKLO(A,B) ((jl128){_mm512_unpacklo_epi64((A).lo,(B).lo), \
                                    _mm512_unpacklo_epi64((A).hi,(B).hi)})
#define JL128_UNPACKHI(A,B) ((jl128){_mm512_unpackhi_epi64((A).lo,(B).lo), \
                                    _mm512_unpackhi_epi64((A).hi,(B).hi)})
#define JL128_PERMUTEX2(A,IDX,B) ((jl128){_mm512_permutex2var_epi64((A).lo,(IDX),(B).lo), \
                                          _mm512_permutex2var_epi64((A).hi,(IDX),(B).hi)})

static void jlproj_collapsmat32_q64_avx2(uint64_t *rowlo, int64_t *rowhi,
                                         const uint8_t *mat, size_t len,
                                         const int64_t alpha[32]) {
  size_t i;
  jl128 va,vb,vc,vd,ve,vf,vg,vh;
  jl128 vk,vl,vm,vn,vo,vp,vq,vr;
  jl128 vs,vt,vu,vv,vw,vx,vy,vz;
  jl128 rw,rx,ry,rz,tu,tv2;
  const jl128 zero = {_mm512_setzero_si512(),_mm512_setzero_si512()};

  va = jl128_load(&alpha[24]);
  vb = jl128_load(&alpha[16]);

  ve = jl128_add(va,vb);
  vf = jl128_sub(va,vb);
  vg = jl128_sub(vb,va);
  vh = jl128_sub(zero,ve);

  vs = JL128_SHUFFLE(ve,vf,0x00);
  vt = JL128_SHUFFLE(ve,vf,0x55);
  vu = JL128_SHUFFLE(ve,vf,0xAA);
  vv = JL128_SHUFFLE(ve,vf,0xFF);
  vw = JL128_SHUFFLE(vg,vh,0x00);
  vx = JL128_SHUFFLE(vg,vh,0x55);
  vy = JL128_SHUFFLE(vg,vh,0xAA);
  vz = JL128_SHUFFLE(vg,vh,0xFF);

  va = jl128_load(&alpha[8]);
  vb = jl128_load(&alpha[0]);

  vc = jl128_add(va,vb);
  vd = jl128_sub(vb,va);

  vk = JL128_SHUFFLE(vc,vd,0x44);
  vl = JL128_SHUFFLE(vc,vd,0xEE);
  va = JL128_SHUFFLE(vk,vk,0x88);
  vb = JL128_SHUFFLE(vk,vk,0xDD);
  vc = JL128_SHUFFLE(vl,vl,0x88);
  vd = JL128_SHUFFLE(vl,vl,0xDD);
  ve = JL128_SHUFFLE(vk,vk,0x22);
  vf = JL128_SHUFFLE(vk,vk,0x77);
  vg = JL128_SHUFFLE(vl,vl,0x22);
  vh = JL128_SHUFFLE(vl,vl,0x77);

  vk = jl128_add(vs,va);
  vl = jl128_add(vt,vb);
  vm = jl128_add(vu,vc);
  vn = jl128_add(vv,vd);
  vo = jl128_sub(vs,ve);
  vp = jl128_sub(vt,vf);
  vq = jl128_sub(vu,vg);
  vr = jl128_sub(vv,vh);
  vs = jl128_add(vw,va);
  vt = jl128_add(vx,vb);
  vu = jl128_add(vy,vc);
  vv = jl128_add(vz,vd);
  vw = jl128_sub(vw,ve);
  vx = jl128_sub(vx,vf);
  vy = jl128_sub(vy,vg);
  vz = jl128_sub(vz,vh);

  va = JL128_UNPACKLO(vk,vo);
  vb = JL128_UNPACKHI(vk,vo);
  vc = JL128_UNPACKLO(vl,vp);
  vd = JL128_UNPACKHI(vl,vp);
  ve = JL128_UNPACKLO(vm,vq);
  vf = JL128_UNPACKHI(vm,vq);
  vg = JL128_UNPACKLO(vn,vr);
  vh = JL128_UNPACKHI(vn,vr);
  vk = JL128_UNPACKLO(vs,vw);
  vl = JL128_UNPACKHI(vs,vw);
  vm = JL128_UNPACKLO(vt,vx);
  vn = JL128_UNPACKHI(vt,vx);
  vo = JL128_UNPACKLO(vu,vy);
  vp = JL128_UNPACKHI(vu,vy);
  vq = JL128_UNPACKLO(vv,vz);
  vr = JL128_UNPACKHI(vv,vz);

#define BLOCK128(T0,T1)                  \
  tu = JL128_PERMUTEX2((T0),ms,(T1));    \
  tv2 = JL128_PERMUTEX2((T0),mt,(T1));   \
  ms = _mm512_srli_epi64(ms,4);          \
  mt = _mm512_srli_epi64(mt,4);          \
  rw = jl128_add(rw,tu);                 \
  ry = jl128_add(ry,tv2);                \
  tu = JL128_PERMUTEX2((T0),ms,(T1));    \
  tv2 = JL128_PERMUTEX2((T0),mt,(T1));   \
  ms = _mm512_srli_epi64(ms,4);          \
  mt = _mm512_srli_epi64(mt,4);          \
  rx = jl128_add(rx,tu);                 \
  rz = jl128_add(rz,tv2)

  for(i=0;i<len/32;i++) {
    __m512i ms = _mm512_load_si512((__m512i*)&mat[256*32/8*i]);
    __m512i mt = _mm512_load_si512((__m512i*)&mat[256*32/8*i + 256*16/8]);
    rw = jl128_load_parts(&rowlo[32*i],&rowhi[32*i]);
    rx = jl128_load_parts(&rowlo[32*i + 8],&rowhi[32*i + 8]);
    ry = jl128_load_parts(&rowlo[32*i + 16],&rowhi[32*i + 16]);
    rz = jl128_load_parts(&rowlo[32*i + 24],&rowhi[32*i + 24]);
    BLOCK128(va,vk);
    BLOCK128(vb,vl);
    BLOCK128(vc,vm);
    BLOCK128(vd,vn);
    BLOCK128(ve,vo);
    BLOCK128(vf,vp);
    BLOCK128(vg,vq);
    BLOCK128(vh,vr);
    jl128_store_parts(&rowlo[32*i],&rowhi[32*i],rw);
    jl128_store_parts(&rowlo[32*i + 8],&rowhi[32*i + 8],rx);
    jl128_store_parts(&rowlo[32*i + 16],&rowhi[32*i + 16],ry);
    jl128_store_parts(&rowlo[32*i + 24],&rowhi[32*i + 24],rz);
  }

#undef BLOCK128
}

#undef JL128_SHUFFLE
#undef JL128_UNPACKLO
#undef JL128_UNPACKHI
#undef JL128_PERMUTEX2
#endif

static void polxvec_jlproj_collapsmat_q64(polx *r, const uint8_t *mat, size_t len, const int64_t alpha[256]) {
  const size_t blocks = (len + 15)/16;

#if defined(_OPENMP)
  const int q64_threads = MIN(JLPROJ_Q64_MAX_THREADS,omp_get_max_threads());
#pragma omp parallel for schedule(static) num_threads(q64_threads)
#endif
  for(size_t block=0; block<blocks; block++) {
    const size_t blen = MIN((size_t)16,len - 16*block);
    const uint8_t *bmat = mat + block*256*16*N/8;
    polx *br = &r[16*block];
    __attribute__((aligned(64))) int64_t raw64[16*N];
    polz t[16];
    size_t i;

#if JLPROJ_Q64_USE_AVX512
    __attribute__((aligned(64))) uint64_t rawlo[16*N];
    __attribute__((aligned(64))) int64_t rawhi[16*N];
    for(i=0;i<blen*N;i++) {
      rawlo[i] = 0;
      rawhi[i] = 0;
    }
    for(i=0;i<256/32;i++)
      jlproj_collapsmat32_q64_avx2(rawlo,rawhi,&bmat[32*16/8*i],blen*N,&alpha[32*i]);
    for(i=0;i<blen*N;i++)
      raw64[i] = jl_q64_reduce_parts(rawlo[i],rawhi[i]);
#else
    __attribute__((aligned(64))) __int128 raw128[16*N];
    for(i=0;i<blen*N;i++)
      raw128[i] = 0;
    for(i=0;i<256/32;i++)
      jlproj_collapsmat32_q64(raw128,&bmat[32*16/8*i],blen*N,&alpha[32*i]);
    for(i=0;i<blen*N;i++)
      raw64[i] = jl_center_i128_modq(raw128[i]);
#endif

    polzvec_fromint64vec(t,blen,1,raw64);
    polzvec_sigmam1(t,t,blen);
    polzvec_topolxvec(br,t,blen);
  }
}
#endif

#if LOGQ <= 64
static void expand_challenge(int64_t alpha[256], const uint8_t buf[256*QBYTES]) {
  int i;
#if QBYTES >= 3 && QBYTES <= 6
  __m512i f;
#endif

#if QBYTES == 3
  const __mmask64 mask = _cvtu64_mask64(0x0B0B0B0B0B0B0B0Bull);
  const __m512i vpermbidx = _mm512_set_epi8(-1,-1,-1,-1,-1,23,22,21,
                                            -1,-1,-1,-1,-1,20,19,18,
                                            -1,-1,-1,-1,-1,17,16,15,
                                            -1,-1,-1,-1,-1,14,13,12,
                                            -1,-1,-1,-1,-1,11,10, 9,
                                            -1,-1,-1,-1,-1, 8, 7, 6,
                                            -1,-1,-1,-1,-1, 5, 4, 3,
                                            -1,-1,-1,-1,-1, 2, 1, 0);
  for(i=0;i<256/8;i++) {
    f = _mm512_loadu_si512((__m512i*)&buf[24*i]);
    f = _mm512_maskz_permutexvar_epi8(mask,vpermbidx,f);
    _mm512_store_si512((__m512i*)&alpha[8*i],f);
  }
#elif QBYTES == 4
  for(i=0;i<256/8;i++) {
    f = _mm512_cvtepu32_epi64(_mm256_load_si256((__m256i*)&buf[32*i]));
    _mm512_store_si512((__m512i*)&alpha[8*i],f);
  }
#elif QBYTES == 5
  const __mmask64 mask = _cvtu64_mask64(0x1F1F1F1F1F1F1F1Full);
  const __m512i vpermbidx = _mm512_set_epi8(-1,-1,-1,39,38,37,36,35,
                                            -1,-1,-1,34,33,32,31,30,
                                            -1,-1,-1,29,28,27,26,25,
                                            -1,-1,-1,24,23,22,21,20,
                                            -1,-1,-1,19,18,17,16,15,
                                            -1,-1,-1,14,13,12,11,10,
                                            -1,-1,-1, 9, 8, 7, 6, 5,
                                            -1,-1,-1, 4, 3, 2, 1, 0);
  for(i=0;i<256/8;i++) {
    f = _mm512_loadu_si512((__m512i*)&buf[40*i]);
    f = _mm512_maskz_permutexvar_epi8(mask,vpermbidx,f);
    _mm512_store_si512((__m512i*)&alpha[8*i],f);
  }
#elif QBYTES == 6
  const __mmask64 mask = _cvtu64_mask64(0x3F3F3F3F3F3F3F3Full);
  const __m512i vpermbidx = _mm512_set_epi8(-1,-1,47,46,45,44,43,42,
                                            -1,-1,41,40,39,38,37,36,
                                            -1,-1,35,34,33,32,31,30,
                                            -1,-1,29,28,27,26,25,24,
                                            -1,-1,23,22,21,20,19,18,
                                            -1,-1,17,16,15,14,13,12,
                                            -1,-1,11,10, 9, 8, 7, 6,
                                            -1,-1, 5, 4, 3, 2, 1, 0);
  for(i=0;i<256/8;i++) {
    f = _mm512_loadu_si512((__m512i*)&buf[48*i]);
    f = _mm512_maskz_permutexvar_epi8(mask,vpermbidx,f);
    _mm512_store_si512((__m512i*)&alpha[8*i],f);
  }
#elif QBYTES == 7
  int j;
  for(i=0;i<256;i++) {
    uint64_t x = 0;
    for(j=0;j<7;j++)
      x |= ((uint64_t)buf[7*i + j]) << (8*j);
    alpha[i] = (int64_t)x;
  }
#elif QBYTES == 8
  for(i=0;i<256;i++) {
    uint64_t x;
    memcpy(&x,&buf[8*i],8);
    alpha[i] = jl_center_u64_modq(x);
  }
#else
#error
#endif
}

void polxvec_jlproj_collapsmat(polx *r, const uint8_t *mat, size_t len, const uint8_t buf[256*QBYTES]) {
#if LOGQ != 64
  const __m512i qoff = _mm512_set1_epi64(QOFF);
  const __m512i mask = _mm512_set1_epi64(((uint64_t)1 << LOGQ) - 1);
  const __m512i mask14 = _mm512_set1_epi64((1 << 14) - 1);
  const __mmask32 oneinfour = _cvtu32_mask32(0x11111111);
#endif
  __attribute__((aligned(64)))
  int64_t alpha[256];

  expand_challenge(alpha,buf);

#if LOGQ == 64
  polxvec_jlproj_collapsmat_q64(r,mat,len,alpha);
  return;
#else
  const size_t blocks = (len + 15)/16;

#if defined(_OPENMP)
  const int jl_threads = MIN(JLPROJ_Q64_MAX_THREADS,omp_get_max_threads());
#pragma omp parallel for schedule(static) num_threads(jl_threads)
#endif
  for(size_t block=0; block<blocks; block++) {
    const size_t blen = MIN((size_t)16,len - 16*block);
    const uint8_t *bmat = mat + block*256*16*N/8;
    polx *br = &r[16*block];
    __attribute__((aligned(64))) int64_t raw[16*N];
    polz t[16];
    __m512i f,g;
    size_t i,j,k;

    for(i=0;i<blen*N;i++)
      raw[i] = 0;
    for(i=0;i<256/32;i++)
      jlproj_collapsmat32(raw,&bmat[32*16/8*i],blen*N,&alpha[32*i]);

    for(i=0;i<blen;i++) {
      for(j=0;j<N/8;j++) {
        f = _mm512_load_si512((__m512i*)&raw[N*i+8*j]);
        g = _mm512_srai_epi64(f,LOGQ);
        f = _mm512_and_si512(f,mask);
        g = _mm512_mul_epi32(g,qoff);
        f = _mm512_add_epi64(f,g);
        for(k=0;k<L-1;k++) {
          g = _mm512_and_si512(f,mask14);
          f = _mm512_srai_epi64(f,14);
          _mm512_mask_compressstoreu_epi16(&t[i].limbs[k].c[8*j],oneinfour,g);
        }
        _mm512_mask_compressstoreu_epi16(&t[i].limbs[L-1].c[8*j],oneinfour,f);
      }
    }
    polzvec_sigmam1(t,t,blen);
    polzvec_topolxvec(br,t,blen);
  }
#endif
}

static inline int64_t _mm512_hsum_epi64(__m512i a) {
  __m256i t;
  __m128i u;
  t = _mm256_add_epi64(_mm512_castsi512_si256(a),_mm512_extracti64x4_epi64(a,1));
  u = _mm_add_epi64(_mm256_castsi256_si128(t),_mm256_extracti64x2_epi64(t,1));
  u = _mm_add_epi64(u,_mm_unpackhi_epi64(u,u));
  return _mm_cvtsi128_si64(u);
}

int64_t jlproj_collapsproj(const int32_t p[256], const uint8_t buf[256*QBYTES]) {
  int i;
  __attribute__((aligned(64)))
  int64_t alpha[256];

  expand_challenge(alpha,buf);

#if LOGQ == 64
  const __int128 q = ((__int128)1 << LOGQ) - QOFF;
  __int128 acc = 0;

  for(i=0;i<256;i++)
    acc += (__int128)p[i] * alpha[i];
  acc %= q;
  if(acc > q/2)
    acc -= q;
  if(acc < -q/2)
    acc += q;
  return (int64_t)acc;
#else
  int64_t r,t;
  __m512i f,g,h,s;
  const __m512i mask = _mm512_set1_epi64(((uint64_t)1 << LOGQ) - 1);
  const __m512i qoff = _mm512_set1_epi64(QOFF);

  s = _mm512_setzero_si512();
  for(i=0;i<256/8;i++) {
    f = _mm512_cvtepi32_epi64(_mm256_loadu_si256((__m256i*)&p[8*i]));
    g = _mm512_load_si512((__m512i*)&alpha[8*i]);
    h = _mm512_mul_epi32(f,g);
    s = _mm512_add_epi64(s,h);
#if LOGQ >= 32
    h = _mm512_srli_epi32(g,31);
    g = _mm512_srli_epi64(g,32);
    g = _mm512_add_epi64(g,h);
    f = _mm512_mul_epi32(f,g);
    g = _mm512_srai_epi64(f,LOGQ-32);
    f = _mm512_slli_epi64(f,32);
    f = _mm512_and_si512(f,mask);
    s = _mm512_add_epi64(s,f);
    g = _mm512_mul_epi32(g,qoff);
    s = _mm512_add_epi64(s,g);
#endif
    f = _mm512_srai_epi64(s,LOGQ);
    s = _mm512_and_si512(s,mask);
    f = _mm512_mul_epi32(f,qoff);
    s = _mm512_add_epi64(s,f);
  }

  r = _mm512_hsum_epi64(s);
  t = r >> LOGQ;
  r &= ((int64_t)1 << LOGQ) - 1;
  r += t*QOFF;
  return r;
#endif
}

#else
#ifndef JLPROJ_Q72_USE_AVX512
#define JLPROJ_Q72_USE_AVX512 1
#endif
static __int128 jl_center_wide(__int128 x) {
  const __int128 q=((__int128)1<<LOGQ)-QOFF;
  x%=q;
  if(x>q/2) x-=q;
  if(x<-q/2) x+=q;
  return x;
}

static void jl_expand_wide(__int128 alpha[256], const uint8_t *buf) {
  for(size_t i=0;i<256;i++) {
    unsigned __int128 x=0;
    for(size_t j=0;j<QBYTES;j++) x|=(unsigned __int128)buf[i*QBYTES+j]<<(8*j);
    alpha[i]=jl_center_wide((__int128)x);
  }
}

void polxvec_jlproj_collapsmat(polx *r, const uint8_t *mat, size_t len,
                              const uint8_t buf[256*QBYTES]) {
  __int128 alpha[256];
  jl_expand_wide(alpha,buf);
#if JLPROJ_Q72_USE_AVX512
  _Alignas(64) int64_t low[256], high[256];
  for(size_t i=0;i<256;i++) {
    low[i]=(uint64_t)alpha[i]&(((uint64_t)1<<48)-1);
    high[i]=(int64_t)((alpha[i]-low[i])/((__int128)1<<48));
  }
#else
  __int128 combo[8][16][8];
  for(size_t group=0;group<8;group++)
    for(size_t nib=0;nib<16;nib++)
      for(size_t row=0;row<8;row++) {
        __int128 sum=0;
        for(size_t b=0;b<4;b++) sum+=(nib & (1u<<b)) ? -alpha[group*32+row+8*b] : alpha[group*32+row+8*b];
        combo[group][nib][row]=sum;
      }
#endif
  const size_t blocks=(len+15)/16;
#ifdef _OPENMP
  const int threads=MIN(JLPROJ_Q64_MAX_THREADS,omp_get_max_threads());
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
  for(size_t block=0;block<blocks;block++) {
    const size_t count=MIN((size_t)16,len-16*block);
    const uint8_t *bmat=mat+block*256*16*N/8;
    polz z[16];
#if JLPROJ_Q72_USE_AVX512
    /* Each signed sum has 256 terms: low < 2^56, high < 2^32.
       Both fit int64_t, so the existing AVX512 kernel needs no carry logic. */
    _Alignas(64) int64_t rowlow[16*N]={0}, rowhigh[16*N]={0};
    for(size_t group=0;group<8;group++) {
      const uint8_t *gmat=bmat+32*16/8*group;
      jlproj_collapsmat32(rowlow,gmat,count*N,&low[32*group]);
      jlproj_collapsmat32(rowhigh,gmat,count*N,&high[32*group]);
    }
    for(size_t i=0;i<count*N;i++) {
      __int128 value=(__int128)rowlow[i]+(__int128)rowhigh[i]*((__int128)1<<48);
      polz_setcoeff_fromint128(&z[i/N],jl_center_wide(value),i%N);
    }
#else
    __int128 raw[16*N]={0};
    for(size_t group=0;group<8;group++) {
      const uint8_t *gmat=bmat+32*16/8*group;
      for(size_t chunk=0;chunk<count*N/32;chunk++)
        for(size_t lane=0;lane<8;lane++) {
          uint64_t s,t;
          memcpy(&s,gmat+256*32/8*chunk+8*lane,8);
          memcpy(&t,gmat+256*32/8*chunk+256*16/8+8*lane,8);
          for(size_t row=0;row<8;row++) {
            raw[32*chunk+lane]+=combo[group][(s>>(8*row))&15][row];
            raw[32*chunk+lane+8]+=combo[group][(s>>(8*row+4))&15][row];
            raw[32*chunk+lane+16]+=combo[group][(t>>(8*row))&15][row];
            raw[32*chunk+lane+24]+=combo[group][(t>>(8*row+4))&15][row];
          }
        }
    }
    for(size_t i=0;i<count*N;i++) polz_setcoeff_fromint128(&z[i/N],jl_center_wide(raw[i]),i%N);
#endif
    polzvec_sigmam1(z,z,count);
    polzvec_topolxvec(&r[16*block],z,count);
  }
}

__int128 jlproj_collapsproj(const int32_t p[256], const uint8_t buf[256*QBYTES]) {
  __int128 alpha[256],sum=0;
  jl_expand_wide(alpha,buf);
  for(size_t i=0;i<256;i++) sum+=(__int128)p[i]*alpha[i];
  return jl_center_wide(sum);
}

static inline int64_t _mm512_hsum_epi64(__m512i a) {
  return _mm512_reduce_add_epi64(a);
}
#endif

uint64_t jlproj_normsq(const int32_t p[256]) {
  int i;
  __m512i f,g,s;

  s = _mm512_setzero_si512();
  for(i=0;i<256/16;i++) {
    f = _mm512_loadu_si512(&p[16*i]);
    g = _mm512_mul_epi32(f,f);
    f = _mm512_srli_epi64(f,32);
    f = _mm512_mul_epi32(f,f);
    f = _mm512_add_epi64(f,g);
    s = _mm512_add_epi64(f,s);
  }

  return _mm512_hsum_epi64(s);
}
