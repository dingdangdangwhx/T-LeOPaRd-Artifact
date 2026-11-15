#ifndef LEOPARD_SHORT_CODEC_H
#define LEOPARD_SHORT_CODEC_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* A block is at most 512 signed coefficients. Rice is used only when it
 * beats min-offset packing; ties have one deterministic representation. */
static uint32_t leopard_short_zigzag(int16_t x) {
  return x < 0 ? (uint32_t)(-2 * (int32_t)x - 1) : 2u * (uint32_t)x;
}

static size_t leopard_short_plan(const int16_t *v, size_t n,
                                 unsigned *mode, int16_t *minv) {
  int32_t lo=v[0], hi=v[0];
  for(size_t i=1;i<n;i++) {
    if(v[i]<lo) lo=v[i];
    if(v[i]>hi) hi=v[i];
  }
  unsigned width=0;
  for(uint32_t span=(uint32_t)(hi-lo);span;span>>=1) width++;
  size_t best=3+(n*width+7)/8;
  *mode=16+width;
  *minv=(int16_t)lo;
  for(unsigned k=0;k<16;k++) {
    size_t bits=n*(k+1);
    for(size_t i=0;i<n;i++) bits+=leopard_short_zigzag(v[i])>>k;
    size_t bytes=1+(bits+7)/8;
    if(bytes<best) { best=bytes; *mode=k; }
  }
  return best;
}

static int leopard_short_encode(uint8_t *buf, size_t len, size_t *used,
                                 const int16_t *v, size_t n) {
  if(n==0 || n>512) return 1;
  unsigned mode;
  int16_t lo;
  size_t bytes=leopard_short_plan(v,n,&mode,&lo);
  if(bytes>len) return 1;
  memset(buf,0,bytes);
  buf[0]=(uint8_t)mode;
  size_t pos=8;
  unsigned width=mode;
  if(mode>=16) {
    buf[1]=(uint8_t)(uint16_t)lo;
    buf[2]=(uint8_t)((uint16_t)lo>>8);
    pos=24;
    width=mode-16;
  }
  for(size_t i=0;i<n;i++) {
    uint32_t x=mode<16 ? leopard_short_zigzag(v[i]) : (uint32_t)((int32_t)v[i]-lo);
    if(mode<16) {
      pos+=x>>width;
      buf[pos>>3]|=(uint8_t)(1u<<(pos&7));
      pos++;
    }
    for(unsigned j=0;j<width;j++,pos++)
      buf[pos>>3]|=(uint8_t)(((x>>j)&1u)<<(pos&7));
  }
  *used=bytes;
  return 0;
}

static int leopard_short_decode(int16_t *v, size_t n, const uint8_t *buf,
                                 size_t len, size_t *used) {
  if(n==0 || n>512 || len==0 || len>SIZE_MAX/8) return 1;
  unsigned mode=buf[0], width=mode;
  size_t pos=8;
  int32_t lo=0;
  if(mode>32) return 1;
  if(mode>=16) {
    if(len<3) return 1;
    uint32_t u=(uint32_t)buf[1]|((uint32_t)buf[2]<<8);
    lo=u<32768 ? (int32_t)u : (int32_t)u-65536;
    pos=24;
    width=mode-16;
  }
  for(size_t i=0;i<n;i++) {
    uint32_t x=0;
    if(mode<16) {
      uint32_t q=0;
      for(;;pos++) {
        if(pos>=len*8) return 1;
        if((buf[pos>>3]>>(pos&7))&1u) { pos++; break; }
        if(++q>(65535u>>width)) return 1;
      }
      x=q<<width;
    }
    if(width>len*8-pos) return 1;
    for(unsigned j=0;j<width;j++,pos++)
      x|=(uint32_t)((buf[pos>>3]>>(pos&7))&1u)<<j;
    int32_t value=mode<16 ? ((x&1u) ? -(int32_t)(x>>1)-1 : (int32_t)(x>>1)) : lo+(int32_t)x;
    if(value < -32768 || value > 32767) return 1;
    v[i]=(int16_t)value;
  }
  size_t bytes=(pos+7)/8;
  if((pos&7) && (buf[bytes-1]>>(pos&7))) return 1;
  unsigned canonical_mode;
  int16_t canonical_lo;
  if(leopard_short_plan(v,n,&canonical_mode,&canonical_lo)!=bytes ||
     canonical_mode!=mode || (mode>=16 && canonical_lo!=lo)) return 1;
  *used=bytes;
  return 0;
}
#endif
