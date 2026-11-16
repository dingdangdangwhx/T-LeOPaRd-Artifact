#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "malloc.h"
#include "fips202.h"
#include "aesctr.h"
#include "data.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "sparsemat.h"
#include "jlproj.h"
#include "labrador.h"
#include "chihuahua.h"
#ifndef CHIHUAHUA_VERIFY_PROFILE
#define CHIHUAHUA_VERIFY_PROFILE 0
#endif
#ifndef CHIHUAHUA_FUSED_COLLAPSE
#define CHIHUAHUA_FUSED_COLLAPSE 1
#endif

static double chihuahua_elapsed_sec(const struct timespec *a, const struct timespec *b) {
  return (double)(b->tv_sec - a->tv_sec) + 1e-9*(double)(b->tv_nsec - a->tv_nsec);
}


static labrador_scalar bytes_to_centered_modq(const uint8_t *buf) {
  size_t i;
  unsigned __int128 x = 0;
  const unsigned __int128 q = ((unsigned __int128)1 << LOGQ) - QOFF;
  unsigned __int128 y;

  for(i=0;i<QBYTES;i++)
    x |= (unsigned __int128)buf[i] << (8*i);
  y = x;
  if(y >= q)
    y -= q;
  if(y > q/2)
    return (labrador_scalar)((__int128)y - (__int128)q);
  return (labrador_scalar)y;
}

polx *init_sparsecnst_half(sparsecnst *cnst, size_t r, size_t nz, size_t buflen, size_t deg,
                           int quadratic, int homogeneous)
{
  void *buf;

  cnst->deg = deg;
  cnst->nz = nz;
  cnst->scale = NULL;
  cnst->a->len = 0;
  if(quadratic) {
    buf = _malloc((r*r+r)*sizeof(size_t));
    cnst->a->rows = (size_t*)buf;
    cnst->a->cols = (size_t*)buf + (r*r+r)/2;

    buf = _aligned_alloc(64,(r*r+r)/2*sizeof(polx));
    cnst->a->coeffs = (polx*)buf;
  }
  else {
    cnst->a->rows = NULL;
    cnst->a->cols = NULL;
    cnst->a->coeffs = NULL;
  }

  if(!homogeneous)
    buflen += MAX(1,deg);

  buf = _malloc(4*nz*sizeof(size_t) + nz*sizeof(polx*) + 63 + buflen*sizeof(polx));
  cnst->idx = (size_t*)buf;
  cnst->off = &cnst->idx[nz];
  cnst->len = &cnst->off[nz];
  cnst->mult = &cnst->len[nz];
  cnst->phi = (polx**)&cnst->mult[nz];
  buf = (void*)(((uintptr_t)&cnst->phi[nz] + 63) & -64ULL);

  if(homogeneous)
    cnst->b = NULL;
  else {
    cnst->b = buf;
    buf = (polx*)buf + MAX(1,deg);
  }

  return buf;
}

void init_sparsecnst_view(sparsecnst *cnst, size_t r, size_t nz, size_t deg,
                          int quadratic, int homogeneous) {
  size_t i;

  (void)init_sparsecnst_half(cnst,r,nz,0,deg,quadratic,homogeneous);
  cnst->scale = _malloc(nz*sizeof(*cnst->scale));
  for(i=0;i<nz;i++)
    cnst->scale[i] = 1;
}

void init_sparsecnst_raw(sparsecnst *cnst, size_t r, size_t nz, const size_t idx[nz], const size_t n[nz], size_t deg,
                         int quadratic, int homogeneous)
{
  size_t i;
  size_t buflen;
  size_t elen[nz];
  polx *buf;

  buflen = 0;
  for(i=0;i<nz;i++) {
    elen[i] = extlen(n[i]*deg,deg);
    buflen += elen[i];
  }

  buf = init_sparsecnst_half(cnst,r,nz,buflen,deg,quadratic,homogeneous);
  for(i=0;i<nz;i++) {
    cnst->idx[i] = idx[i];
    cnst->off[i] = 0;
    cnst->len[i] = n[i]*deg;
    cnst->mult[i] = 1;
    cnst->phi[i] = (polx*)buf;
    buf += elen[i];
  }
}

int set_sparsecnst_raw(sparsecnst *cnst, uint8_t h[16], size_t nz, const size_t idx[nz],
                       const size_t n[nz], size_t deg, int64_t *phi, int64_t *b)
{
  size_t j,k;
  polz t[deg];
  __attribute__((aligned(16)))
  uint8_t hashbuf[deg*N*QBYTES];
  shake128incctx shakectx;

  if(nz != cnst->nz) {
    fprintf(stderr,"ERROR in set_sparsecnst_raw(): Mismatch in number of affected vectors\n");
    return 1;
  }
  if(deg != cnst->deg) {
    fprintf(stderr,"ERROR in set_sparsecnst_raw(): Mismatch in extension degree\n");
    return 2;
  }
  if((b && !cnst->b) || (!b && cnst->b)) {
    fprintf(stderr,"ERROR in set_sparsecnst_raw(): Mismatch in homogeneity\n");
    return 3;
  }
  for(j=0;j<nz;j++) {
    if(idx[j] != cnst->idx[j]) {
      fprintf(stderr,"ERROR in set_sparsecnst_raw(): Mismatch in index of %zu-th affected vector\n",j);
      return 4;
    }
    if(cnst->off[j])  //FIXME
      return 5;
    if(cnst->mult[j] != 1)
      return 6;
  }

  cnst->a->len = 0;

  shake128_inc_init(&shakectx);
  shake128_inc_absorb(&shakectx,h,16);

  if(b) {
    polzvec_fromint64vec(t,1,deg,b);
    polzvec_topolxvec(cnst->b,t,deg);
    polzvec_bitpack(hashbuf,t,deg);
    shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
  }

  for(j=0;j<nz;j++) {
    for(k=0;k<n[j];k++) {
      polzvec_fromint64vec(t,1,deg,phi);
      polzvec_topolxvec(&cnst->phi[j][k*deg],t,deg);
      polzvec_bitpack(hashbuf,t,deg);
      shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
      phi += deg*N;
    }
  }

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(h,16,&shakectx);
  return 0;
}

void free_sparsecnst(sparsecnst *cnst) {
  free(cnst->idx);
  free(cnst->a->rows);
  free(cnst->a->coeffs);
  free(cnst->scale);
  cnst->idx = cnst->a->rows = NULL;
  cnst->a->coeffs = NULL;
  cnst->scale = NULL;
}

static int64_t sparsecnst_scale(const sparsecnst *cnst, size_t i) {
  return cnst->scale ? cnst->scale[i] : 1;
}

static labrador_scalar centered_scalar_mul(labrador_scalar a, int64_t b) {
  return scalar_mulq(a,b);
}

/* Multiplying two full-width q residues directly can exceed the CRT
 * reconstruction headroom at LOGQ=64.  Evaluate the scalar in small signed
 * limbs and canonicalize after every Horner step instead. */
static void polxvec_full_scalar_scale_add(polx *dst, const polx *src,
                                          size_t len, labrador_scalar scalar) {
#if LOGQ == 72 && !defined(CHIHUAHUA_SCALAR_HORNER)
  polxvec_scale_add_modq(dst,src,len,scalar);
#else
  enum { LIMB_BITS = 22 };
  const int64_t base = (int64_t)1 << LIMB_BITS;
  int64_t digits[(LOGQ + LIMB_BITS - 1)/LIMB_BITS + 1];
  __int128 x = scalar;
  size_t ndigits = 0, i;
  polx *tmp;

  do {
    digits[ndigits++] = (int64_t)(x % base);
    x /= base;
  } while(x != 0);

  tmp = _aligned_alloc(64,len*sizeof(*tmp));
  polxvec_setzero(tmp,len);
  for(i=ndigits;i>0;i--) {
    if(i != ndigits) {
      polxvec_scale(tmp,tmp,len,base);
      polxvec_refresh(tmp,len);
    }
    if(digits[i-1] != 0) {
      polxvec_scale_add(tmp,src,len,digits[i-1]);
      polxvec_refresh(tmp,len);
    }
  }
  polxvec_add(dst,dst,tmp,len);
  polxvec_refresh(dst,len);
  free(tmp);
#endif
}

static void free_sparsecnst_array(sparsecnst *cnst, size_t k) {
  size_t i;

  if(!cnst)
    return;
  for(i=0;i<k;i++)
    free_sparsecnst(&cnst[i]);
  free(cnst);
}

static void refresh_sparsecnst_touched_phi(const proof *pi,
                                           constraint *ocnst,
                                           const sparsecnst *icnst,
                                           size_t k, size_t n) {
  size_t i,j,u,v;

  for(i=0;i<k;i++) {
    if(icnst[i].deg != 0)
      continue;
    for(j=0;j<icnst[i].nz;j++) {
      polx (*phi)[n] = (polx(*)[n])ocnst->phi;

      u = v = 0;
      while(u < icnst[i].idx[j]) {
        phi += pi->nu[u];
        v = pi->nu[u] ? 0 : v + pi->n[u];
        u += 1;
      }
      polxvec_refresh(&(*phi)[v + icnst[i].off[j]],icnst[i].len[j]);
    }
  }
}

static size_t proof_rowbase(const proof *pi, size_t idx) {
  size_t i, base = 0;

  for(i=0;i<idx;i++)
    base += pi->nu[i];
  return base;
}

static void sparsemat_add_coeff(sparsemat *a, size_t row, size_t col, const polx *coeff) {
  size_t i;

  if(row > col) {
    size_t tmp = row;
    row = col;
    col = tmp;
  }

  for(i=0;i<a->len;i++) {
    if(a->rows[i] == row && a->cols[i] == col) {
      polx_add(&a->coeffs[i],&a->coeffs[i],coeff);
      return;
    }
  }

  a->rows[a->len] = row;
  a->cols[a->len] = col;
  polxvec_copy(&a->coeffs[a->len],coeff,1);
  a->len += 1;
}

static void sparsemat_scale_add_rankchunks(sparsemat *dst, const sparsemat *src, const proof *pi, labrador_scalar s) {
  size_t i,j;
  polx coeff[1];

  for(i=0;i<src->len;i++) {
    size_t row = src->rows[i];
    size_t col = src->cols[i];
    size_t base_row = proof_rowbase(pi,row);
    size_t base_col = proof_rowbase(pi,col);
    size_t chunks = MIN(pi->nu[row],pi->nu[col]);

    for(j=0;j<chunks;j++) {
#if LOGQ > 64
      polxvec_setzero(coeff,1);
      polxvec_full_scalar_scale_add(coeff,&src->coeffs[i],1,s);
#else
      polx_scale(coeff,&src->coeffs[i],s);
#endif
      sparsemat_add_coeff(dst,base_row+j,base_col+j,coeff);
    }
  }
}

static void sparsemat_polx_mul_add_rankchunks(sparsemat *dst, const sparsemat *src, const proof *pi, const polx *c) {
  size_t i,j;
  polx coeff[1];

  for(i=0;i<src->len;i++) {
    size_t row = src->rows[i];
    size_t col = src->cols[i];
    size_t base_row = proof_rowbase(pi,row);
    size_t base_col = proof_rowbase(pi,col);
    size_t chunks = MIN(pi->nu[row],pi->nu[col]);

    for(j=0;j<chunks;j++) {
      polx_mul(coeff,c,&src->coeffs[i]);
      sparsemat_add_coeff(dst,base_row+j,base_col+j,coeff);
    }
  }
}

void sparsecnst_eval(polx *b, const sparsecnst *cnst, polx *sx[], const witness *wt) {
  const size_t r = wt->r;
  const size_t *n = wt->n;

  size_t i,j,k;
  size_t deg2 = MAX(1,cnst->deg);
  polx t[deg2];

  if(!*sx) {
    k = 0;
    for(i=0;i<r;i++)
      k += n[i];
    sx[0] = _aligned_alloc(64,k*sizeof(polx));
    for(i=1;i<r;i++)
      sx[i] = &sx[i-1][n[i-1]];
    for(i=0;i<r;i++)
      polxvec_frompolyvec(sx[i],wt->s[i],n[i]);
  }

  polxvec_setzero(b,deg2);

  /* quadratic term */
  for(i=0;i<cnst->a->len;i++) {
    j = cnst->a->rows[i];
    k = cnst->a->cols[i];
    polxvec_sprod(t,sx[j],sx[k],MIN(n[j],n[k]));  // TODO: store?
    polx_mul(t,t,&cnst->a->coeffs[i]);
    polx_add(b,b,t);
  }

  /* linear term */
  for(i=0;i<cnst->nz;i++) {
    j = cnst->idx[i];
    k = cnst->off[i];
    polxvec_mul_extension(t,cnst->phi[i],&sx[j][k],cnst->len[i],deg2/cnst->mult[i],cnst->mult[i]);
    if(sparsecnst_scale(cnst,i) != 1)
      polxvec_scale(t,t,deg2,sparsecnst_scale(cnst,i));
    polxvec_add(b,b,t,deg2);
  }
}

int sparsecnst_check(const sparsecnst *cnst, polx *sx[], const witness *wt) {
  int ret;
  size_t deg2 = MAX(1,cnst->deg);
  polx b[deg2];

  sparsecnst_eval(b,cnst,sx,wt);
  if(cnst->b)
    polxvec_sub(b,b,cnst->b,deg2);

  if(cnst->deg)
    ret = polxvec_iszero(b,deg2);
  else
    ret = polx_iszero_constcoeff(b);

  return ret;
}

int init_prncplstmnt_raw(prncplstmnt *st, size_t r, const size_t n[r],
                         uint64_t betasq, size_t k, int quadratic)
{
  size_t i;
  void *buf;

  if(betasq > JLMAXNORMSQ) {
    fprintf(stderr,"ERROR in init_prcplstmnt_raw(): Total witness norm too big for JL projection\n");
    return 1;
  }

  buf = _malloc(r*sizeof(size_t) + k*sizeof(sparsecnst));
  st->n = (size_t*)buf;
  st->cnst = (sparsecnst*)&st->n[r];

  st->r = r;
  st->k = k;
  st->quadratic = quadratic;
  st->betasq = betasq;
  for(i=0;i<r;i++)
    st->n[i] = n[i];
  for(i=0;i<k;i++)
    st->cnst[i].idx = NULL;

  st->frontend_data = NULL;
  st->free_frontend_data = NULL;
  st->build_postcommit_cnst = NULL;

  memset(st->h,0,16);
  return 0;
}

int set_prncplstmnt_lincnst_raw(prncplstmnt *st, size_t i, size_t nz, const size_t idx[nz],
                                const size_t n[nz], size_t deg, int64_t *phi, int64_t *b)
{
  size_t j,k;

  if(i >= st->k) {
    fprintf(stderr,"ERROR in set_prcplstmnt_lincnst_raw(): Constraint %zu does not exist\n",i);
    return 1;
  }

  sparsecnst *cnst = &st->cnst[i];
  if(cnst->idx) {
    fprintf(stderr,"ERROR in set_prncplstmnt_lincnst_raw(): Constraint has already been set\n");
    return 2;
  }
  for(j=0;j<nz;j++) {
    k = idx[j];
    if(k >= st->r) {
      fprintf(stderr,"ERROR in set_prncplstmnt_lincnst_raw(): Witness vector %zu does not exist\n",k);
      return 3;
    }
    if(n[j]*deg != st->n[k]) {
      fprintf(stderr,"ERROR in set_prncplstmnt_lincnst_raw(): Mismatch in witness vector length\n");
      return 4;
    }
  }

  init_sparsecnst_raw(cnst,st->r,nz,idx,n,deg,0,b == NULL);
  (void)set_sparsecnst_raw(cnst,st->h,nz,idx,n,deg,phi,b);  // errors can not happen
  return 0;
}

void free_prncplstmnt(prncplstmnt *st) {
  size_t i;

  if(!st->n) return;
  for(i=0;i<st->k;i++)
    free_sparsecnst(&st->cnst[i]);
  if(st->free_frontend_data)
    st->free_frontend_data(st->frontend_data);
  free(st->n);  // one buffer for n,cnst
  st->n = NULL;
  st->frontend_data = NULL;
  st->free_frontend_data = NULL;
}

void print_prncplstmnt_pp(const prncplstmnt *st) {
  size_t i;

  printf("Chihuahua statement ");
  for(i=0;i<5;i++)
    printf("%02hhX",st->h[i]);
  printf(":\n");

  printf("  Witness multiplicity: %zu\n",st->r);
  printf("  Witness ranks: ");
  for(i=0;i<st->r;i++) {
    printf("%zu",st->n[i]);
    if(i<st->r-1) printf(", ");
  }
  printf("\n");

  printf("  Number of dot-product constraints: %zu\n",st->k);
  printf("  Norm constraint: %.2f\n",sqrt(st->betasq));
  printf("\n");
}

void collaps_sparsecnst(constraint *ocnst, statement *ost, const proof *pi, const sparsecnst *icnst, size_t k) {
  size_t i,j,u,v;
  size_t challenge_idx = 0;
  const size_t n = ost->n;
  labrador_scalar s;
  polx (*phi)[n];

  j = 0;
  for(i=0;i<k;i++)
    if(icnst[i].deg == 0)
      j += 1;

  __attribute__((aligned(16)))
  uint8_t hashbuf[16+j*QBYTES];

  shake128(hashbuf,sizeof(hashbuf),ost->h,16);
  memcpy(ost->h,hashbuf,16);

  for(i=0;i<k;i++,icnst++) {
    if(icnst->deg) continue;

    s = bytes_to_centered_modq(&hashbuf[16+challenge_idx++*QBYTES]);

    if(icnst->b)
      polxvec_full_scalar_scale_add(ocnst->b,icnst->b,1,s);

    for(j=0;j<icnst->nz;j++) {
      phi = (polx(*)[n])ocnst->phi;
      u = v = 0;
      while(u < icnst->idx[j]) {
        phi += pi->nu[u];
        v = (pi->nu[u]) ? 0 : v + pi->n[u];
        u += 1;
      }
      polxvec_full_scalar_scale_add(
          &(*phi)[v+icnst->off[j]],icnst->phi[j],icnst->len[j],
          centered_scalar_mul(s,sparsecnst_scale(icnst,j)));
    }

    sparsemat_scale_add_rankchunks(ocnst->a,icnst->a,pi,s);
  }
}

/*
 * Check a single homogeneous equation directly, without a possibly zero scale.
 * Preserve transcript evolution but avoid a redundant full-width scale.
 */
static void collaps_single_homogeneous_sparsecnst(
    constraint *ocnst, statement *ost, const proof *pi,
    const sparsecnst *icnst) {
  uint8_t hashbuf[16+QBYTES];
  const size_t n = ost->n;
  size_t j, u, v;
  polx (*phi)[n];

  shake128(hashbuf,sizeof(hashbuf),ost->h,16);
  memcpy(ost->h,hashbuf,16);
  for(j=0;j<icnst->nz;j++) {
    phi = (polx(*)[n])ocnst->phi;
    u = v = 0;
    while(u < icnst->idx[j]) {
      phi += pi->nu[u];
      v = pi->nu[u] ? 0 : v + pi->n[u];
      u += 1;
    }
    polxvec_add(&(*phi)[v+icnst->off[j]],
                &(*phi)[v+icnst->off[j]],
                icnst->phi[j],icnst->len[j]);
  }
}


static void reduce_phi_lanes(polx *phi, size_t len) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(len >= 4096)
#endif
  for(size_t off=0;off<len;off+=256)
    polxvec_reduce(&phi[off],MIN(256,len-off));
}

typedef struct {
  size_t offset, len;
  const polx *coeff;
} aggregate_term;

/* Partition by destination, not by constraint: overlapping contributions
 * retain their original order, while independent output chunks run in parallel. */
static int aggregate_degree_one(statement *ost, const proof *pi,
                                const sparsecnst *cnst, size_t k,
                                const uint8_t seed[16]) {
  enum { BLOCK = 128 };
  size_t count=0, t=0;
  const size_t length=ost->r*ost->n, chunks=(length+BLOCK-1)/BLOCK;
  if(length<4096) return 0;
  for(size_t i=0;i<k;i++) {
    if(cnst[i].deg==0) continue;
    if(cnst[i].deg!=1) return 0;
    for(size_t j=0;j<cnst[i].nz;j++)
      if(cnst[i].mult[j]!=1) return 0;
    count+=cnst[i].nz;
  }
  if(count<8) return 0;
  aggregate_term *terms=_malloc(count*sizeof(*terms));
  polx *factors=_aligned_alloc(64,count*sizeof(*factors));
  size_t *offsets=_malloc((chunks+1)*sizeof(*offsets));
  memset(offsets,0,(chunks+1)*sizeof(*offsets));
  for(size_t i=0;i<k;i++) {
    if(cnst[i].deg==0) continue;
    polx alpha;
    polxvec_quarternary(&alpha,1,seed,i);
    if(cnst[i].b) polxvec_sprod_add(ost->cnst->b,&alpha,cnst[i].b,1);
    for(size_t j=0;j<cnst[i].nz;j++,t++) {
      size_t base=0, v=0;
      for(size_t u=0;u<cnst[i].idx[j];u++) {
        base+=pi->nu[u]*ost->n;
        v=pi->nu[u] ? 0 : v+pi->n[u];
      }
      terms[t]=(aggregate_term){base+v+cnst[i].off[j],cnst[i].len[j],cnst[i].phi[j]};
      const int64_t scale=sparsecnst_scale(&cnst[i],j);
      if(scale==1) factors[t]=alpha;
      else polx_scale(&factors[t],&alpha,scale);
      if(terms[t].len)
        for(size_t b=terms[t].offset/BLOCK;b<=(terms[t].offset+terms[t].len-1)/BLOCK;b++)
          offsets[b+1]++;
    }
    sparsemat_polx_mul_add_rankchunks(ost->cnst->a,cnst[i].a,pi,&alpha);
  }
  for(size_t b=0;b<chunks;b++) offsets[b+1]+=offsets[b];
  size_t *entries=_malloc(offsets[chunks]*sizeof(*entries));
  size_t *cursor=_malloc(chunks*sizeof(*cursor));
  memcpy(cursor,offsets,chunks*sizeof(*cursor));
  for(t=0;t<count;t++)
    if(terms[t].len)
      for(size_t b=terms[t].offset/BLOCK;b<=(terms[t].offset+terms[t].len-1)/BLOCK;b++)
        entries[cursor[b]++]=t;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for(size_t b=0;b<chunks;b++) {
    const size_t first=b*BLOCK, last=MIN(first+BLOCK,length);
    for(size_t j=offsets[b];j<offsets[b+1];j++) {
      const size_t id=entries[j];
      const size_t lo=MAX(first,terms[id].offset), hi=MIN(last,terms[id].offset+terms[id].len);
      polxvec_polx_mul_add(ost->cnst->phi+lo,factors+id,
                          terms[id].coeff+lo-terms[id].offset,hi-lo);
    }
  }
  free(cursor); free(entries); free(offsets); free(factors); free(terms);
  return 1;
}

void aggregate_sparsecnst(statement *ost, const proof *pi, const sparsecnst *cnst, size_t k) {
  size_t i,j,u,v;
  __attribute__((aligned(16)))
  uint8_t hashbuf[32];
  polx (*phi)[ost->n];

  shake128(hashbuf,32,ost->h,16);
  memcpy(ost->h,hashbuf,16);
  reduce_phi_lanes(ost->cnst->phi,ost->r*ost->n);
#ifndef CHIHUAHUA_SERIAL_AGGREGATE
  if(aggregate_degree_one(ost,pi,cnst,k,&hashbuf[16])) {
    ost->cnst->a->coeffs=realloc(ost->cnst->a->coeffs,ost->cnst->a->len*sizeof(polx));
    return;
  }
#endif

  for(i=0;i<k;i++,cnst++) {
    if(cnst->deg == 0) continue;
    polx alpha[cnst->deg];
    polxvec_quarternary(alpha,cnst->deg,&hashbuf[16],i);

    if(cnst->b)
      polxvec_sprod_add(ost->cnst->b,alpha,cnst->b,cnst->deg);

    for(j=0;j<cnst->nz;j++) {
      const int64_t scale = sparsecnst_scale(cnst,j);
      const polx *collapse_alpha = alpha;
      polx scaled_alpha[cnst->deg];

      phi = (polx(*)[ost->n])ost->cnst->phi;
      u = v = 0;
      while(u < cnst->idx[j]) {
        phi += pi->nu[u];
        v = (pi->nu[u]) ? 0 : v + pi->n[u];
        u += 1;
      }
      if(scale != 1) {
        polxvec_scale(scaled_alpha,alpha,cnst->deg,scale);
        collapse_alpha = scaled_alpha;
      }
      if(CHIHUAHUA_FUSED_COLLAPSE && cnst->deg == 1 && cnst->mult[j] == 1) {
        const size_t len = cnst->len[j];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(len >= 4096)
#endif
        for(size_t off=0;off<len;off+=256)
          polxvec_polx_mul_add(&(*phi)[v+cnst->off[j]+off],collapse_alpha,
                               &cnst->phi[j][off],MIN(256,len-off));
      } else {
        polxvec_collaps_add_extension(&(*phi)[v+cnst->off[j]],collapse_alpha,cnst->phi[j],cnst->len[j],
                                      cnst->deg/cnst->mult[j],cnst->mult[j]);
        /* Pointwise-add kernels leave unreduced signed 16-bit lanes. */
        reduce_phi_lanes(&(*phi)[v+cnst->off[j]],cnst->len[j]);
      }
    }

    sparsemat_polx_mul_add_rankchunks(ost->cnst->a,cnst->a,pi,alpha);
  }

  ost->cnst->a->coeffs = realloc(ost->cnst->a->coeffs,ost->cnst->a->len*sizeof(polx));
}

static int has_coefficient_constraint(const sparsecnst *cnst, size_t k) {
  for(size_t i=0;i<k;i++)
    if(cnst[i].deg == 0) return 1;
  return 0;
}

int principle_prove(statement *ost, witness *owt, proof *pi, const prncplstmnt *ist, const witness *iwt, int tail) {
  int ret;
  size_t i;
  constraint cnst[1] = {};
  sparsecnst *postcommit = NULL;
  size_t postcommit_k = 0;
  void *buf = NULL;
  struct timespec t0, t1;
  double t_commit = 0.0, t_project = 0.0, t_lift = 0.0, t_aggregate = 0.0, t_amortize = 0.0;
  double t_lift_collapse = 0.0, t_lift_eval = 0.0, t_refresh = 0.0;

  ret = init_proof(pi,iwt,ist->quadratic,tail);
  if(ret)  // commitments not secure (1/2)
    return ret;
  init_statement(ost,pi,ist->h);
  init_witness(owt,ost);
  printf("Predicted witness norm: %.2f\n\n",sqrt(pi->normsq));

  {
    buf = _aligned_alloc(64,ost->r*ost->n*(sizeof(polx)+JL_BYTES));
    polx (*sx)[ost->n] = (polx(*)[ost->n])buf;
    uint8_t (*jlmat)[ost->n][JL_BYTES] = (uint8_t(*)[ost->n][JL_BYTES])sx[ost->r];
    clock_gettime(CLOCK_MONOTONIC,&t0);
    commit(ost,owt,pi,sx,iwt);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_commit = chihuahua_elapsed_sec(&t0,&t1);
    if(ist->build_postcommit_cnst) {
      ret = ist->build_postcommit_cnst(&postcommit,&postcommit_k,ist,ost->h);
      if(ret) {
        ret += 100;
        goto err;
      }
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    ret = project(ost,pi,jlmat,iwt);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_project = chihuahua_elapsed_sec(&t0,&t1);
    if(ret) {
      ret += 10;
      goto err;
    }

    init_constraint(cnst,ost);
    for(i=0;i<LIFTS;i++) {
      clock_gettime(CLOCK_MONOTONIC,&t0);
      collaps_jlproj(cnst,ost,pi,jlmat);
      if(has_coefficient_constraint(postcommit,postcommit_k)) {
        if(postcommit_k == 1 && postcommit[0].deg == 0 &&
           postcommit[0].b == NULL && postcommit[0].a->len == 0 &&
           postcommit[0].scale == NULL)
          collaps_single_homogeneous_sparsecnst(
              cnst,ost,pi,&postcommit[0]);
        else
          collaps_sparsecnst(cnst,ost,pi,postcommit,postcommit_k);
        refresh_sparsecnst_touched_phi(pi,cnst,postcommit,postcommit_k,ost->n);
      }
      collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      t_lift_collapse += chihuahua_elapsed_sec(&t0,&t1);
      clock_gettime(CLOCK_MONOTONIC,&t0);
      lift_aggregate_zqcnst(ost,pi,i,cnst,sx);
      clock_gettime(CLOCK_MONOTONIC,&t1);
      t_lift_eval += chihuahua_elapsed_sec(&t0,&t1);
    }
    t_lift = t_lift_collapse + t_lift_eval;
    free_constraint(cnst);
    if(has_coefficient_constraint(postcommit,postcommit_k)) {
      polx_refresh(ost->cnst->b);
      polxvec_refresh(ost->cnst->phi,ost->r*ost->n);
      polxvec_refresh(ost->cnst->a->coeffs,ost->cnst->a->len);
    }

    clock_gettime(CLOCK_MONOTONIC,&t0);
    aggregate_sparsecnst(ost,pi,ist->cnst,ist->k);
    if(postcommit_k)
      aggregate_sparsecnst(ost,pi,postcommit,postcommit_k);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_aggregate = chihuahua_elapsed_sec(&t0,&t1);
    clock_gettime(CLOCK_MONOTONIC,&t0);
    amortize(ost,owt,pi,sx);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_amortize = chihuahua_elapsed_sec(&t0,&t1);
    free(buf);
    buf = NULL;
  }

  clock_gettime(CLOCK_MONOTONIC,&t0);
  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_refresh = chihuahua_elapsed_sec(&t0,&t1);
  printf("Principle prove breakdown: commit=%.4fs, project=%.4fs, lift=%.4fs (collapse=%.4fs, eval=%.4fs), aggregate=%.4fs, amortize=%.4fs, refresh=%.4fs\n",
         t_commit, t_project, t_lift, t_lift_collapse, t_lift_eval, t_aggregate, t_amortize, t_refresh);
  free_sparsecnst_array(postcommit,postcommit_k);
  return 0;

err:
  free(buf);
  free_constraint(cnst);
  free_proof(pi);
  free_statement(ost);
  free_witness(owt);
  free_sparsecnst_array(postcommit,postcommit_k);
  return ret;
}

int principle_reduce(statement *ost, const proof *pi, const prncplstmnt *ist) {
  size_t i;
  int ret;
  uint8_t (*jlmat)[ost->n][JL_BYTES];
  constraint cnst[1] = {};
  sparsecnst *postcommit = NULL;
  size_t postcommit_k = 0;
#if CHIHUAHUA_VERIFY_PROFILE
  struct timespec t0, t1;
  double t_commit = 0.0, t_project = 0.0, t_init = 0.0;
  double t_jl = 0.0, t_sparse = 0.0, t_lift = 0.0;
  double t_agg = 0.0, t_amort = 0.0, t_refresh = 0.0;
#endif

  init_statement(ost,pi,ist->h);
  jlmat = _aligned_alloc(64,ost->r*ost->n*JL_BYTES);

#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  reduce_commit(ost,pi);
  if(ist->build_postcommit_cnst) {
    ret = ist->build_postcommit_cnst(&postcommit,&postcommit_k,ist,ost->h);
    if(ret) {
      ret += 100;
      goto err;
    }
  }
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_commit = chihuahua_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  ret = reduce_project(ost,jlmat,pi,pi->r,ist->betasq);
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_project = chihuahua_elapsed_sec(&t0,&t1);
#endif
  if(ret) goto err;  // projection too long

#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  init_constraint(cnst,ost);
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_init = chihuahua_elapsed_sec(&t0,&t1);
#endif
  for(i=0;i<LIFTS;i++) {
#if CHIHUAHUA_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
    collaps_jlproj(cnst,ost,pi,jlmat);
#if CHIHUAHUA_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_jl += chihuahua_elapsed_sec(&t0,&t1);
    clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
    if(has_coefficient_constraint(postcommit,postcommit_k)) {
      if(postcommit_k == 1 && postcommit[0].deg == 0 &&
         postcommit[0].b == NULL && postcommit[0].a->len == 0 &&
         postcommit[0].scale == NULL)
        collaps_single_homogeneous_sparsecnst(
            cnst,ost,pi,&postcommit[0]);
      else
        collaps_sparsecnst(cnst,ost,pi,postcommit,postcommit_k);
      refresh_sparsecnst_touched_phi(pi,cnst,postcommit,postcommit_k,ost->n);
    }
    collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
#if CHIHUAHUA_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_sparse += chihuahua_elapsed_sec(&t0,&t1);
    clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
    reduce_lift_aggregate_zqcnst(ost,pi,i,cnst);
#if CHIHUAHUA_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t1);
    t_lift += chihuahua_elapsed_sec(&t0,&t1);
#endif
  }
  free_constraint(cnst);
  free(jlmat);
  jlmat = NULL;
  if(has_coefficient_constraint(postcommit,postcommit_k)) {
    polx_refresh(ost->cnst->b);
    polxvec_refresh(ost->cnst->phi,ost->r*ost->n);
    polxvec_refresh(ost->cnst->a->coeffs,ost->cnst->a->len);
  }

#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  aggregate_sparsecnst(ost,pi,ist->cnst,ist->k);
  if(postcommit_k)
    aggregate_sparsecnst(ost,pi,postcommit,postcommit_k);
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_agg = chihuahua_elapsed_sec(&t0,&t1);
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  ret = reduce_amortize(ost,pi);
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_amort = chihuahua_elapsed_sec(&t0,&t1);
#endif
  if(ret) {  // commitments not secure (1/2)
    ret += 10;
    goto err;
  }

#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
#if CHIHUAHUA_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  t_refresh = chihuahua_elapsed_sec(&t0,&t1);
  printf("Chihuahua principle_reduce profile: commit=%.6fs project=%.6fs init=%.6fs jl=%.6fs sparse=%.6fs lift=%.6fs aggregate=%.6fs amortize=%.6fs refresh=%.6fs\n",
         t_commit,t_project,t_init,t_jl,t_sparse,t_lift,t_agg,t_amort,t_refresh);
#endif
  free_sparsecnst_array(postcommit,postcommit_k);
  return 0;

err:
  free_statement(ost);
  free(jlmat);
  free_constraint(cnst);
  free_sparsecnst_array(postcommit,postcommit_k);
  return ret;
}

int principle_verify(const prncplstmnt *st, const witness *wt) {
  int ret;
  size_t i;
  uint64_t normsq = 0;
  polx *sx[wt->r];

  if(wt->r != st->r) {
    fprintf(stderr,"ERROR in principle_verify(): Mismatch in witness multiplicity\n");
    return 1;
  }
  for(i=0;i<wt->r;i++)
    if(wt->n[i] != st->n[i]) {
      fprintf(stderr,"ERROR in principle_verify(): Mismatch in length of witness vector %zu\n",i);
      return 2;
    }
  for(i=0;i<wt->r;i++)
    normsq += polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  if(normsq > st->betasq) {
    fprintf(stderr,"ERROR in principle_verify(): Total witness vector norm too big\n");
    return 3;
  }

  *sx = NULL;
  for(i=0;i<st->k;i++) {
    ret = !sparsecnst_check(&st->cnst[i],sx,wt);
    if(ret) {
      fprintf(stderr,"ERROR in principle_verify(): Sparse dot-product constraint %zu does not hold\n",i);
      ret = 10+i;
      goto end;
    }
  }

end:
  free(*sx);
  return ret;
}
