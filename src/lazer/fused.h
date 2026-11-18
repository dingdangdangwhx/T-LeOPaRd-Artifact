#include "threading.h"

static int
nizk2_sparse_accumulate_inplace(spolymat_t dst, spolymat_t src, int_srcptr scale)
{
  const unsigned int nd=dst->nelems, ns=src->nelems;
  unsigned int i=0,j=0,missing=0;
  while(j<ns) {
    uint32_t a=i<nd ? ((uint32_t)spolymat_get_row(dst,i)<<16)|spolymat_get_col(dst,i) : UINT32_MAX;
    uint32_t b=((uint32_t)spolymat_get_row(src,j)<<16)|spolymat_get_col(src,j);
    if(a<b) i++;
    else { if(b<a) missing++; else i++; j++; }
  }
  if(nd+missing>dst->nelems_max) return 0;
  poly_ptr *target=calloc(ns ? ns : 1,sizeof(*target));
  if(!target) abort();
  spolymat_fromcrt(dst);
  spolymat_fromcrt(src);
  i=0; j=0;
  while(j<ns) {
    uint32_t a=i<nd ? ((uint32_t)spolymat_get_row(dst,i)<<16)|spolymat_get_col(dst,i) : UINT32_MAX;
    uint32_t b=((uint32_t)spolymat_get_row(src,j)<<16)|spolymat_get_col(src,j);
    if(a<b) { i++; continue; }
    if(a==b) target[j]=spolymat_get_elem(dst,i++);
    else {
      target[j]=spolymat_insert_elem(dst,b>>16,b&65535);
      poly_set_zero(target[j]);
    }
    j++;
  }
#ifdef _OPENMP
#pragma omp parallel if(nizk2_parallel_allowed())
#endif
  {
    poly_t product;
    poly_alloc(product,dst->ring);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int x=0;x<ns;x++) {
      poly_ptr add=spolymat_get_elem(src,x);
      if(scale) { poly_scale(product,scale,add); add=product; }
      poly_add(target[x],target[x],add,0);
    }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int x=0;x<dst->nelems;x++) {
      poly_ptr p=spolymat_get_elem(dst,x);
      poly_mod(p,p);
    }
    poly_free(product);
  }
  if(missing) spolymat_sort(dst);
  free(target);
  return 1;
}

/* One-source accumulation without three full sparse matrix temporaries. */
static int
nizk2_sparse_accumulate(spolymat_t dst, spolymat_t src, int_srcptr scale)
{
  if(!src || dst == src || dst->nelems + src->nelems < 256
     || getenv("NIZK2_DISABLE_FUSED")) return 0;
#if LAZER_AFFINE_FUSED >= 2
  return nizk2_sparse_accumulate_inplace(dst,src,scale);
#endif
  const size_t cap = dst->nelems + src->nelems;
  poly_ptr *old = calloc(cap,sizeof(*old)), *add = calloc(cap,sizeof(*add));
  if(!old || !add) abort();
  spolymat_t tmp;
  spolymat_alloc(tmp,dst->ring,dst->nrows,dst->ncols,cap);
  spolymat_fromcrt(dst);
  spolymat_fromcrt(src);
  unsigned int i=0,j=0,k=0;
  while(i<dst->nelems || j<src->nelems) {
    uint32_t a=i<dst->nelems ? ((uint32_t)spolymat_get_row(dst,i)<<16)|spolymat_get_col(dst,i) : UINT32_MAX;
    uint32_t b=j<src->nelems ? ((uint32_t)spolymat_get_row(src,j)<<16)|spolymat_get_col(src,j) : UINT32_MAX;
    uint32_t key=a<b ? a : b;
    spolymat_insert_elem(tmp,key>>16,key&65535);
    if(a<=b) old[k]=spolymat_get_elem(dst,i++);
    if(b<=a) add[k]=spolymat_get_elem(src,j++);
    k++;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int x=0;x<k;x++) {
    poly_ptr p=spolymat_get_elem(tmp,x);
    if(add[x]) {
      if(scale) poly_scale(p,scale,add[x]);
      else poly_set(p,add[x]);
      if(old[x]) poly_add(p,p,old[x],0);
    } else poly_set(p,old[x]);
    poly_mod(p,p);
  }
  tmp->sorted=1;
  spolymat_set(dst,tmp);
  spolymat_free(tmp);
  free(add); free(old);
  return 1;
}
