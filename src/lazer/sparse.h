#ifdef _OPENMP
#include "threading.h"
#endif

/* Sparse indices are assembled serially; each polynomial then has one owner. */
static int affine_sparse_apply(spolymat_t dst, spolymat_t src,
                              int_srcptr scalar, unsigned int rotation, int op) {
  if(dst==src || src->nelems<256 || getenv("LNP_DISABLE_AFFINE_SPARSE")) return 0;
  dst->nelems=0;
  for(unsigned int j=0;j<src->nelems;j++)
    spolymat_insert_elem(dst,spolymat_get_row(src,j),spolymat_get_col(src,j));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int j=0;j<src->nelems;j++) {
    poly_ptr a=spolymat_get_elem(dst,j), b=spolymat_get_elem(src,j);
    if(op==0) poly_set(a,b);
    else if(op==1) poly_mod(a,b);
    else if(op==2) poly_scale(a,scalar,b);
    else poly_lrot(a,b,rotation);
  }
  dst->sorted=src->sorted;
  return 1;
}

static int affine_sparse_add(spolymat_t dst, spolymat_t a, spolymat_t b, int crt) {
  if(dst==a || dst==b || a->nelems+b->nelems<512 || getenv("LNP_DISABLE_AFFINE_SPARSE")) return 0;
  size_t cap=a->nelems+b->nelems;
  poly_ptr *av=calloc(cap,sizeof(*av)), *bv=calloc(cap,sizeof(*bv));
  if(!av || !bv) abort();
  unsigned int i=0,j=0,k=0;
  dst->nelems=0;
  while(i<a->nelems || j<b->nelems) {
    uint32_t ai=i<a->nelems ? ((uint32_t)spolymat_get_row(a,i)<<16)|spolymat_get_col(a,i) : UINT32_MAX;
    uint32_t bj=j<b->nelems ? ((uint32_t)spolymat_get_row(b,j)<<16)|spolymat_get_col(b,j) : UINT32_MAX;
    uint32_t key=ai<bj ? ai : bj;
    spolymat_insert_elem(dst,key>>16,key&65535);
    if(ai<=bj) av[k]=spolymat_get_elem(a,i++);
    if(bj<=ai) bv[k]=spolymat_get_elem(b,j++);
    k++;
  }
  /* Inputs can be the same matrix: pre-convert once, before lazy writes. */
  if(crt) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
    for(unsigned int x=0;x<a->nelems;x++) _tocrt(spolymat_get_elem(a,x));
    if(a!=b) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
      for(unsigned int x=0;x<b->nelems;x++) _tocrt(spolymat_get_elem(b,x));
    }
  } else {
    spolymat_fromcrt(a);
    if(a!=b) spolymat_fromcrt(b);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int x=0;x<k;x++) {
    poly_ptr out=spolymat_get_elem(dst,x);
    if(av[x] && bv[x]) poly_add(out,av[x],bv[x],crt);
    else poly_set(out,av[x] ? av[x] : bv[x]);
  }
  dst->sorted=1;
  free(bv); free(av);
  return 1;
}

static int affine_sparse_scale2(spolymat_t dst, poly_t scalar, spolymat_t src) {
  if(dst==src || src->nelems<256 || getenv("LNP_DISABLE_AFFINE_SPARSE")) return 0;
  _tocrt(scalar);
  dst->nelems=0;
  for(unsigned int j=0;j<src->nelems;j++)
    spolymat_insert_elem(dst,spolymat_get_row(src,j),spolymat_get_col(src,j));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int j=0;j<src->nelems;j++)
    poly_mul(spolymat_get_elem(dst,j),scalar,spolymat_get_elem(src,j));
  dst->sorted=src->sorted;
  return 1;
}
