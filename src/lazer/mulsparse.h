#ifdef _OPENMP
#include "threading.h"
#endif
static int affine_sparse_matvec(polyvec_t out, spolymat_t matrix, polyvec_t vec) {
  if(out==vec || matrix->nelems<1024 || getenv("LNP_DISABLE_AFFINE_MATVEC")) return 0;
  unsigned int nr=out->nelems;
  unsigned int *starts=calloc(nr+1,sizeof(*starts)), *cursor=calloc(nr,sizeof(*cursor));
  unsigned int *indices=malloc(matrix->nelems*sizeof(*indices));
  if(!starts || !cursor || !indices) abort();
  for(unsigned int i=0;i<matrix->nelems;i++) starts[spolymat_get_row(matrix,i)+1]++;
  for(unsigned int r=0;r<nr;r++) { starts[r+1]+=starts[r]; cursor[r]=starts[r]; }
  for(unsigned int i=0;i<matrix->nelems;i++) indices[cursor[spolymat_get_row(matrix,i)]++]=i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int i=0;i<matrix->nelems;i++) _tocrt(spolymat_get_elem(matrix,i));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nizk2_parallel_allowed())
#endif
  for(unsigned int j=0;j<vec->nelems;j++) _tocrt(polyvec_get_elem(vec,j));
  polyvec_set_zero(out);
#ifdef _OPENMP
#pragma omp parallel if(nizk2_parallel_allowed())
#endif
  {
    poly_t product; poly_alloc(product,matrix->ring);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int r=0;r<nr;r++) {
      poly_ptr dst=polyvec_get_elem(out,r);
      for(unsigned int j=starts[r];j<starts[r+1];j++) {
        unsigned int i=indices[j];
        poly_mul(product,spolymat_get_elem(matrix,i),polyvec_get_elem(vec,spolymat_get_col(matrix,i)));
        poly_add(dst,dst,product,1);
      }
    }
    poly_free(product);
  }
  free(indices); free(cursor); free(starts);
  return 1;
}
