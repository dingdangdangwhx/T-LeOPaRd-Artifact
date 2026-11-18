#include "threading.h"
/* Hermitian symmetry halves public Gram dot products. No relation changes. */
static void layout_affine_gram(spolymat_t R2, spolyvec_t r1, polymat_ptr Es,
    polymat_ptr Em, polyvec_ptr v, polymat_ptr oEs, polymat_ptr oEm,
    polyvec_ptr ov, unsigned int m1, unsigned int l, unsigned int Z,
    unsigned int nr, polyring_srcptr ring) {
  const unsigned int nc=m1+l;
  unsigned int *active=malloc(nc*sizeof(*active)), na=0;
  poly_ptr *values=calloc(nc*nr,sizeof(*values)), *autos=calloc(nc*nr,sizeof(*autos));
  poly_ptr *vv=calloc(nr,sizeof(*vv));
  if(!active || !values || !autos || !vv) abort();
  for(unsigned int c=0;c<nc;c++) {
    polymat_ptr mat=c<m1 ? Es : Em, aut=c<m1 ? oEs : oEm;
    if(!mat) continue;
    unsigned int col=c<m1 ? c : c-m1, used=0;
    for(unsigned int r=0;r<nr;r++) {
      poly_ptr x=polymat_get_elem(mat,r,col);
      if(!x->crt && __lnp_public_poly_is_zero(x)) continue;
      values[c*nr+r]=x; autos[c*nr+r]=polymat_get_elem(aut,r,col); used=1;
    }
    if(used) active[na++]=c;
  }
  if(v) for(unsigned int r=0;r<nr;r++) {
    poly_ptr x=polyvec_get_elem(v,r);
    if(x->crt || !__lnp_public_poly_is_zero(x)) vv[r]=x;
  }
  /* All lazy CRT conversions finish before shared inputs are read in parallel. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nc*nr>=4096 && nizk2_parallel_allowed())
#endif
  for(unsigned int c=0;c<nc;c++) for(unsigned int r=0;r<nr;r++)
    if(values[c*nr+r]) { _tocrt(values[c*nr+r]); _tocrt(autos[c*nr+r]); }
  for(unsigned int r=0;r<nr;r++) if(vv[r]) _tocrt(vv[r]);
  unsigned int *ja=malloc(nc*nc*sizeof(*ja)), *jb=malloc(nc*nc*sizeof(*jb)), nj=0;
  if(!ja || !jb) abort();
  for(unsigned int a=0;a<na;a++) for(unsigned int b=a;b<na;b++)
    for(unsigned int r=0;r<nr;r++) if(values[active[a]*nr+r] && values[active[b]*nr+r]) {
      ja[nj]=a; jb[nj++]=b; break;
    }
  polyvec_t gram,mirror,left,right;
  polyvec_alloc(gram,ring,nj); polyvec_alloc(mirror,ring,nj);
  polyvec_alloc(left,ring,na); polyvec_alloc(right,ring,na);
  polyvec_set_zero(gram); polyvec_set_zero(left);
#ifdef _OPENMP
#pragma omp parallel if(nj>=128 && nizk2_parallel_allowed())
#endif
  {
    poly_t product; poly_alloc(product,ring);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int j=0;j<nj;j++) {
      unsigned int c=active[ja[j]], cc=active[jb[j]];
      poly_ptr dst=polyvec_get_elem(gram,j);
      for(unsigned int r=0;r<nr;r++) if(values[c*nr+r] && values[cc*nr+r]) {
        poly_mul(product,autos[c*nr+r],values[cc*nr+r]);
        poly_add(dst,dst,product,1);
      }
      if(c!=cc) poly_auto(polyvec_get_elem(mirror,j),dst);
    }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int a=0;a<na;a++) {
      unsigned int c=active[a]; poly_ptr dst=polyvec_get_elem(left,a);
      for(unsigned int r=0;r<nr;r++) if(values[c*nr+r] && vv[r]) {
        poly_mul(product,autos[c*nr+r],vv[r]); poly_add(dst,dst,product,1);
      }
      poly_auto(polyvec_get_elem(right,a),dst);
    }
    poly_free(product);
  }
  for(unsigned int j=0;j<nj;j++) {
    unsigned int c=active[ja[j]], cc=active[jb[j]];
    unsigned int pc=c<m1 ? c : c+Z, qc=cc<m1 ? cc : cc+Z, row,col;
    poly_ptr val=polyvec_get_elem(gram,j);
    if(__lnp_public_poly_is_zero(val)) continue;
    _diag(&row,&col,2*pc+1,2*qc);
    poly_set(spolymat_insert_elem(R2,row,col),val);
    if(c!=cc) {
      _diag(&row,&col,2*qc+1,2*pc);
      poly_set(spolymat_insert_elem(R2,row,col),polyvec_get_elem(mirror,j));
    }
  }
  if(v) for(unsigned int a=0;a<na;a++) {
    unsigned int c=active[a], pc=c<m1 ? c : c+Z;
    poly_ptr val=polyvec_get_elem(left,a);
    if(!__lnp_public_poly_is_zero(val)) poly_set(spolyvec_insert_elem(r1,2*pc+1),val);
    val=polyvec_get_elem(right,a);
    if(!__lnp_public_poly_is_zero(val)) poly_set(spolyvec_insert_elem(r1,2*pc),val);
  }
  polyvec_free(right); polyvec_free(left); polyvec_free(mirror); polyvec_free(gram);
  free(jb); free(ja); free(vv); free(autos); free(values); free(active);
  (void)ov;
}
