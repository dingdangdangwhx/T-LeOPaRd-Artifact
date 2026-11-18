#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../src/protocol/nizk2.c"

#define CHECK(c) do { if(!(c)) { fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c); return 1; } } while(0)

static leopard_coeff centered128(__int128 x) {
  const __int128 q=(unsigned __int128)LEOPARD_MOD_Q;
  x%=q;
  if(x>q/2) x-=q;
  if(x< -q/2) x+=q;
  return (leopard_coeff)x;
}

static int extraction_test(void) {
  polyring_srcptr ring=params->tbox_params->tbox->ring;
  leopard_pp pp={{1},{2}};
  uint8_t seed[32]={17};
  polymat_t A,B,S,E,G;
  polyvec_t r,ea,eb,k,ca,cb,tmp,u,row,srow,erow;
  polymat_alloc(A,ring,COM_RANK,COM_RANK);
  polymat_alloc(B,ring,2*N,COM_RANK);
  polymat_alloc(S,ring,2*N,COM_RANK);
  polymat_alloc(E,ring,2*N,COM_RANK);
  polymat_alloc(G,ring,2*N,N);
  leopard_server_matrix_a1(A,&pp);
  for(unsigned int i=0;i<2*N;i++) {
    polymat_get_row(srow,S,i); sample_polyvec_bnd(srow,-1,1,seed,100+i);
    polymat_get_row(erow,E,i); sample_polyvec_bnd(erow,-1,1,seed,200+i);
    polymat_get_row(row,B,i); polyvec_mul2(row,srow,A); polyvec_add(row,row,erow,0);
  }
  leopard_commitment_encoding(G);
  polyvec_alloc(r,ring,COM_RANK); polyvec_alloc(ea,ring,COM_RANK);
  polyvec_alloc(eb,ring,2*N); polyvec_alloc(k,ring,N);
  polyvec_alloc(ca,ring,COM_RANK); polyvec_alloc(cb,ring,2*N);
  polyvec_alloc(tmp,ring,2*N); polyvec_alloc(u,ring,2*N);
  leopard_coeff *kc=alloc_coeffs(N), *uc=alloc_coeffs(2*N);
  const leopard_coeff half=(leopard_coeff)((leopard_ucoeff)LEOPARD_MOD_Q/2);
  const leopard_coeff a=(leopard_coeff)leopard_double_scale();
  for(unsigned int trial=0;trial<8;trial++) {
    unsigned int mult=1+trial%2;
    sample_polyvec_bnd(r,0,mult,seed,1000+trial*4);
    sample_polyvec_bnd(ea,0,mult,seed,1001+trial*4);
    sample_polyvec_bnd(eb,0,mult,seed,1002+trial*4);
    leopard_uniform_key(k,seed,1003+trial*4);
    root_to_raw_polyvec(kc,k);
    kc[0]=half; kc[1]=-half; kc[2]=0; kc[3]=-1;
    raw_to_root_polyvec(k,kc);
    polyvec_mul(ca,A,r); polyvec_add(ca,ca,ea,0);
    polyvec_mul(cb,B,r); polyvec_add(cb,cb,eb,0);
    polyvec_mul(tmp,G,k); polyvec_add(cb,cb,tmp,0);
    polyvec_mul(tmp,S,ca); polyvec_sub(u,cb,tmp,0);
    root_to_raw_polyvec(uc,u);
    for(unsigned int i=0;i<N*D;i++) {
      leopard_coeff w=centered128((__int128)uc[N*D+i]-(__int128)a*uc[i]);
      __int128 v=-(__int128)w;
      leopard_coeff e0=(leopard_coeff)(v>=0 ? (v+a/2)/a : -((-v+a/2)/a));
      CHECK(centered128((__int128)uc[i]-e0)==kc[i]);
    }
  }
  free(uc); free(kc);
  polyvec_free(u); polyvec_free(tmp); polyvec_free(cb); polyvec_free(ca);
  polyvec_free(k); polyvec_free(eb); polyvec_free(ea); polyvec_free(r);
  polymat_free(G); polymat_free(E); polymat_free(S); polymat_free(B); polymat_free(A);
  puts("double-encryption extraction: random and boundary keys, single/aggregate PASS");
  return 0;
}

int main(void) {
  lazer_init();
  CHECK(NIZK2_DOUBLE_COMMITMENT && params->nm_indices==N);
  CHECK(extraction_test()==0);
  leopard_pp pp={{1},{2}};
  leopard_pk_data pk={0}; leopard_sk_data sk={0};
  leopard_request_data req={0};
  polyring_srcptr ring=params->tbox_params->tbox->ring;
  CHECK(leopard_nizk2_keygen_data(&pk,&sk,&pp)==0);
  CHECK(sk.ei_rows==2*N && pk.ck_rows==COM_RANK+2*N);
  {
    const leopard_coeff half=(leopard_coeff)((leopard_ucoeff)LEOPARD_MOD_Q/2);
    const leopard_coeff target[]={half,-half,0,-1};
    const leopard_coeff scale=(leopard_coeff)leopard_double_scale();
    for(unsigned int h=0;h<4;h++) {
      leopard_coeff delta=centered128((__int128)target[h]-sk.k_coeffs[h]);
      sk.k_coeffs[h]=target[h];
      pk.ck_coeffs[COM_RANK*D+h]=centered128((__int128)pk.ck_coeffs[COM_RANK*D+h]+delta);
      pk.ck_coeffs[(COM_RANK+N)*D+h]=centered128((__int128)pk.ck_coeffs[(COM_RANK+N)*D+h]+(__int128)scale*delta);
    }
  }
  {
    char path[]="/tmp/nizk2-double-pk-XXXXXX";
    int fd=mkstemp(path); CHECK(fd>=0); close(fd);
    leopard_pk_data copy={0};
    CHECK(leopard_write_pk_file(path,&pk)==0);
    CHECK(leopard_read_pk_file(path,&copy)==0);
    CHECK(copy.ck_rows==pk.ck_rows);
    CHECK(!memcmp(copy.ck_coeffs,pk.ck_coeffs,COM_ROWS*D*sizeof(leopard_coeff)));
    leopard_pk_data_clear(&copy); unlink(path);
  }
  req.cx_rows=H; req.cx_cols=N; req.cx_coeffs=alloc_coeffs(H*N);
  for(unsigned int i=0;i<H*N*D;i++) req.cx_coeffs[i]=(leopard_coeff)(i%7)-3;
  req.ar_seed[0]=9;
  {
    leopard_coeff *v0=alloc_coeffs(M+N), *v1=alloc_coeffs(M+N);
    leopard_coeff *e0=alloc_coeffs(M+N), *e1=alloc_coeffs(M+N);
    CHECK(leopard_compute_vk_from_ar_seed(v0,e0,&pp,&sk,req.ar_seed,0)==0);
    CHECK(leopard_compute_vk_from_ar_seed(v1,e1,&pp,&sk,req.ar_seed,0)==0);
    CHECK(memcmp(e0,e1,(M+N)*D*sizeof(leopard_coeff))!=0);
    free(e1); free(e0); free(v1); free(v0);
    puts("replicated holders use independent private error coins PASS");
  }
  for(unsigned int mult=1;mult<=2;mult++) {
    if(mult==2) {
      for(unsigned int i=0;i<COM_ROWS*D;i++) pk.ck_coeffs[i]=centered128((__int128)2*pk.ck_coeffs[i]);
      for(unsigned int i=0;i<N*D;i++) sk.k_coeffs[i]=centered128((__int128)2*sk.k_coeffs[i]);
      for(unsigned int i=0;i<COM_RANK*D;i++) { sk.r_coeffs[i]*=2; sk.e_coeffs[i]*=2; }
      for(unsigned int i=0;i<2*N*D;i++) sk.ei_coeffs[i]*=2;
    }
    polymat_t A; polyvec_t s,t,vk,ux,res;
    polymat_alloc(A,ring,REL_ROWS,WIT_ROWS);
    polyvec_alloc(s,ring,WIT_ROWS); polyvec_alloc(t,ring,REL_ROWS);
    polyvec_alloc(vk,ring,M+N); polyvec_alloc(ux,ring,H); polyvec_alloc(res,ring,REL_ROWS);
    CHECK(leopard_build_reply_instance(A,s,t,vk,ux,&pp,&pk,&sk,&req,NULL,NULL,NULL)==0);
    polyvec_mul(res,A,s); polyvec_add(res,res,t,0);
    poly_t zero; poly_alloc(zero,ring); poly_set_zero(zero);
    for(unsigned int i=0;i<REL_ROWS;i++) CHECK(poly_eq(polyvec_get_elem(res,i),zero));
    uint8_t *proof=calloc(lin_params_get_prooflen(params),1); CHECK(proof);
    double start=get_sec();
    size_t len=prover(proof,s,A,t,pp.nizk2_pp);
    printf("mult=%u prove=%.6f bytes=%zu\n",mult,get_sec()-start,len); fflush(stdout);
    CHECK(verifier(proof,A,t,pp.nizk2_pp)==1);
    if(mult==1 && getenv("NIZK2_L2_NEGATIVE_TESTS")) {
      const unsigned int badrows[]={0,COM_RANK,2*COM_RANK,2*COM_RANK+N};
      polyvec_t saved_t; polyvec_alloc(saved_t,ring,REL_ROWS); polyvec_set(saved_t,t);
      for(unsigned int j=0;j<4;j++) {
        poly_ptr x=polyvec_get_elem(s,badrows[j]);
        poly_fromcrt(x); poly_mod(x,x); poly_redc(x,x);
        INT_T(saved,ring->q->nlimbs); int_set(saved,poly_get_coeff(x,0));
        int_set_i64(poly_get_coeff(x,0),256);
        polyvec_mul(t,A,s); polyvec_neg_self(t); polyvec_mod(t,t); polyvec_redc(t,t);
        size_t badlen=prover(proof,s,A,t,pp.nizk2_pp);
        CHECK(badlen>0 && verifier(proof,A,t,pp.nizk2_pp)==0);
        poly_fromcrt(x); int_set(poly_get_coeff(x,0),saved);
        polyvec_set(t,saved_t);
      }
      polyvec_free(saved_t);
      len=prover(proof,s,A,t,pp.nizk2_pp);
      CHECK(verifier(proof,A,t,pp.nizk2_pp)==1);
      puts("L2 opening: excessive r/e_A/e_0/e_1 norm with valid algebra REJECTED");
    }
    const unsigned int rows[]={0,COM_RANK,COM_RANK+N,COM_ROWS,REL_ROWS-1};
    for(unsigned int j=0;j<sizeof(rows)/sizeof(rows[0]);j++) {
      poly_ptr v=polyvec_get_elem(t,rows[j]);
      poly_add(v,v,polymat_get_elem(A,0,COM_RANK),0);
      CHECK(verifier(proof,A,t,pp.nizk2_pp)==0);
      poly_sub(v,v,polymat_get_elem(A,0,COM_RANK),0);
    }
    if(mult==1) {
      /* A new proof for an impossible exact row must also be rejected. */
      polyvec_t row; polymat_get_row(row,A,REL_ROWS-1); polyvec_set_zero(row);
      poly_set_one(polyvec_get_elem(t,REL_ROWS-1));
      prover(proof,s,A,t,pp.nizk2_pp);
      CHECK(verifier(proof,A,t,pp.nizk2_pp)==0);
    }
    free(proof); poly_free(zero); polyvec_free(res); polyvec_free(ux); polyvec_free(vk);
    polyvec_free(t); polyvec_free(s); polymat_free(A);
    puts("full-modulus key proof and all exact relation negative tests PASS");
  }
  free(req.cx_coeffs); leopard_pk_data_clear(&pk); leopard_sk_data_clear(&sk);
  return 0;
}
