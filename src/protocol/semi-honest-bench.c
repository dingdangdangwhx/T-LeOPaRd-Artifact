/* Proof-free protocol benchmark: shared client fast path and RNS arithmetic. */
#define LEOPARD_PROOF_FREE 1
#include "nizk1.c"
#include "gaussian_lnp.h"
#include "wire.h"

#if LOGQ == 56
#define SH_ONLINE_LOG 13
#elif LOGQ == 64
#define SH_ONLINE_LOG 20
#else
#error The maintained proof-free benchmark profiles use q56 or q64.
#endif

enum { SH_COMPONENTS = 3, SH_SERVERS = 3 };
static const unsigned int holder_masks[SH_COMPONENTS] = {3,5,6};

static double now_sec(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec + 1e-9*(double)t.tv_nsec;
}

static void canonical(polx *v, size_t count) {
  leopard_coeff *coeffs = checked_malloc(count*N*sizeof(*coeffs),"canonical coefficients");
  export_polx_coeffs(coeffs,v,count);
  import_polx_coeffs(v,coeffs,count);
  free(coeffs);
}

static void binary_vector(polx *out, size_t count) {
  poly *small = alloc_polyvec(count);
  uint8_t *bits = checked_malloc((count*N+7)/8,"binary coins");
  randombytes(bits,(count*N+7)/8);
  for(size_t i=0;i<count;i++)
    for(size_t j=0;j<N;j++) {
      size_t index=i*N+j;
      small[i].vec->c[j]=(bits[index/8]>>(index%8))&1;
    }
  polxvec_frompolyvec(out,small,count);
  free(bits);
  free(small);
}

static void gaussian_vector(polx *out, size_t count, unsigned int level,
                            const uint8_t *fixed_seed) {
  uint8_t seed[16];
  int64_t *coeffs=checked_malloc(count*N*sizeof(*coeffs),"Gaussian coefficients");
  if(fixed_seed != NULL) memcpy(seed,fixed_seed,sizeof(seed));
  else randombytes(seed,sizeof(seed));
  gaussian_lnp_i64(coeffs,count*N,seed,0,level);
  polxvec_fromint64vec(out,count,1,coeffs);
  free(coeffs);
}

static void matrix_vector(polx *out, const polx *matrix, const polx *v,
                          size_t rows, size_t cols) {
  for(size_t i=0;i<rows;i++)
    polxvec_sprod(out+i,matrix+i*cols,v,cols);
}

static void full_sprod(polx *out, const polx *a, const polx *b, size_t len) {
  /* Bound each full-q integer convolution before CRT reconstruction. */
  const size_t chunk=LOGQ==64 ? 4 : 32;
  polx part[1];
  polxvec_setzero(out,1);
  for(size_t i=0;i<len;i+=chunk) {
    size_t count=len-i<chunk ? len-i : chunk;
    polxvec_sprod(part,a+i,b+i,count);
    canonical(part,1);
    polx_add(out,out,part);
  }
}

static void output_hash(uint8_t out[32], const polx *z, const uint8_t tag[16],
                        const uint8_t bits[80], const uint8_t pk_hash[32]) {
  uint8_t rounded[64];
  shake128incctx ctx;
  leopard_round_polx_p(rounded,z);
  shake128_inc_init(&ctx);
  absorb_bytes(&ctx,"T-LeOPaRd-proof-free-output",26);
  absorb_bytes(&ctx,tag,16);
  absorb_bytes(&ctx,bits,80);
  absorb_bytes(&ctx,rounded,sizeof(rounded));
  absorb_bytes(&ctx,pk_hash,32);
  shake128_inc_finalize(&ctx);
  shake128_inc_squeeze(out,32,&ctx);
}

static void offset_and_eval_noise(polx *offset, polx *eval_noise,
                                  const uint8_t pk_hash[32],
                                  const uint8_t secret[16],
                                  const uint8_t tag[16], const uint8_t bits[80]) {
  shake128incctx ctx;
  uint8_t seed[16], bytes[QBYTES];
  leopard_coeff coeffs[64];
  const unsigned __int128 q=((unsigned __int128)1<<LOGQ)-QOFF;
  shake128_inc_init(&ctx);
  absorb_bytes(&ctx,"T-LeOPaRd-ROpk",14);
  absorb_bytes(&ctx,pk_hash,32);
  absorb_bytes(&ctx,tag,16);
  absorb_bytes(&ctx,bits,80);
  shake128_inc_finalize(&ctx);
  for(size_t i=0;i<N;i++) {
    unsigned __int128 value;
    do {
      value=0;
      shake128_inc_squeeze(bytes,sizeof(bytes),&ctx);
      for(size_t j=0;j<sizeof(bytes);j++) value|=(unsigned __int128)bytes[j]<<(8*j);
    } while(value>=q);
    coeffs[i]=(leopard_coeff)(value>q/2 ? (__int128)value-(__int128)q : (__int128)value);
  }
  import_polx_coeffs(offset,coeffs,1);
  if(eval_noise==NULL) return;
  shake128_inc_init(&ctx);
  absorb_bytes(&ctx,"T-LeOPaRd-ROE",13);
  absorb_bytes(&ctx,secret,16);
  absorb_bytes(&ctx,tag,16);
  absorb_bytes(&ctx,bits,80);
  shake128_inc_finalize(&ctx);
  shake128_inc_squeeze(seed,sizeof(seed),&ctx);
  gaussian_vector(eval_noise,1,SH_ONLINE_LOG,seed);
}

static int run_once(size_t run, int report) {
  const size_t m=DEMO_M, ell=UMC_ELL, rows=m+ell, rank=ell;
  const unsigned int quorums[3]={5,3,6};
  leopard_pp pp;
  leopard_preq_data pre_public={0};
  leopard_pre_state_data pre={0};
  leopard_nizk1_fast_public_profile fast={0}, measured_fast={0};
  uint8_t secret[16], seed_hash[32], pk_hash[32], tag[16], bits[80];
  uint8_t input[10]={0x5a,0xc3,0x19,0x80,0x27,0x71,0x04,0xbe,0xd2,0x61};
  const uint8_t tag_bytes[]="POPRF-tag-00001";
  polx *keys=alloc_polxvec(3*m), *commitments=alloc_polxvec(3*(rank+2*m));
  polx *master=alloc_polxvec(m), *A=alloc_polxvec(rank*rank);
  polx *B=alloc_polxvec(2*m*rank), *r=alloc_polxvec(rank), *e=alloc_polxvec(rank+2*m);
  polx *vk=alloc_polxvec(3*3*rows), *es=alloc_polxvec(3*3*rows);
  polx *ar=alloc_polxvec(m*rows), *ar_t=alloc_polxvec(m*rows);
  polx *cx=alloc_polxvec(m), *bx=alloc_polxvec(m), *mask=alloc_polxvec(rows);
  polx *key_sum=alloc_polxvec(m), *vsum=alloc_polxvec(rows), *esum=alloc_polxvec(rows);
  polx z[1], offset[1], eval_noise[1], eps[1], u[1], usum[1], epssum[1], tmp[1], exact[1];
  double times[7]={0}, server_times[2]={0}, started;
  /* State for the three test epochs, per server: 0 absent, 1 ready, 2 spent. */
  unsigned char epoch_state[3][SH_SERVERS]={{0}};
  int rc=1;

  started=now_sec();
  randombytes((uint8_t *)&pp,sizeof(pp));
  times[0]=now_sec()-started;

  started=now_sec();
  randombytes(secret,sizeof(secret));
  shake128(seed_hash,sizeof(seed_hash),secret,sizeof(secret));
  leopard_public_almostuniform(A,rank*rank,pp.nizk1_seed,7000,leopard_public_threads());
  leopard_public_almostuniform(B,2*m*rank,pp.nizk1_seed,8000,leopard_public_threads());
  polxvec_setzero(master,m);
  for(size_t c=0;c<3;c++) {
    uint8_t seed[16];
    polx *commit=commitments+c*(rank+2*m), *key=keys+c*m;
    randombytes(seed,sizeof(seed));
    polxvec_almostuniform(key,m,seed,0);
    canonical(key,m);
    binary_vector(r,rank);
    binary_vector(e,rank+2*m);
    matrix_vector(commit,A,r,rank,rank);
    matrix_vector(commit+rank,B,r,2*m,rank);
    for(size_t i=0;i<rank+2*m;i++) polx_add(commit+i,commit+i,e+i);
    const unsigned __int128 q=((unsigned __int128)1<<LOGQ)-QOFF;
    const uint64_t scale=(uint64_t)sqrtl((long double)q);
    uint64_t a=scale;
    while((unsigned __int128)a*a>q) a--;
    while((unsigned __int128)(a+1)*(a+1)<=q) a++;
    for(size_t i=0;i<m;i++) {
      polx_add(commit+rank+i,commit+rank+i,key+i);
      polx_scale_add(commit+rank+m+i,key+i,a);
      polx_add(master+i,master+i,key+i);
    }
  }
  canonical(master,m);
  shake128incctx hash;
  shake128_inc_init(&hash);
  absorb_bytes(&hash,"T-LeOPaRd-PK",12);
  for(size_t i=0;i<3*(rank+2*m);i++) absorb_polx(&hash,commitments+i);
  absorb_bytes(&hash,seed_hash,sizeof(seed_hash));
  shake128_inc_finalize(&hash);
  shake128_inc_squeeze(pk_hash,sizeof(pk_hash),&hash);
  times[1]=now_sec()-started;

  for(size_t qi=0;qi<3;qi++) {
    uint8_t epoch[8]={0}, material[sizeof("T-LeOPaRd-epoch-v1")-1+32+8], seed[16];
    epoch[7]=(uint8_t)(1+qi);
    memcpy(material,"T-LeOPaRd-epoch-v1",sizeof(material)-40);
    memcpy(material+sizeof(material)-40,pk_hash,32);
    memcpy(material+sizeof(material)-8,epoch,8);
    started=now_sec();
    for(size_t server=0;server<3;server++) {
      if(epoch_state[qi][server] != 0) goto cleanup;
      shake128(seed,16,material,sizeof(material));
      leopard_public_almostuniform(ar,m*rows,seed,NONCE_AR_BASE,leopard_public_threads());
      canonical(ar,m*rows);
      for(size_t i=0;i<rows;i++)
        for(size_t k=0;k<m;k++) polxvec_copy(ar_t+i*m+k,ar+k*rows+i,1);
      for(size_t c=0;c<3;c++) if(holder_masks[c]&(1u<<server)) {
        size_t index=(server*3+c)*rows;
        gaussian_vector(es+index,rows,3,NULL);
        for(size_t k=0;k<rows;k++) full_sprod(vk+index+k,ar_t+k*m,keys+c*m,m);
        for(size_t i=0;i<rows;i++) polx_add(vk+index+i,vk+index+i,es+index+i);
        canonical(vk+index,rows);
      }
      epoch_state[qi][server]=1;
    }
    if(qi==0) times[2]=(now_sec()-started)/3;
    uint8_t finalized[32], evaluated[32], repeated[32];
    started=now_sec();
    shake128(seed,16,material,sizeof(material));
    leopard_client_epoch(seed,1);
    leopard_pre_state_data_clear(&pre);
    leopard_preq_data_clear(&pre_public);
    if(leopard_nizk1_preproc_client_data(&pre_public,&pre,&pp,1)) goto cleanup;
    leopard_tag_hash(tag,tag_bytes,sizeof(tag_bytes)-1);
    leopard_expand_input_bits(bits,input);
    if(leopard_fast_public_compute_preprocessed(&fast,&pre,0,tag,bits,NULL,cx)) goto cleanup;
    canonical(cx,m);
    import_polx_coeffs(mask,pre.rhat_coeffs,rows);
    double request_time=now_sec()-started;

    polxvec_setzero(usum,1); polxvec_setzero(epssum,1);
    polxvec_setzero(vsum,rows); polxvec_setzero(esum,rows);
    size_t active=0;
    for(size_t server=0;server<3;server++) if(quorums[qi]&(1u<<server)) {
      started=now_sec();
      if(epoch_state[qi][server] != 1) goto cleanup;
      polxvec_setzero(key_sum,m);
      for(size_t c=0;c<3;c++) {
        unsigned int intersection=quorums[qi]&holder_masks[c];
        unsigned int owner=31u-(unsigned int)__builtin_clz(intersection);
        if(owner!=server) continue;
        for(size_t k=0;k<m;k++) polx_add(key_sum+k,key_sum+k,keys+c*m+k);
      }
      canonical(key_sum,m);
      gaussian_vector(eps,1,SH_ONLINE_LOG,NULL);
      full_sprod(u,cx,key_sum,m);
      polx_add(u,u,eps);
      /* Coefficient export is the existing reply-encoding preparation. */
      leopard_coeff coeffs[N];
      export_polx_coeffs(coeffs,u,1);
      epoch_state[qi][server]=2;
      double elapsed=now_sec()-started;
      if(qi==0) server_times[active]=elapsed;
      active++;
      polx_add(usum,usum,u); polx_add(epssum,epssum,eps);
    }
    if(active!=2) goto cleanup;

    started=now_sec();
    for(size_t c=0;c<3;c++) {
      unsigned int intersection=quorums[qi]&holder_masks[c];
      unsigned int owner=31u-(unsigned int)__builtin_clz(intersection);
      size_t index=(owner*3+c)*rows;
      for(size_t k=0;k<rows;k++) polx_add(vsum+k,vsum+k,vk+index+k);
    }
    polxvec_sprod(tmp,mask,vsum,rows);
    polx_sub(z,usum,tmp);
    offset_and_eval_noise(offset,NULL,pk_hash,NULL,tag,bits);
    polx_add(z,z,offset);
    output_hash(finalized,z,tag,bits,pk_hash);
    double finalize_time=now_sec()-started;

    started=now_sec();
    if(leopard_compute_bxt_fast(bx,tag,bits)) goto cleanup;
    canonical(bx,m);
    full_sprod(exact,bx,master,m);
    offset_and_eval_noise(offset,eval_noise,pk_hash,secret,tag,bits);
    polx_add(tmp,exact,offset); polx_add(tmp,tmp,eval_noise);
    output_hash(evaluated,tmp,tag,bits,pk_hash);
    double eval_time=now_sec()-started;
    offset_and_eval_noise(offset,eval_noise,pk_hash,secret,tag,bits);
    polx_add(tmp,exact,offset); polx_add(tmp,tmp,eval_noise);
    output_hash(repeated,tmp,tag,bits,pk_hash);
    /* Check cancellation before rounding, separately from output agreement. */
    for(size_t c=0;c<3;c++) {
      unsigned int intersection=quorums[qi]&holder_masks[c];
      unsigned int owner=31u-(unsigned int)__builtin_clz(intersection);
      size_t index=(owner*3+c)*rows;
      for(size_t k=0;k<rows;k++) polx_add(esum+k,esum+k,es+index+k);
    }
    polxvec_sprod(tmp,mask,esum,rows);
    polx_sub(tmp,epssum,tmp); polx_add(tmp,tmp,exact); polx_add(tmp,tmp,offset);
    leopard_coeff actual[N], expected[N];
    export_polx_coeffs(actual,z,1); export_polx_coeffs(expected,tmp,1);
    if(memcmp(actual,expected,sizeof(actual))) {
      fprintf(stderr,"exact cancellation mismatch quorum=%u\n",quorums[qi]);
      goto cleanup;
    }
    if(memcmp(finalized,evaluated,32)||memcmp(evaluated,repeated,32)) {
      fprintf(stderr,"rounding/output mismatch quorum=%u\n",quorums[qi]);
      goto cleanup;
    }
    if(qi==0) {
      times[3]=request_time;
      measured_fast=fast;
      times[4]=(server_times[0]+server_times[1])/2;
      times[5]=finalize_time; times[6]=eval_time;
    }
  }
  if(report) {
    printf("{\"run\":%zu,\"logq\":%d,\"m\":%zu,\"ell\":%zu,\"p\":%d,\"sigma_es\":12.4,\"sigma_online\":%.1f,\"correctness\":\"PASS\",\"quorums\":3,\"proofs_executed\":0,\"seconds\":[",
           run,LOGQ,m,ell,LEOPARD_ROUND_P,1.55*(1u<<SH_ONLINE_LOG));
    for(size_t i=0;i<7;i++) printf("%s%.9f",i?",":"",times[i]);
    printf("],\"request_breakdown\":[%.9f,%.9f,%.9f,%.9f,%.9f],\"bytes\":[48,%zu,%zu,%zu,%zu,0,0]}\n",
           measured_fast.total_sec,measured_fast.selector_material_sec,measured_fast.slot_load_sec,measured_fast.chain_sec,measured_fast.cx_sec,
           3*(rank+2*m)*N*QBYTES+32,
           2*rows*N*QBYTES+8,m*N*QBYTES+24+pre.cr_rows*N*10/8,(size_t)N*QBYTES+28);
    fflush(stdout);
  }
  rc=0;
cleanup:
  leopard_preq_data_clear(&pre_public); leopard_pre_state_data_clear(&pre);
  free(esum); free(vsum); free(key_sum); free(mask); free(bx); free(cx);
  free(ar_t); free(ar); free(es); free(vk); free(e); free(r); free(B); free(A);
  free(master); free(commitments); free(keys);
  if(rc) fprintf(stderr,"proof-free correctness check failed (run %zu)\n",run);
  return rc;
}

int main(int argc, char **argv) {
  size_t runs=argc>1 ? strtoul(argv[1],NULL,10) : 3;
  if(runs==0 || runs>1000 || LEOPARD_ROUND_P!=4) return 1;
  for(size_t i=0;i<=runs;i++) if(run_once(i,i!=0)) return 1;
  return 0;
}
