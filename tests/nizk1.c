#include "../src/protocol/nizk1.c"

static void check(int ok, const char *name) {
  if(!ok) { fprintf(stderr,"FAIL: %s\n",name); exit(1); }
  printf("PASS: %s\n",name);
}

static int direct_probe(size_t count) {
  /* This is a small backend relation probe, not a T-OPRF benchmark. */
  check(count>0 && count<=64,"direct probe resource guard");
  witness wt={0}; prncplstmnt st={0}; composite cp={0};
  size_t *n=calloc(count,sizeof(*n));
  polx **sx=calloc(count,sizeof(*sx));
  uint8_t seed[16]={19};
  for(size_t i=0;i<count;i++) n[i]=1;
  init_witness_raw(&wt,count,n);
  uint64_t norm=0;
  for(size_t i=0;i<count;i++) {
    polyvec_quarternary(wt.s[i],1,seed,i);
    wt.normsq[i]=polyvec_sprodz(wt.s[i],wt.s[i],1);
    norm+=wt.normsq[i];
  }
  check(!init_prncplstmnt_raw(&st,count,n,norm+16*N,1,1),"direct statement init");
  init_sparsecnst_view(st.cnst,count,0,1,1,0);
  st.cnst->a->len=count;
  polxvec_quarternary(st.cnst->a->coeffs,count,seed,100);
  for(size_t i=0;i<count;i++) st.cnst->a->rows[i]=st.cnst->a->cols[i]=i;
  sparsecnst_eval(st.cnst->b,st.cnst,sx,&wt);
  shake128incctx hash;
  shake128_inc_init(&hash);
  absorb_bytes(&hash,"direct-weighted-route-probe",sizeof("direct-weighted-route-probe"));
  absorb_u64(&hash,count);
  absorb_polxvec(&hash,st.cnst->a->coeffs,count);
  absorb_polxvec(&hash,st.cnst->b,1);
  shake128_inc_finalize(&hash);
  shake128_inc_squeeze(st.h,sizeof(st.h),&hash);
  free(sx[0]); free(sx); free(n);
  check(!principle_verify(&st,&wt),"direct weighted quadratic relation");
  struct timespec a,b,c;
  clock_gettime(CLOCK_MONOTONIC,&a);
  int proved=composite_prove_principle(&cp,&st,&wt);
  clock_gettime(CLOCK_MONOTONIC,&b);
  check(!proved,"direct weighted full proof");
  int verified=composite_verify_principle(&cp,&st);
  clock_gettime(CLOCK_MONOTONIC,&c);
  check(!verified,"direct weighted full verification");
  size_t rr=0;
  for(size_t i=0;i<cp.pi[0]->r;i++) rr+=cp.pi[0]->nu[i];
  printf("DIRECT count=%zu prove=%.6f verify=%.6f first_r=%zu gram=%zu\n",
         count,elapsed_sec(&a,&b),elapsed_sec(&b,&c),rr,rr*(rr+1)/2);
  polx one; polx_monomial(&one,1,0);
  polx_add(st.cnst->b,st.cnst->b,&one);
  check(composite_verify_principle(&cp,&st)!=0,"direct changed equation rejected");
  free_composite(&cp); free_prncplstmnt(&st); free_witness(&wt); free_comkey();
  return 0;
}

int main(int argc,char **argv) {
  setvbuf(stdout,NULL,_IONBF,0);
  configure_openmp_threads();
  if(argc>1) return direct_probe(strtoul(argv[1],NULL,10));
  outer_instance inst={0}; witness wt={0}; prncplstmnt st={0}; composite cp={0};
  check(!generate_inner_instance(&inst,2),"ARP inner proof and verifier");
  selector_lnp22_ctx *ctx=&inst.ctx;
  const size_t row=RBIT_CONST_ROWS;
  polx saved=ctx->rbit_const_h[row], bad;
  polx_monomial(&bad,RBIT_ARP_LIMIT+1,0);
  ctx->rbit_const_h[row]=bad;
  check(validate_rbit_const_transcript(ctx)!=0,"positive ARP overflow rejected");
  polx_monomial(&ctx->rbit_const_h[row],-RBIT_ARP_LIMIT-1,0);
  check(validate_rbit_const_transcript(ctx)!=0,"negative ARP overflow rejected");
  polz z={0}; polz_setcoeff_fromint128(&z,(__int128)1<<70,0);
  polz_topolx(&ctx->rbit_const_h[row],&z);
  check(validate_rbit_const_transcript(ctx)!=0,"full-width ARP overflow rejected");
  ctx->rbit_const_h[row]=saved;
  check(!validate_rbit_const_transcript(ctx),"valid ARP transcript restored");
  polx lhs,rhs,delta;
  polx_monomial(&delta,1,0);
  folded_linear_eval(&lhs,ctx,ctx->s1,inst.s2);
  compute_public_fold_rhs(&rhs,ctx);
  polx_sub(&rhs,&rhs,&lhs);
  check(polx_iszero(&rhs),"ARP projection equations included in linear fold");
  polx_add(&ctx->rbit_const_h[row],&saved,&delta);
  check(!validate_rbit_const_transcript(ctx),"small ARP tamper still passes range check alone");
  compute_linear_constants(ctx);
  compute_public_fold_rhs(&rhs,ctx);
  polx_sub(&rhs,&rhs,&lhs);
  check(!polx_iszero(&rhs),"projection relation rejects small ARP tamper");
  ctx->rbit_const_h[row]=saved;
  compute_linear_constants(ctx);
  saved=ctx->rbit_const_h[0];
  polx_add(&ctx->rbit_const_h[0],&saved,&delta);
  check(validate_rbit_const_transcript(ctx)!=0,"constantness check remains enforced");
  ctx->rbit_const_h[0]=saved;
  uint8_t oldseed[16]; memcpy(oldseed,ctx->rbit_const_seed,16);
  saved=ctx->rbit_const_tg[row]; polx one; polx_monomial(&one,1,0);
  polx_add(&ctx->rbit_const_tg[row],&saved,&one);
  derive_rbit_const_seed(ctx);
  check(memcmp(oldseed,ctx->rbit_const_seed,16)!=0,"ARP matrix bound to mask commitment");
  ctx->rbit_const_tg[row]=saved; derive_rbit_const_seed(ctx);
  check(!memcmp(oldseed,ctx->rbit_const_seed,16),"ARP challenge restored");
  polx pair[2], product, value, phi;
  polx_monomial(&pair[0],3,0);
  const __int128 fixture=(__int128)3161622474400LL*1000000000+472163928;
  polz_setcoeff_fromint128(&z,fixture,0); polz_topolx(&pair[1],&z);
  polxvec_sprod(&value,pair,pair,2);
  polx_sub(&value,&value,&pair[0]); polx_sub(&value,&value,&pair[1]);
  check(polx_iszero(&value),"modular non-bit fixture passes bare unweighted sum");
  saved=ctx->rbit_const_h[row];
  polxvec_setzero(&value,1);
  for(size_t i=0;i<2;i++) {
    sample_rbit_const_phi(&phi,ctx,row,i);
    polx_mul(&product,&phi,&pair[i]); polx_add(&value,&value,&product);
  }
  ctx->rbit_const_h[row]=value;
  check(validate_rbit_const_transcript(ctx)!=0,"fixed-mask projection rejects modular fixture");
  ctx->rbit_const_h[row]=saved;
  uint64_t budget=0;
  leopard_outer_statement_profile profile={0};
  check(!leopard_prepare_outer_witness(&wt,&inst,&budget),"ARP outer witness");
  check(!leopard_prepare_outer_statement(&st,ctx,inst.tA,inst.w,inst.c,
          inst.vlin,inst.tq,inst.vq,budget,0,&profile),"ARP outer statement");
  check(!principle_verify(&st,&wt),"ARP outer algebraic checks");
  if(OUTER_PACK_LINEAR_TAIL) {
    check(wt.r==OUT_R-3 && st.r==OUT_R-3,"linear-tail block count reduced");
    uint64_t sum=0;
    for(size_t i=0;i<wt.r;i++) {
      check(wt.n[i]==st.n[i],"packed witness/statement dimensions agree");
      sum+=wt.normsq[i];
    }
    check(sum+16*N<=st.betasq,"packed witness retains original joint norm budget");
    const size_t tail=ctx->rhat_len+ctx->umc_eua_len+ctx->umc_eub_len+
                      ctx->umc_r_len+ctx->umc_va_len+ctx->umc_vb_len;
    const size_t offsets[]={0,tail,2*tail,2*tail+ctx->m2};
    for(size_t i=0;i<4;i++) {
      int16_t *coeff=&wt.s[OUT_T0][offsets[i]].vec->c[0];
      const int16_t saved_coeff=*coeff;
      *coeff=(int16_t)(saved_coeff==32767 ? saved_coeff-1 : saved_coeff+1);
      check(principle_verify(&st,&wt)!=0,"each packed segment remains constrained");
      *coeff=saved_coeff;
    }
    check(!principle_verify(&st,&wt),"packed witness restored");
  }
  check(!composite_prove_principle(&cp,&st,&wt),"ARP complete outer proof");
  check(!composite_verify_principle(&cp,&st),"ARP complete outer verification");
  free_composite(&cp); free_prncplstmnt(&st); free_witness(&wt);
  free_outer_instance(&inst); free_comkey();
  return 0;
}
