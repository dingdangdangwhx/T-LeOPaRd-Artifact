/*
 * Specialized outer frontend for the direct selector-chain relation.
 *
 * The inner response no longer contains Z_t = x_t B_{t+1}.  Instead, the
 * outer witness contains a decomposed auxiliary
 *
 *   P_t = <directD_t, zB_{t+1}> - beta_t zX_t
 *
 * for t < L-1, and P_{L-1} = -beta_{L-1} zX_{L-1}.  The beta values are
 * uniform ring challenges derived after t_A.  Consequently -<zX,P>
 * contains both the selector transition products and the randomized
 * coefficient-bit equation, without adding a separate witness block.
 *
 * P is laid out like the compact direct X vector, allowing Labrador's
 * whole-vector dot product to express -<zX,P>.  Every P
 * coordinate has its own degree-one constraint.  Chihuahua combines these
 * constraints only after committing to the outer witness, so the folding
 * challenges cannot be anticipated by the prover.
 */

#if LNP22_USE_NATIVE_FOLD_BINARY
static int build_direct_raut_postcommit_constraints(
    sparsecnst **out, size_t *out_k, const prncplstmnt *st,
    const uint8_t postcommit_h[16]) {
  const size_t raut_len = st->n[OUT_RAUT0];
  sparsecnst *cnst;
  polx *alpha, *salpha;
  size_t round;

  *out = NULL;
  *out_k = 0;
  if(raut_len == 0 || st->n[OUT_RAUT1] != raut_len ||
     st->n[OUT_T0] < raut_len || st->n[OUT_T1] < raut_len)
    return 1;

  cnst = _malloc(OUTER_RAUT_ALPHA_ROUNDS*sizeof(*cnst));
  memset(cnst,0,OUTER_RAUT_ALPHA_ROUNDS*sizeof(*cnst));
  alpha = alloc_polxvec(2*raut_len);
  salpha = &alpha[raut_len];
  for(round=0;round<OUTER_RAUT_ALPHA_ROUNDS;round++) {
    polx *buf;

    polxvec_almostuniform(alpha,raut_len,postcommit_h,
                          0x52415554ULL + round);
    polxvec_sigmam1(salpha,alpha,raut_len);
    buf = init_sparsecnst_half(&cnst[round],st->r,4,4*raut_len,0,0,1);
    cnst[round].idx[0] = OUT_T0;
    cnst[round].off[0] = 0;
    cnst[round].len[0] = raut_len;
    cnst[round].mult[0] = 1;
    cnst[round].phi[0] = buf;
    buf += raut_len;
    cnst[round].idx[1] = OUT_T1;
    cnst[round].off[1] = 0;
    cnst[round].len[1] = raut_len;
    cnst[round].mult[1] = 1;
    cnst[round].phi[1] = buf;
    buf += raut_len;
    cnst[round].idx[2] = OUT_RAUT0;
    cnst[round].off[2] = 0;
    cnst[round].len[2] = raut_len;
    cnst[round].mult[2] = 1;
    cnst[round].phi[2] = buf;
    buf += raut_len;
    cnst[round].idx[3] = OUT_RAUT1;
    cnst[round].off[3] = 0;
    cnst[round].len[3] = raut_len;
    cnst[round].mult[3] = 1;
    cnst[round].phi[3] = buf;

    polxvec_copy(cnst[round].phi[0],salpha,raut_len);
    polxvec_scale(cnst[round].phi[1],salpha,raut_len,
                  OUTER_XY_DECOMP_BASE);
    polxvec_neg(cnst[round].phi[2],alpha,raut_len);
    polxvec_scale(cnst[round].phi[3],alpha,raut_len,
                  -OUTER_XY_DECOMP_BASE);
  }
  free(alpha);
  *out = cnst;
  *out_k = OUTER_RAUT_ALPHA_ROUNDS;
  return 0;
}
#endif

static void compute_direct_outer_p(polx *p, const selector_lnp22_ctx *ctx,
                                   const polx *z1) {
  const polx *zB = &z1[ctx->off[LNP22_B_BLOCK]];
  const polx *zX = &z1[ctx->off[LNP22_X_BLOCK]];
  polx weighted[1];
  size_t t;

  for(t=0;t<ctx->Lchain-1;t++) {
    polxvec_sprod(&p[t],&ctx->directD[t*ctx->gdigits],
                  &zB[(t+1)*ctx->gdigits],ctx->gdigits);
    polx_mul(weighted,&ctx->xbit_weight[t],&zX[t]);
    polx_sub(&p[t],&p[t],weighted);
    polx_refresh(&p[t]);
  }
  for(;t<ctx->x_len;t++) {
    polx_mul(&p[t],&ctx->xbit_weight[t],&zX[t]);
    polx_neg(&p[t],&p[t]);
    polx_refresh(&p[t]);
  }
}

static size_t direct_x_permutation(const selector_lnp22_ctx *ctx, size_t pos) {
  (void)ctx;
  return pos;
}

static void permute_direct_x_polx(polx *x, const selector_lnp22_ctx *ctx) {
  polx *tmp = alloc_polxvec(ctx->x_len);
  size_t i;

  for(i=0;i<ctx->x_len;i++)
    polxvec_copy(&tmp[i],&x[direct_x_permutation(ctx,i)],1);
  polxvec_copy(x,tmp,ctx->x_len);
  free(tmp);
}

static int decompose_direct_x_two_digits(poly *dst0, poly *dst1,
                                         const int64_t *src,
                                         const selector_lnp22_ctx *ctx,
                                         int64_t base) {
  size_t pos, j;
  int invalid = 0;

  for(pos=0;pos<ctx->x_len;pos++) {
    const size_t srcpos = direct_x_permutation(ctx,pos);
    for(j=0;j<N;j++) {
      const int64_t x = src[srcpos*N+j];
      const int64_t q = balanced_round_div_i64(x,base);
      const int64_t rem = x - base*q;

      invalid |= outer_abs_i64(rem) > 32767 || outer_abs_i64(q) > 32767;
      dst0[pos].vec[0].c[j] = (int16_t)rem;
      dst1[pos].vec[0].c[j] = (int16_t)q;
    }
  }
  if(invalid)
    fprintf(stderr,"ERROR: direct X decomposition digit out of int16 range\n");
  return invalid;
}

static void set_direct_product_coeff(polx *out, int64_t left,
                                     int64_t right, int sign) {
  set_const_coeff(out,right);
  polx_scale(out,out,(sign < 0) ? -left : left);
}

static void build_direct_outer_statement(
    prncplstmnt *st, selector_lnp22_ctx *ctx,
    const polx *tA, const polx *w, const polx *vlin,
    const polx *tq, const polx *vq, const polx *c,
    const polx *a3, const polx *ufold_public,
    uint64_t betasq, int transfer_ctx_storage) {
  const size_t r = OUT_R;
  const size_t p_constraints = ctx->x_len;
  const size_t opening_idx = p_constraints;
  const size_t linear_idx = opening_idx + 1;
  const size_t quadratic_idx = linear_idx + 1;
  size_t n[OUT_R];
  size_t tail_len, binary_len, i, pos;
  shake128incctx shakectx;
  polx *openfold = NULL, *open1 = NULL, *open2 = NULL, *openrhs = NULL;
  polx *linear1 = NULL, *linear2 = NULL;
  polx *quadB = NULL, *quadX = NULL, *quadbinary = NULL, *quadz2 = NULL;
  polx *pscale = NULL;
  polx *bitweight = ctx->xbit_weight;
  polx linrhs[1], quadrhs[1];
#if LNP22_USE_NATIVE_FOLD_BINARY
  polx tmpq[1], c2[1], ones[1];
#endif
  sparsecnst *cnst;

  tail_len = ctx->rhat_len + ctx->umc_eua_len + ctx->umc_eub_len +
             ctx->umc_r_len + ctx->umc_va_len + ctx->umc_vb_len;
  binary_len = LNP22_USE_NATIVE_FOLD_BINARY ? ctx->rhat_len :
                                              ctx->rhat_bits_len;

  n[OUT_B0] = ctx->b_len;
  n[OUT_B1] = ctx->b_len;
  n[OUT_X0] = ctx->x_len;
  n[OUT_X1] = ctx->x_len;
  n[OUT_P0] = p_constraints;
  n[OUT_P1] = p_constraints;
  n[OUT_P2] = p_constraints;
  n[OUT_P3] = p_constraints;
  n[OUT_P4] = p_constraints;
  n[OUT_T0] = tail_len;
  n[OUT_T1] = tail_len;
#if LNP22_USE_NATIVE_FOLD_BINARY
  n[OUT_RAUT0] = binary_len;
  n[OUT_RAUT1] = binary_len;
#else
  n[OUT_RBIT0] = binary_len;
  n[OUT_RBIT1] = binary_len;
#endif
  n[OUT_Z20] = ctx->m2;
  n[OUT_Z21] = ctx->m2;

  if(init_prncplstmnt_raw(st,r,n,betasq,p_constraints + 3,1)) {
    fprintf(stderr,"ERROR: failed to initialize direct outer statement\n");
    exit(1);
  }
#if LNP22_USE_NATIVE_FOLD_BINARY
  st->build_postcommit_cnst = build_direct_raut_postcommit_constraints;
#endif

  shake128_inc_init(&shakectx);
  absorb_outer_compact_statement(&shakectx,ctx,tA,w,vlin,tq,vq,c);
  absorb_u64(&shakectx,LNP22_USE_NATIVE_FOLD_BINARY ?
             OUTER_RAUT_ALPHA_ROUNDS : 0);

  openfold = alloc_polxvec(ctx->kmsis);
  open1 = alloc_polxvec(ctx->m1);
  open2 = alloc_polxvec(ctx->m2);
  openrhs = alloc_polxvec(1);
  quadB = alloc_polxvec(ctx->b_len);
  quadX = alloc_polxvec(ctx->x_len);
  quadbinary = alloc_polxvec(binary_len);
  quadz2 = alloc_polxvec(ctx->m2);
  pscale = alloc_polxvec(OUTER_DIRECT_P_DIGITS);
  if(!openfold || !open1 || !open2 || !openrhs || !quadB || !quadX ||
     !quadbinary || !quadz2 || !pscale) {
    fprintf(stderr,"ERROR: direct outer statement allocation failure\n");
    exit(1);
  }

  derive_outer_openfold(openfold,ctx,tA,w,vlin,tq,vq,c);
  compute_openfold_coeffs(open1,open2,openrhs,ctx,tA,w,c,openfold);
  polx_mul(linrhs,c,ufold_public);
  polx_add(linrhs,linrhs,vlin);

  for(i=0;i<ctx->b_len;i++)
    polx_mul(&quadB[i],c,&ctx->direct_lin_B[i]);
  for(i=0;i<ctx->x_len;i++) {
    polx_mul(&quadX[i],c,&bitweight[i]);
    polx_neg(&quadX[i],&quadX[i]);
  }
#if LNP22_USE_NATIVE_FOLD_BINARY
  for(i=0;i<ctx->x_len;i++) {
    polx_mul(tmpq,c,&ctx->native_binary_lin1[i]);
    polx_add(&quadX[i],&quadX[i],tmpq);
  }
  polx_allones(ones);
  polx_mul(tmpq,c,ctx->native_binary_scale);
  polx_mul(tmpq,tmpq,ones);
  polx_neg(tmpq,tmpq);
  for(i=0;i<binary_len;i++)
    polxvec_copy(&quadbinary[i],tmpq,1);
#else
  for(i=0;i<binary_len;i++)
    polx_neg(&quadbinary[i],c);
#endif
  polxvec_copy(quadz2,a3,ctx->m2);
#if LNP22_USE_NATIVE_FOLD_BINARY
  for(i=0;i<ctx->m2;i++) {
    polx_mul(tmpq,c,&ctx->native_binary_lin2[i]);
    polx_add(&quadz2[i],&quadz2[i],tmpq);
  }
#endif
  polx_mul(quadrhs,c,tq);
  polx_add(quadrhs,quadrhs,vq);
#if LNP22_USE_NATIVE_FOLD_BINARY
  polx_mul(c2,c,c);
  polx_mul(tmpq,c2,ctx->native_binary_uquad);
  polx_add(quadrhs,quadrhs,tmpq);
#endif
#if LNP22_USE_NATIVE_FOLD_BINARY
  polxvec_refresh(quadB,ctx->b_len);
  polxvec_refresh(quadX,ctx->x_len);
  polxvec_refresh(quadbinary,binary_len);
  polxvec_refresh(quadz2,ctx->m2);
  polx_refresh(quadrhs);
#endif

  for(i=0;i<OUTER_DIRECT_P_DIGITS;i++)
    set_const_coeff(&pscale[i],outer_pow_i64(OUTER_DIRECT_P_BASE,(int)i));

  linear1 = ctx->lincoeff;
  linear2 = ctx->lincoeff2;
  permute_direct_x_polx(&open1[ctx->off[LNP22_X_BLOCK]],ctx);
  permute_direct_x_polx(&linear1[ctx->off[LNP22_X_BLOCK]],ctx);
  permute_direct_x_polx(quadX,ctx);

  /*
   * P well-formedness proves
   *   P_t - <directD_t,zB_{t+1}> + beta_t zX_t = 0
   * for t < L-1, and P_{L-1} + beta_{L-1} zX_{L-1} = 0.
   */
  for(pos=0;pos<p_constraints;pos++) {
    const int has_transition = pos < ctx->Lchain-1;
    const size_t nz = OUTER_DIRECT_P_DIGITS + 2 +
                      (has_transition ? 2 : 0);
    size_t k = 0, pd;

    cnst = &st->cnst[pos];
    init_sparsecnst_view(cnst,r,nz,1,0,0);
    polxvec_setzero(cnst->b,1);
    for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++)
      append_cnst_offset_view(cnst,&k,OUT_P0 + pd,pos,1,&pscale[pd],1);
    if(has_transition) {
      const size_t boff = (pos+1)*ctx->gdigits;
      append_cnst_offset_view(cnst,&k,OUT_B0,boff,ctx->gdigits,
                              &ctx->directD[pos*ctx->gdigits],-1);
      append_cnst_offset_view(cnst,&k,OUT_B1,boff,ctx->gdigits,
                              &ctx->directD[pos*ctx->gdigits],
                              -OUTER_BZ_DECOMP_BASE);
    }
    append_cnst_offset_view(cnst,&k,OUT_X0,pos,1,&bitweight[pos],1);
    append_cnst_offset_view(cnst,&k,OUT_X1,pos,1,&bitweight[pos],
                            OUTER_XY_DECOMP_BASE);
  }

  cnst = &st->cnst[opening_idx];
  init_sparsecnst_view(cnst,r,
                       LNP22_USE_NATIVE_FOLD_BINARY ? 8 : 10,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,ctx->b_len,
                   &open1[ctx->off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,ctx->b_len,
                   &open1[ctx->off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,
                   &open1[ctx->off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,
                   &open1[ctx->off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,
                   &open1[ctx->off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,
                   &open1[ctx->off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
#if !LNP22_USE_NATIVE_FOLD_BINARY
  append_cnst_view(cnst,&i,OUT_RBIT0,binary_len,
                   &open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,binary_len,
                   &open1[ctx->off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
#endif
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,open2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,open2,OUTER_Z2_DECOMP_BASE);
  polxvec_copy(cnst->b,openrhs,1);

  cnst = &st->cnst[linear_idx];
  init_sparsecnst_view(cnst,r,
                       LNP22_USE_NATIVE_FOLD_BINARY ? 8 : 10,1,0,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,ctx->b_len,
                   &linear1[ctx->off[LNP22_B_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_B1,ctx->b_len,
                   &linear1[ctx->off[LNP22_B_BLOCK]],OUTER_BZ_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,
                   &linear1[ctx->off[LNP22_X_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,
                   &linear1[ctx->off[LNP22_X_BLOCK]],OUTER_XY_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_T0,tail_len,
                   &linear1[ctx->off[LNP22_RHAT_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_T1,tail_len,
                   &linear1[ctx->off[LNP22_RHAT_BLOCK]],OUTER_XY_DECOMP_BASE);
#if !LNP22_USE_NATIVE_FOLD_BINARY
  append_cnst_view(cnst,&i,OUT_RBIT0,binary_len,
                   &linear1[ctx->off[LNP22_RHAT_BITS_BLOCK]],1);
  append_cnst_view(cnst,&i,OUT_RBIT1,binary_len,
                   &linear1[ctx->off[LNP22_RHAT_BITS_BLOCK]],OUTER_XY_DECOMP_BASE);
#endif
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,linear2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,linear2,OUTER_Z2_DECOMP_BASE);
  polxvec_copy(cnst->b,linrhs,1);

  cnst = &st->cnst[quadratic_idx];
  init_sparsecnst_view(cnst,r,8,1,1,0);
  i = 0;
  append_cnst_view(cnst,&i,OUT_B0,ctx->b_len,quadB,1);
  append_cnst_view(cnst,&i,OUT_B1,ctx->b_len,quadB,OUTER_BZ_DECOMP_BASE);
  append_cnst_view(cnst,&i,OUT_X0,ctx->x_len,quadX,1);
  append_cnst_view(cnst,&i,OUT_X1,ctx->x_len,quadX,OUTER_XY_DECOMP_BASE);
#if LNP22_USE_NATIVE_FOLD_BINARY
  append_cnst_view(cnst,&i,OUT_RAUT0,binary_len,quadbinary,1);
  append_cnst_view(cnst,&i,OUT_RAUT1,binary_len,quadbinary,
                   OUTER_XY_DECOMP_BASE);
#else
  append_cnst_view(cnst,&i,OUT_RBIT0,binary_len,quadbinary,1);
  append_cnst_view(cnst,&i,OUT_RBIT1,binary_len,quadbinary,
                   OUTER_XY_DECOMP_BASE);
#endif
  append_cnst_view(cnst,&i,OUT_Z20,ctx->m2,quadz2,1);
  append_cnst_view(cnst,&i,OUT_Z21,ctx->m2,quadz2,OUTER_Z2_DECOMP_BASE);
  cnst->a->len = 2*OUTER_DIRECT_P_DIGITS +
                 (LNP22_USE_NATIVE_FOLD_BINARY ? 4 : 3);
  {
    size_t idx = 0, xa, pd;

    for(xa=0;xa<2;xa++) {
      const int64_t xscale = outer_pow_i64(OUTER_XY_DECOMP_BASE,(int)xa);
      for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++) {
        const int64_t pscale_i = outer_pow_i64(OUTER_DIRECT_P_BASE,(int)pd);
        cnst->a->rows[idx] = OUT_X0 + xa;
        cnst->a->cols[idx] = OUT_P0 + pd;
        set_direct_product_coeff(&cnst->a->coeffs[idx],xscale,pscale_i,-1);
        idx++;
      }
    }
#if LNP22_USE_NATIVE_FOLD_BINARY
    for(size_t ta=0;ta<2;ta++) {
      for(size_t ra=0;ra<2;ra++) {
        cnst->a->rows[idx] = OUT_T0 + ta;
        cnst->a->cols[idx] = OUT_RAUT0 + ra;
        set_direct_product_coeff(
            &cnst->a->coeffs[idx],
            outer_pow_i64(OUTER_XY_DECOMP_BASE,(int)ta),
            outer_pow_i64(OUTER_XY_DECOMP_BASE,(int)ra),1);
        polx_mul(&cnst->a->coeffs[idx],&cnst->a->coeffs[idx],
                 ctx->native_binary_scale);
        polx_refresh(&cnst->a->coeffs[idx]);
        idx++;
      }
    }
#else
    cnst->a->rows[idx] = OUT_RBIT0;
    cnst->a->cols[idx] = OUT_RBIT0;
    set_const_coeff(&cnst->a->coeffs[idx],1);
    idx++;
    cnst->a->rows[idx] = OUT_RBIT0;
    cnst->a->cols[idx] = OUT_RBIT1;
    set_const_coeff(&cnst->a->coeffs[idx],2*OUTER_XY_DECOMP_BASE);
    idx++;
    cnst->a->rows[idx] = OUT_RBIT1;
    cnst->a->cols[idx] = OUT_RBIT1;
    set_const_coeff(&cnst->a->coeffs[idx],
                    outer_pow_i64(OUTER_XY_DECOMP_BASE,2));
#endif
  }
  polxvec_copy(cnst->b,quadrhs,1);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);

  free(openfold);
  free(openrhs);
  {
    polx *owned_directD = NULL;
    polx *owned[11];

    if(transfer_ctx_storage) {
      ctx->lincoeff = NULL;
      ctx->lincoeff2 = NULL;
      owned_directD = ctx->directD;
      ctx->directD = NULL;
      bitweight = ctx->xbit_weight;
      ctx->xbit_weight = NULL;
    } else {
      linear1 = NULL;
      linear2 = NULL;
    }
    owned[0] = open1;
    owned[1] = open2;
    owned[2] = linear1;
    owned[3] = linear2;
    owned[4] = quadB;
    owned[5] = quadX;
    owned[6] = quadbinary;
    owned[7] = quadz2;
    owned[8] = pscale;
    owned[9] = owned_directD;
    owned[10] = transfer_ctx_storage ? bitweight : NULL;
    attach_outer_frontend_storage(st,owned);
  }
}

static void build_selector_lnp22_outer_statement_public(
    prncplstmnt *st, selector_lnp22_ctx *ctx,
    const polx *tA, const polx *w, const polx *vlin,
    const polx *tq, const polx *vq, const polx *c,
    const polx *a3, const polx *ufold_public, uint64_t betasq) {
  derive_direct_chain_fold(ctx,tA);
  compute_direct_public_chain_coeffs(ctx);
  build_direct_outer_statement(st,ctx,tA,w,vlin,tq,vq,c,a3,ufold_public,
                               betasq,0);
}

static void prepare_selector_lnp22_outer_direct(
    prncplstmnt *st, witness *wt, uint64_t seed_offset,
    outer_instance *public_transcript) {
  outer_instance inst = {};
  struct timespec prepare_start, inner_done, prepare_done;
  const size_t r = OUT_R;
  size_t n[OUT_R], i;
  size_t tail_len, binary_len;
  size_t p_len;
  uint64_t betasq = 0;
  polx *p = NULL;
  polx *raut = NULL;
  __int128 *pcoeffs = NULL;
  int64_t *rautcoeffs = NULL;
  poly *pdst[OUTER_DIRECT_P_DIGITS];

  clock_gettime(CLOCK_MONOTONIC,&prepare_start);
  if(generate_inner_instance(&inst,seed_offset)) {
    fprintf(stderr,"ERROR: failed to generate direct inner LNP22 instance\n");
    exit(1);
  }
  clock_gettime(CLOCK_MONOTONIC,&inner_done);
  printf("Selector outer inner generation: batch=%d, rejection rounds=%zu\n",
         LNP22_BATCH_ROUNDS,inst.rejection_rounds + 1);

  tail_len = inst.ctx.rhat_len + inst.ctx.umc_eua_len +
             inst.ctx.umc_eub_len + inst.ctx.umc_r_len +
             inst.ctx.umc_va_len + inst.ctx.umc_vb_len;
  binary_len = LNP22_USE_NATIVE_FOLD_BINARY ? inst.ctx.rhat_len :
                                              inst.ctx.rhat_bits_len;
  p_len = inst.ctx.x_len;
  n[OUT_B0] = inst.ctx.b_len;
  n[OUT_B1] = inst.ctx.b_len;
  n[OUT_X0] = inst.ctx.x_len;
  n[OUT_X1] = inst.ctx.x_len;
  n[OUT_P0] = p_len;
  n[OUT_P1] = p_len;
  n[OUT_P2] = p_len;
  n[OUT_P3] = p_len;
  n[OUT_P4] = p_len;
  n[OUT_T0] = tail_len;
  n[OUT_T1] = tail_len;
#if LNP22_USE_NATIVE_FOLD_BINARY
  n[OUT_RAUT0] = binary_len;
  n[OUT_RAUT1] = binary_len;
#else
  n[OUT_RBIT0] = binary_len;
  n[OUT_RBIT1] = binary_len;
#endif
  n[OUT_Z20] = inst.ctx.m2;
  n[OUT_Z21] = inst.ctx.m2;
  init_witness_raw(wt,r,n);

  p = alloc_polxvec(p_len);
  pcoeffs = _malloc(p_len*N*sizeof(*pcoeffs));
  compute_direct_outer_p(p,&inst.ctx,inst.z1);
  polxvec_centered_to_i128(pcoeffs,p,p_len);
#if LNP22_USE_NATIVE_FOLD_BINARY
  raut = alloc_polxvec(binary_len);
  rautcoeffs = _malloc(binary_len*N*sizeof(*rautcoeffs));
  for(i=0;i<binary_len;i++)
    polx_sigmam1(&raut[i],
                 &inst.z1[inst.ctx.off[LNP22_RHAT_BLOCK] + i]);
  polxvec_centered_to_i64(rautcoeffs,raut,binary_len);
#endif
#if DIRECT_OUTER_VALIDATE
  {
    polx lhs[1], tmp[1], rhs[1];

    polxvec_sprod(lhs,inst.ctx.direct_lin_B,
                  &inst.z1[inst.ctx.off[LNP22_B_BLOCK]],inst.ctx.b_len);
    polx_refresh(lhs);
    polx_mul(lhs,lhs,inst.c);
    polxvec_setzero(tmp,1);
    for(i=0;i<p_len;i++) {
      polx product[1];

      polx_mul(product,
               &inst.z1[inst.ctx.off[LNP22_X_BLOCK] + i],
               &p[i]);
      polx_add(tmp,tmp,product);
    }
    polx_sub(lhs,lhs,tmp);
    polxvec_sprod(tmp,inst.ctx.xbit_weight,
                  &inst.z1[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len);
    polx_mul(tmp,inst.c,tmp);
    polx_sub(lhs,lhs,tmp);
#if LNP22_USE_NATIVE_FOLD_BINARY
    rhat_native_binary_verifier_value(
        tmp,&inst.ctx,&inst.z1[inst.ctx.off[LNP22_RHAT_BLOCK]],inst.c);
    polx_add(lhs,lhs,tmp);
    polxvec_sprod(tmp,inst.ctx.native_binary_lin1,
                  &inst.z1[inst.ctx.off[LNP22_X_BLOCK]],inst.ctx.x_len);
    polx_mul(tmp,inst.c,tmp);
    polx_add(lhs,lhs,tmp);
    polxvec_sprod(tmp,inst.ctx.native_binary_lin2,inst.z2,inst.ctx.m2);
    polx_mul(tmp,inst.c,tmp);
    polx_add(lhs,lhs,tmp);
#else
    rhat_bits_verifier_value(
        tmp,&inst.ctx,&inst.z1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],inst.c);
    polx_add(lhs,lhs,tmp);
#endif
    inner_product_a3(tmp,inst.a3,inst.z2,inst.ctx.m2);
    polx_add(lhs,lhs,tmp);
    polx_mul(rhs,inst.c,inst.tq);
    polx_add(rhs,rhs,inst.vq);
#if LNP22_USE_NATIVE_FOLD_BINARY
    polx_mul(tmp,inst.c,inst.c);
    polx_mul(tmp,tmp,inst.uquad_public);
    polx_add(rhs,rhs,tmp);
#endif
    polx_sub(lhs,lhs,rhs);
    if(!polx_iszero(lhs)) {
      fprintf(stderr,"ERROR: direct outer conceptual quadratic relation failed\n");
      exit(1);
    }
  }
#endif
  for(i=0;i<OUTER_DIRECT_P_DIGITS;i++)
    pdst[i] = wt->s[OUT_P0+i];

  if(decompose_segment_two_digits(
         wt->s[OUT_B0],wt->s[OUT_B1],
         &inst.z1coeffs[inst.ctx.off[LNP22_B_BLOCK]*N],
         inst.ctx.b_len,OUTER_BZ_DECOMP_BASE,"B") ||
     decompose_direct_x_two_digits(
         wt->s[OUT_X0],wt->s[OUT_X1],
         &inst.z1coeffs[inst.ctx.off[LNP22_X_BLOCK]*N],
         &inst.ctx,OUTER_XY_DECOMP_BASE) ||
     decompose_segment_five_digits(
         pdst,pcoeffs,p_len,OUTER_DIRECT_P_BASE,"P") ||
     decompose_segment_two_digits(
         wt->s[OUT_T0],wt->s[OUT_T1],
         &inst.z1coeffs[inst.ctx.off[LNP22_RHAT_BLOCK]*N],
         tail_len,OUTER_XY_DECOMP_BASE,"tail") ||
#if LNP22_USE_NATIVE_FOLD_BINARY
     decompose_segment_two_digits(
         wt->s[OUT_RAUT0],wt->s[OUT_RAUT1],rautcoeffs,
         binary_len,OUTER_XY_DECOMP_BASE,"RAUT") ||
#else
     decompose_segment_two_digits(
         wt->s[OUT_RBIT0],wt->s[OUT_RBIT1],
         &inst.z1coeffs[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]*N],
         binary_len,OUTER_XY_DECOMP_BASE,"RBIT") ||
#endif
     decompose_segment_two_digits(
         wt->s[OUT_Z20],wt->s[OUT_Z21],inst.z2coeffs,
         inst.ctx.m2,OUTER_Z2_DECOMP_BASE,"z2")) {
    fprintf(stderr,"ERROR: failed to decompose direct outer witness\n");
    exit(1);
  }

  for(i=0;i<r;i++) {
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    betasq += wt->normsq[i];
  }
  betasq += 16*N;

  build_direct_outer_statement(st,&inst.ctx,inst.tA,inst.w,inst.vlin,
                               inst.tq,inst.vq,inst.c,inst.a3,
                               inst.ufold_public,betasq,1);

#if defined(DIRECT_OUTER_DEBUG) && DIRECT_OUTER_DEBUG
  {
    const sparsecnst *qcnst = &st->cnst[inst.ctx.x_len + 2];
    polx *sx[OUT_R] = {0};
    polx eval[1], matrix[1], part[1], expected[1], diff[1];
    size_t j;

    sparsecnst_eval(eval,qcnst,sx,wt);
    {
      size_t xa, pd, idx = 0;
      int coeffs_ok = 1;

      for(xa=0;xa<2;xa++) {
        for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++) {
          polx left[1], right[1], want[1], delta[1];
          set_const_coeff(left,outer_pow_i64(OUTER_XY_DECOMP_BASE,(int)xa));
          set_const_coeff(right,outer_pow_i64(OUTER_DIRECT_P_BASE,(int)pd));
          polx_mul(want,left,right);
          polx_neg(want,want);
          polx_sub(delta,&qcnst->a->coeffs[idx],want);
          coeffs_ok &= polx_iszero(delta);
          idx++;
        }
      }
      fprintf(stderr,"DEBUG direct outer: XP coefficients=%s\n",
              coeffs_ok ? "ok" : "mismatch");
    }
    {
      size_t xa, pd, idx = 0;

      for(xa=0;xa<2;xa++) {
        for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++) {
          polx got[1], want[1], delta[1];
          polx *scaled_x = alloc_polxvec(inst.ctx.x_len);
          polx *scaled_p = alloc_polxvec(inst.ctx.x_len);
          polxvec_sprod(got,sx[OUT_X0+xa],sx[OUT_P0+pd],
                        inst.ctx.x_len);
          polx_mul(got,got,&qcnst->a->coeffs[idx]);
          polxvec_scale(scaled_x,sx[OUT_X0+xa],inst.ctx.x_len,
                        outer_pow_i64(OUTER_XY_DECOMP_BASE,(int)xa));
          polxvec_scale(scaled_p,sx[OUT_P0+pd],inst.ctx.x_len,
                        outer_pow_i64(OUTER_DIRECT_P_BASE,(int)pd));
          polxvec_sprod(want,scaled_x,scaled_p,inst.ctx.x_len);
          polx_neg(want,want);
          polx_sub(delta,got,want);
          if(!polx_iszero(delta))
            fprintf(stderr,"DEBUG direct outer: pair xa=%zu pd=%zu mismatch\n",
                    xa,pd);
          free(scaled_x);
          free(scaled_p);
          idx++;
        }
      }
    }
    {
      polx *recon = alloc_polxvec(inst.ctx.x_len);
      polx *tmpvec = alloc_polxvec(inst.ctx.x_len);
      size_t pd;

      polxvec_copy(recon,sx[OUT_X0],inst.ctx.x_len);
      polxvec_scale(tmpvec,sx[OUT_X1],inst.ctx.x_len,
                    OUTER_XY_DECOMP_BASE);
      polxvec_add(recon,recon,tmpvec,inst.ctx.x_len);
      fprintf(stderr,"DEBUG direct outer: X reconstruction=%s\n",
              polxvec_equal_slow(
                  recon,&inst.z1[inst.ctx.off[LNP22_X_BLOCK]],
                  inst.ctx.x_len) ? "ok" : "mismatch");

      polxvec_setzero(recon,inst.ctx.x_len);
      for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++) {
        polxvec_scale(tmpvec,sx[OUT_P0+pd],inst.ctx.x_len,
                      outer_pow_i64(OUTER_DIRECT_P_BASE,(int)pd));
        polxvec_add(recon,recon,tmpvec,inst.ctx.x_len);
      }
      fprintf(stderr,"DEBUG direct outer: P reconstruction=%s\n",
              polxvec_equal_slow(recon,p,inst.ctx.x_len) ?
              "ok" : "mismatch");
      free(recon);
      free(tmpvec);
    }
    polxvec_setzero(matrix,1);
    for(j=0;j<2*OUTER_DIRECT_P_DIGITS;j++) {
      polxvec_sprod(part,sx[qcnst->a->rows[j]],sx[qcnst->a->cols[j]],
                    MIN(wt->n[qcnst->a->rows[j]],wt->n[qcnst->a->cols[j]]));
      polx_mul(part,part,&qcnst->a->coeffs[j]);
      polx_add(matrix,matrix,part);
    }
    polxvec_sprod(expected,&inst.z1[inst.ctx.off[LNP22_X_BLOCK]],
                  p,inst.ctx.x_len);
    polx_neg(expected,expected);
    polx_sub(diff,matrix,expected);
    fprintf(stderr,"DEBUG direct outer: XP matrix=%s\n",
            polx_iszero(diff) ? "ok" : "mismatch");
    {
      polx *rx = alloc_polxvec(inst.ctx.x_len);
      polx *rp = alloc_polxvec(inst.ctx.x_len);
      polx *tv = alloc_polxvec(inst.ctx.x_len);
      size_t pd;

      polxvec_copy(rx,sx[OUT_X0],inst.ctx.x_len);
      polxvec_scale(tv,sx[OUT_X1],inst.ctx.x_len,OUTER_XY_DECOMP_BASE);
      polxvec_add(rx,rx,tv,inst.ctx.x_len);
      polxvec_setzero(rp,inst.ctx.x_len);
      for(pd=0;pd<OUTER_DIRECT_P_DIGITS;pd++) {
        polxvec_scale(tv,sx[OUT_P0+pd],inst.ctx.x_len,
                      outer_pow_i64(OUTER_DIRECT_P_BASE,(int)pd));
        polxvec_add(rp,rp,tv,inst.ctx.x_len);
      }
      polxvec_sprod(expected,rx,rp,inst.ctx.x_len);
      polx_neg(expected,expected);
      polx_sub(diff,matrix,expected);
      fprintf(stderr,"DEBUG direct outer: XP reconstructed dot=%s\n",
              polx_iszero(diff) ? "ok" : "mismatch");
      free(rx);
      free(rp);
      free(tv);
    }

    polxvec_setzero(matrix,1);
    for(j=2*OUTER_DIRECT_P_DIGITS;j<qcnst->a->len;j++) {
      polxvec_sprod(part,sx[qcnst->a->rows[j]],sx[qcnst->a->cols[j]],
                    MIN(wt->n[qcnst->a->rows[j]],wt->n[qcnst->a->cols[j]]));
      polx_mul(part,part,&qcnst->a->coeffs[j]);
      polx_add(matrix,matrix,part);
    }
#if LNP22_USE_NATIVE_FOLD_BINARY
    polxvec_sprod(expected,&inst.z1[inst.ctx.off[LNP22_RHAT_BLOCK]],
                  raut,inst.ctx.rhat_len);
    polx_mul(expected,inst.ctx.native_binary_scale,expected);
#else
    polxvec_sprod(expected,
                  &inst.z1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],
                  &inst.z1[inst.ctx.off[LNP22_RHAT_BITS_BLOCK]],
                  inst.ctx.rhat_bits_len);
#endif
    polx_sub(diff,matrix,expected);
    fprintf(stderr,"DEBUG direct outer: bit matrix=%s\n",
            polx_iszero(diff) ? "ok" : "mismatch");
    free(sx[0]);
  }
#endif

  free(p);
  free(pcoeffs);
  free(raut);
  free(rautcoeffs);
  if(public_transcript != NULL) {
    *public_transcript = inst;
    memset(&inst,0,sizeof(inst));
  } else {
    free_outer_instance(&inst);
  }
  clock_gettime(CLOCK_MONOTONIC,&prepare_done);
  printf("Selector outer prepare breakdown: inner-generation=%.4fs, outer-frontend=%.4fs\n",
         elapsed_sec(&prepare_start,&inner_done),
         elapsed_sec(&inner_done,&prepare_done));
}

static void prepare_selector_lnp22_outer(prncplstmnt *st, witness *wt,
                                         uint64_t seed_offset) {
  prepare_selector_lnp22_outer_direct(st,wt,seed_offset,NULL);
}
