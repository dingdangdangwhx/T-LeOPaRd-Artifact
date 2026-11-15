#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "malloc.h"
#include "jlproj.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"
#include "pack.h"

#ifndef PACK_VERIFY_PROFILE
#define PACK_VERIFY_PROFILE 0
#endif
#ifndef PACK_LEVEL_PROFILE
#define PACK_LEVEL_PROFILE 0
#endif
#ifndef PACK_PREFLIGHT_STOP
#define PACK_PREFLIGHT_STOP 1
#endif
#ifndef PACK_MIN_GAIN_KB
#define PACK_MIN_GAIN_KB 0.0
#endif
#ifndef PACK_RELATION_DEBUG
#define PACK_RELATION_DEBUG 0
#endif

static double pack_elapsed_sec(const struct timespec *a, const struct timespec *b) {
  return (double)(b->tv_sec - a->tv_sec) + 1e-9*(double)(b->tv_nsec - a->tv_nsec);
}

typedef struct {
  double initial_prove_s;
  double initial_size_s;
  double recursive_prove_s;
  double recursive_size_s;
  size_t recursive_levels;
} pack_profile;

static double proof_size_kb(const proof *pi) {
  double s;
  const comparams *cpp = pi->cpp;

  s  = (log2(sqrt(jlproj_normsq(pi->p)))-4+2.05)*256;
  s += (cpp->u1len+cpp->u2len+LIFTS)*N*LOGQ;
  s /= 8192;
  return s;
}

static double witness_size_kb(const witness *wt) {
  size_t i;
  double s = 0;

  for(i=0;i<wt->r;i++)
    s += (log2((double)wt->normsq[i]/(N*wt->n[i]))/2+2.05)*N*wt->n[i];
  s /= 8192;
  return s;
}

static int pack_output_statement_secure(const statement *st) {
  const comparams *cpp = st->cpp;
  const double inner_norm = 6*T*SLACK*
      ldexp(1,(cpp->f-1)*cpp->b)*sqrt(st->betasq);

#if PACK_LEVEL_PROFILE
  size_t need_inner=1,need_outer=1;
  while(need_inner<=32 && !sis_secure(need_inner,inner_norm)) need_inner++;
  while(need_outer<=32 && !sis_secure(need_outer,2*SLACK*sqrt(st->betasq))) need_outer++;
  printf("Pack actual norm: %.4g, rank=%zu/%zu required=%zu/%zu f=%zu b=%zu fg=%zu bg=%zu\n",
         sqrt(st->betasq),cpp->kappa,cpp->kappa1,need_inner,need_outer,cpp->f,cpp->b,cpp->fg,cpp->bg);
#endif

  if(!sis_secure(cpp->kappa,inner_norm))
    return 0;
  if(!st->tail && !sis_secure(cpp->kappa1,2*SLACK*sqrt(st->betasq)))
    return 0;
  return 1;
}

static int predict_candidate_size_kb(double *proof_kb, double *witness_kb,
                                     const statement *st, const witness *wt,
                                     int tail) {
  proof pi = {};
  const int quadratic = st->cpp->fg != 0;
  size_t i;
  size_t out_r = 0;
  size_t out_n = 0;
  size_t joined = 0;
  size_t out_m;
  size_t witness_len;
  uint64_t input_normsq = 0;
  double projected_norm;
  double bits;
  int ret;

  ret = init_proof(&pi,wt,quadratic,tail);
  if(ret)
    return ret;

  for(i=0;i<pi.r;i++) {
    joined += pi.n[i];
    if(pi.nu[i]) {
      out_r += pi.nu[i];
      out_n = MAX(out_n,(joined + pi.nu[i] - 1)/pi.nu[i]);
      joined = 0;
    }
  }
  out_m = out_r*pi.cpp->fu*pi.cpp->kappa;
  out_m += (pi.cpp->fu + pi.cpp->fg)*(out_r*out_r + out_r)/2;
  witness_len = pi.cpp->f*out_n + (tail ? 0 : out_m);

  if(witness_len == 0 || pi.normsq == 0) {
    free_proof(&pi);
    return 1;
  }
  bits = (log2((double)pi.normsq/(N*witness_len))/2 + 2.05)*N*witness_len;
  *witness_kb = bits/8192;

  for(i=0;i<wt->r;i++)
    input_normsq += wt->normsq[i];
  projected_norm = sqrt(256.0*(double)input_normsq);
  bits = (log2(projected_norm) - 4 + 2.05)*256;
  bits += (pi.cpp->u1len + pi.cpp->u2len + LIFTS)*N*LOGQ;
  *proof_kb = bits/8192;

  free_proof(&pi);
  return 0;
}

void free_composite(composite *p) {
  size_t i;

  free_witness(&p->owt);
  for(i=0;i<p->l;i++) {
    free_proof(p->pi[i]);
    free(p->pi[i]);
    p->pi[i] = NULL;
  }
  p->l = 0;
}

static int composite_prove(composite *p, statement *tst, witness *twt, double *twtsize, pack_profile *prof) {
  int ret;
  size_t i = 0;
  double pisize;
  clock_t t0;

  while(p->l < 16) {
#if PACK_PREFLIGHT_STOP || PACK_LEVEL_PROFILE
    {
      double predicted_proof = 0.0;
      double predicted_witness = 0.0;
      const int prediction_ok = !predict_candidate_size_kb(&predicted_proof,&predicted_witness,
                                                            &tst[i],&twt[i],0);
#if PACK_LEVEL_PROFILE
      if(prediction_ok)
        printf("Pack preflight[%zu]: proof=%.3f KB witness=%.3f KB total=%.3f KB current=%.3f KB tail=0\n",
               p->l,predicted_proof,predicted_witness,predicted_proof + predicted_witness,twtsize[i]);
#endif
      if(PACK_PREFLIGHT_STOP && prediction_ok &&
         twtsize[i] - predicted_proof - predicted_witness <= PACK_MIN_GAIN_KB) {
        p->pi[p->l] = _malloc(sizeof(proof));
        break;
      }
    }
#endif
    p->pi[p->l] = _malloc(sizeof(proof));
    t0 = clock();
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],0);
    if(prof)
      prof->recursive_prove_s += (double)(clock() - t0) / CLOCKS_PER_SEC;
    if(ret) return ret;
#if PACK_RELATION_DEBUG
    fprintf(stderr,
            "DEBUG Pack candidate[%zu] params: r=%zu f=%zu b=%zu fu=%zu bu=%zu fg=%zu bg=%zu kappa=%zu kappa1=%zu\n",
            p->l,p->pi[p->l]->r,p->pi[p->l]->cpp->f,p->pi[p->l]->cpp->b,
            p->pi[p->l]->cpp->fu,p->pi[p->l]->cpp->bu,
            p->pi[p->l]->cpp->fg,p->pi[p->l]->cpp->bg,
            p->pi[p->l]->cpp->kappa,p->pi[p->l]->cpp->kappa1);
    ret = verify(&tst[i^1],&twt[i^1]);
    fprintf(stderr,"DEBUG Pack candidate[%zu] relation: %s (%d)\n",
            p->l,ret ? "FAILED" : "ok",ret);
    if(ret)
      return 40 + (int)p->l;
#endif
    if(!pack_output_statement_secure(&tst[i^1])) {
#if PACK_LEVEL_PROFILE
      printf("Pack candidate[%zu]: rejected because actual norm is not SIS-secure\n",
             p->l);
#endif
      free_proof(p->pi[p->l]);
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }
    t0 = clock();
    pisize = proof_size_kb(p->pi[p->l]);
    twtsize[i^1] = witness_size_kb(&twt[i^1]);
    if(prof) {
      prof->recursive_size_s += (double)(clock() - t0) / CLOCKS_PER_SEC;
      prof->recursive_levels += 1;
    }
#if PACK_LEVEL_PROFILE
    printf("Pack candidate[%zu]: proof=%.3f KB next-witness=%.3f KB current-witness=%.3f KB gain=%.3f KB tail=0\n",
           p->l,pisize,twtsize[i^1],twtsize[i],twtsize[i] - pisize - twtsize[i^1]);
#endif
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }

    free_statement(&tst[i]);
    free_witness(&twt[i]);
    p->size += pisize;
    p->l += 1;
    i ^= 1;
  }

  if(p->l < 16) {
#if PACK_LEVEL_PROFILE
    {
      double predicted_proof = 0.0;
      double predicted_witness = 0.0;
      if(!predict_candidate_size_kb(&predicted_proof,&predicted_witness,&tst[i],&twt[i],1))
        printf("Pack preflight[%zu]: proof=%.3f KB witness=%.3f KB total=%.3f KB current=%.3f KB tail=1\n",
               p->l,predicted_proof,predicted_witness,predicted_proof + predicted_witness,twtsize[i]);
    }
#endif
    t0 = clock();
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],1);
    if(prof)
      prof->recursive_prove_s += (double)(clock() - t0) / CLOCKS_PER_SEC;
    if(ret) return ret;
#if PACK_RELATION_DEBUG
    ret = verify(&tst[i^1],&twt[i^1]);
    fprintf(stderr,"DEBUG Pack tail candidate[%zu] relation: %s (%d)\n",
            p->l,ret ? "FAILED" : "ok",ret);
    if(ret)
      return 60 + (int)p->l;
#endif
    const int tail_secure = pack_output_statement_secure(&tst[i^1]);
    t0 = clock();
    pisize = proof_size_kb(p->pi[p->l]);
    twtsize[i^1] = witness_size_kb(&twt[i^1]);
    if(prof) {
      prof->recursive_size_s += (double)(clock() - t0) / CLOCKS_PER_SEC;
      prof->recursive_levels += 1;
    }
#if PACK_LEVEL_PROFILE
    printf("Pack candidate[%zu]: proof=%.3f KB next-witness=%.3f KB current-witness=%.3f KB gain=%.3f KB tail=1\n",
           p->l,pisize,twtsize[i^1],twtsize[i],twtsize[i] - pisize - twtsize[i^1]);
#endif
    if(!tail_secure || pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
    }
    else {
      p->size += pisize;
      p->l += 1;
      i ^= 1;
    }

    free_statement(&tst[i^1]);
    free_witness(&twt[i^1]);
  }

  free_statement(&tst[i]);
  p->owt = twt[i];
  p->size += twtsize[i];
  return 0;
}

int composite_prove_principle(composite *p, const prncplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;
  clock_t t0;
  pack_profile prof = {};

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  t0 = clock();
  ret = principle_prove(tst,twt,p->pi[p->l],st,wt,0);
  prof.initial_prove_s = (double)(clock() - t0) / CLOCKS_PER_SEC;
  if(ret)
    goto err;
  t0 = clock();
  p->size = proof_size_kb(p->pi[p->l]);
  twtsize[0] = witness_size_kb(twt);
  prof.initial_size_s = (double)(clock() - t0) / CLOCKS_PER_SEC;
  p->l += 1;
  ret = composite_prove(p,tst,twt,twtsize,&prof);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Chihuahua Pack members: %zu\n",p->l);
  printf("Chihuahua Pack size: %.2f KB\n",p->size);
  printf("Chihuahua Pack time breakdown: principle=%.4fs, recursive=%.4fs, size-only=%.4fs, recursive levels=%zu\n",
         prof.initial_prove_s,
         prof.recursive_prove_s,
         prof.initial_size_s + prof.recursive_size_s,
         prof.recursive_levels);
  printf("Chihuahua Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_simple(composite *p, commitment *com, const smplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  ret = simple_prove(tst,twt,p->pi[p->l],com,st,wt,0);
  if(ret)
    goto err;
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize,NULL);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Dachshund Pack members: %zu\n",p->l);
  printf("Dachshund Pack size: %.2f KB\n",p->size);
  printf("Dachshund Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_polcom(composite *p, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y) {
  int ret;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  polcom_eval(&twt[1],ppi,ctx,x,y);
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;
  p->size = print_polcomprf_pp(ppi);
  print_prncplstmnt_pp(tst0);
  twtsize[1] = print_witness_pp(&twt[1]);

  ret = principle_prove(tst,twt,p->pi[p->l],tst0,&twt[1],0);
  if(ret) {
    ret += 10;
    goto err;
  }
  p->size += print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  free_prncplstmnt(tst0);
  free_witness(&twt[1]);
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize,NULL);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Greyhound Pack members: %zu\n",p->l);
  printf("Greyhound Pack size: %.2f KB\n",p->size);
  printf("Greyhound Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

static int composite_verify(const composite *p, statement *tst) {
  size_t i,j;
  int ret;
#if PACK_VERIFY_PROFILE
  struct timespec t0, t1;
#endif

  i = 0;
  for(j=1;j<p->l;j++) {
#if PACK_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
    ret = reduce(&tst[i^1],p->pi[j],&tst[i]);
#if PACK_VERIFY_PROFILE
    clock_gettime(CLOCK_MONOTONIC,&t1);
    printf("Pack verify profile: reduce[%zu]=%.6fs\n",j,pack_elapsed_sec(&t0,&t1));
#endif
    free_statement(&tst[i]);
    i ^= 1;
    if(ret)  // projection too long or commitments not secure (1/2/3)
      return ret + 4*j;
  }

#if PACK_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t0);
#endif
  ret = verify(&tst[i],&p->owt);
#if PACK_VERIFY_PROFILE
  clock_gettime(CLOCK_MONOTONIC,&t1);
  printf("Pack verify profile: final verify=%.6fs\n",pack_elapsed_sec(&t0,&t1));
#endif
  if(ret)
    return ret + 100;

  return 0;
}

int composite_verify_principle(const composite *p, const prncplstmnt *st) {
  int ret;
  statement tst[2] = {};
  struct timespec t_start, t_end;

  clock_gettime(CLOCK_MONOTONIC,&t_start);
#if PACK_VERIFY_PROFILE
  {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC,&t0);
    ret = principle_reduce(tst,p->pi[0],st);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    printf("Pack verify profile: principle_reduce=%.6fs\n",pack_elapsed_sec(&t0,&t1));
  }
#else
  ret = principle_reduce(tst,p->pi[0],st);
#endif
  if(ret)  // projection too long or commitments not secure (1/2/3)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  clock_gettime(CLOCK_MONOTONIC,&t_end);
  printf("Chihuahua Pack verification time: %.4fs\n",pack_elapsed_sec(&t_start,&t_end));
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_simple(const composite *p, const commitment *com, const smplstmnt *st) {
  int ret;
  statement tst[2] = {};
  clock_t t;

  t = clock();
  ret = simple_reduce(tst,p->pi[0],com,st);
  if(ret)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;
  printf("Dachshund Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_polcom(const composite *p, const polcomprf *ppi) {
  int ret;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  clock_t t;

  t = clock();
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;

  ret = principle_reduce(tst,p->pi[0],tst0);
  free_prncplstmnt(tst0);
  if(ret) {
    ret += 10;
    goto err;
  }

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = clock() - t;
  printf("Greyhound Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}
