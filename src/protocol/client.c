#define _POSIX_C_SOURCE 200809L

#include <limits.h>
#include <sys/random.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "src/core/labrador.h"

#include "common.h"
#include "nizk1.h"
#include "wire.h"

#define LEOPARD_BATCH_MAX 2
#define LEOPARD_BUNDLE_FLAG (UINT64_C(1) << 63)

typedef struct {
  uint64_t *masks;
  size_t count;
  size_t cap;
} leopard_mask_list;

static double elapsed_sec_ts(const struct timespec *start,
                             const struct timespec *end) {
  return (double)(end->tv_sec - start->tv_sec) +
         1e-9 * (double)(end->tv_nsec - start->tv_nsec);
}

static size_t file_size_bytes(const char *path) {
  struct stat st;

  if(stat(path,&st) != 0)
    return 0;
  return (size_t)st.st_size;
}

static unsigned int packed_bits_u64(leopard_ucoeff span) {
  unsigned int bits = 0;

  while(span != 0) {
    bits++;
    span >>= 1;
  }
  return bits;
}

static size_t packed_payload_bytes(size_t count, unsigned int bits) {
  if(bits == 0 || count == 0)
    return 0;
  return (count * (size_t)bits + 7) / 8;
}

static size_t packed_i64_wire_bytes(const leopard_coeff *v, size_t count) {
  leopard_coeff minv;
  leopard_coeff maxv;
  unsigned int bits;
  size_t i;

  if(count == 0)
    return 0;
  minv = v[0];
  maxv = v[0];
  for(i=1;i<count;i++) {
    if(v[i] < minv)
      minv = v[i];
    if(v[i] > maxv)
      maxv = v[i];
  }
  bits = packed_bits_u64((leopard_ucoeff)maxv - (leopard_ucoeff)minv);
  return LEOPARD_COEFF_WIRE_BYTES + sizeof(uint8_t) + packed_payload_bytes(count,bits);
}

static int build_component_path(char *out, size_t out_len,
                                const char *dir, const char *prefix,
                                size_t server_idx, size_t component_idx) {
  return snprintf(out,out_len,"%s/%s_%zu_%zu.bin",dir,prefix,server_idx,component_idx) <
         (int)out_len ? 0 : 1;
}

static int build_bundle_path(char *out, size_t out_len,
                             const char *dir, const char *prefix,
                             size_t server_idx, size_t bundle_idx) {
  return snprintf(out,out_len,"%s/%s_%zu_%zu.bin",dir,prefix,server_idx,bundle_idx) <
         (int)out_len ? 0 : 1;
}

static int leopard_mask_list_push(leopard_mask_list *list, uint64_t mask) {
  uint64_t *tmp;

  if(list->count == list->cap) {
    size_t new_cap = (list->cap == 0) ? 8 : 2 * list->cap;

    tmp = realloc(list->masks,new_cap * sizeof(*tmp));
    if(tmp == NULL)
      return 1;
    list->masks = tmp;
    list->cap = new_cap;
  }
  list->masks[list->count++] = mask;
  return 0;
}

static void leopard_mask_list_clear(leopard_mask_list *list) {
  free(list->masks);
  memset(list,0,sizeof(*list));
}

static int leopard_popcount64(uint64_t x) {
  int c = 0;

  while(x != 0) {
    x &= x - 1;
    c++;
  }
  return c;
}

static int leopard_max_bit_index(uint64_t mask) {
  int idx = -1;

  while(mask != 0) {
    idx++;
    mask >>= 1;
  }
  return idx;
}

static int leopard_min_bit_index(uint64_t mask) {
  int idx = 0;

  if(mask == 0)
    return -1;
  while((mask & 1u) == 0) {
    idx++;
    mask >>= 1;
  }
  return idx;
}

static uint64_t leopard_bundle_id(const uint64_t *masks, size_t count) {
  uint64_t h = 0x9e3779b97f4a7c15ULL ^ (uint64_t)count;
  size_t i;

  if(count == 1)
    return masks[0];
  for(i=0;i<count;i++) {
    h ^= masks[i] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h *= 0xff51afd7ed558ccdULL;
  }
  return LEOPARD_BUNDLE_FLAG | (h & ~LEOPARD_BUNDLE_FLAG);
}

static int leopard_component_owner(uint64_t component_mask, uint64_t q_mask) {
  if(NIZK2_DOUBLE_COMMITMENT) {
    uint64_t holders=component_mask & q_mask;
    return holders ? 63-__builtin_clzll(holders) : -1;
  }
  return leopard_min_bit_index(component_mask & q_mask);
}

static int leopard_enum_masks_rec(size_t nservers, size_t choose,
                                  size_t next, size_t chosen,
                                  uint64_t current, int require_idx,
                                  leopard_mask_list *out) {
  size_t i;

  if(chosen == choose) {
    if(require_idx < 0 || ((current >> require_idx) & 1u) != 0)
      return leopard_mask_list_push(out,current);
    return 0;
  }
  if(next >= nservers)
    return 0;
  for(i=next;i + (choose - chosen) <= nservers;i++) {
    if(leopard_enum_masks_rec(nservers,choose,i + 1,chosen + 1,
                              current | (UINT64_C(1) << i),
                              require_idx,out) != 0) {
      return 1;
    }
  }
  return 0;
}

static int leopard_collect_global_masks(size_t nservers, size_t threshold,
                                        leopard_mask_list *out) {
  if(nservers == 0 || threshold == 0 || threshold > nservers || nservers > 63)
    return 1;
  memset(out,0,sizeof(*out));
  return leopard_enum_masks_rec(nservers,nservers - threshold + 1,0,0,0,-1,out);
}

static int leopard_collect_local_masks(size_t nservers, size_t threshold,
                                       size_t server_idx,
                                       leopard_mask_list *out) {
  if(server_idx >= nservers)
    return 1;
  memset(out,0,sizeof(*out));
  return leopard_enum_masks_rec(nservers,nservers - threshold + 1,0,0,0,
                                (int)server_idx,out);
}

static int leopard_find_local_component_index(size_t nservers, size_t threshold,
                                              size_t server_idx,
                                              uint64_t component_mask,
                                              size_t *out_idx) {
  leopard_mask_list local = {0};
  size_t i;
  int ret = 1;

  if(leopard_collect_local_masks(nservers,threshold,server_idx,&local) != 0)
    return 1;
  for(i=0;i<local.count;i++) {
    if(local.masks[i] == component_mask) {
      *out_idx = i;
      ret = 0;
      break;
    }
  }
  leopard_mask_list_clear(&local);
  return ret;
}

static int leopard_validate_q_mask(size_t nservers, size_t threshold,
                                   uint64_t q_mask) {
  const uint64_t valid_mask =
      (nservers >= 64) ? UINT64_MAX : ((UINT64_C(1) << nservers) - UINT64_C(1));

  if(nservers == 0 || threshold == 0 || threshold > nservers || nservers > 63)
    return 1;
  if((q_mask & ~valid_mask) != 0)
    return 1;
  if((size_t)leopard_popcount64(q_mask) != threshold)
    return 1;
  return 0;
}

static int leopard_add_i64(leopard_coeff **out, size_t *out_count,
                           const leopard_coeff *in, size_t count) {
  size_t i;

  if(*out == NULL) {
    *out = malloc(count * sizeof(**out));
    if(*out == NULL)
      return 1;
    memcpy(*out,in,count * sizeof(**out));
    *out_count = count;
    return 0;
  }
  if(*out_count != count)
    return 1;
  for(i=0;i<count;i++)
    (*out)[i] += in[i];
  return 0;
}

static int leopard_add_prep(leopard_prep_data *out, const leopard_prep_data *in,
                            uint64_t bundle_id) {
  size_t count = in->slot_count * in->vk_rows * LEOPARD_RING_DEG;
  size_t out_count = out->slot_count * out->vk_rows * LEOPARD_RING_DEG;

  if(out->vk_coeffs != NULL &&
     (out->slot_count != in->slot_count || out->vk_rows != in->vk_rows)) {
    return 1;
  }
  if(leopard_add_i64(&out->vk_coeffs,&out_count,in->vk_coeffs,count) != 0)
    return 1;
  out->slot_count = in->slot_count;
  out->vk_rows = in->vk_rows;
  out->component_mask = bundle_id;
  return 0;
}

static int leopard_load_finalize_bundle(leopard_prep_data *out_prep,
                                        const char *prep_dir,
                                        size_t server_idx,
                                        const size_t *local_indices,
                                        const uint64_t *masks,
                                        size_t count) {
  uint64_t bid = leopard_bundle_id(masks,count);
  size_t bi;
  char path[PATH_MAX];

  memset(out_prep,0,sizeof(*out_prep));
  for(bi=0;bi<count;bi++) {
    leopard_prep_data prep = {0};

    if(build_component_path(path,sizeof(path),prep_dir,"prep",server_idx,local_indices[bi]) != 0 ||
       leopard_read_prep_file(path,&prep) != 0 ||
       prep.component_mask != masks[bi] ||
       leopard_add_prep(out_prep,&prep,bid) != 0) {
      leopard_prep_data_clear(&prep);
      return 1;
    }
    leopard_prep_data_clear(&prep);
  }
  return 0;
}

static int leopard_parse_q_spec(const char *spec,
                                size_t nservers,
                                size_t threshold,
                                uint64_t *out_q_mask) {
  uint64_t q_mask = 0;

  if(spec == NULL || spec[0] == '\0')
    return 1;
  if(strchr(spec,',') != NULL) {
    char *copy = NULL;
    char *tok;
    char *saveptr = NULL;

    copy = strdup(spec);
    if(copy == NULL)
      return 1;
    for(tok = strtok_r(copy,",",&saveptr); tok != NULL;
        tok = strtok_r(NULL,",",&saveptr)) {
      char *endptr = NULL;
      unsigned long long idx;

      if(*tok == '\0') {
        free(copy);
        return 1;
      }
      idx = strtoull(tok,&endptr,10);
      if(endptr == tok || *endptr != '\0' || idx >= nservers) {
        free(copy);
        return 1;
      }
      q_mask |= UINT64_C(1) << idx;
    }
    free(copy);
  } else {
    char *endptr = NULL;

    q_mask = strtoull(spec,&endptr,0);
    if(endptr == spec || *endptr != '\0')
      return 1;
  }

  if(leopard_validate_q_mask(nservers,threshold,q_mask) != 0)
    return 1;
  *out_q_mask = q_mask;
  return 0;
}

static void usage(const char *prog) {
  fprintf(stderr,
          "usage:\n"
          "  %s setup <pp.bin>\n"
          "  %s preproc-client <pp.bin> <slot_count> <preq.bin> <pre_state.bin>\n"
          "  %s request <pp.bin> <pre_state.bin> <nservers> <threshold> <q_spec> <tag_hex> <x_hex> <req.bin> <st.bin>\n"
          "  %s verify-request <pp.bin> <tag_hex> <req.bin> [cached-preq.bin]\n"
          "  %s finalize <st.bin> <prep_dir> <rep_dir>\n"
          "  %s eval <pp.bin> <tag_hex> <x_hex> <master.bin>\n",
          prog,prog,prog,prog,prog,prog);
}

int main(int argc, char **argv) {
  struct timespec t0, t1;
  double sec;
  int rc = 1;

  if(argc < 2) {
    usage(argv[0]);
    return 1;
  }

  const int epoch_request = !strcmp(argv[1],"epoch-request");
  const int epoch_verify = !strcmp(argv[1],"epoch-verify");
  if(!strcmp(argv[1],"epoch-setup") || !strcmp(argv[1],"setup")) {
    leopard_pp pp;

    if(argc != 3) {
      usage(argv[0]);
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    leopard_pp_init(&pp);
    const char *bench_seed=getenv("LEOPARD_BENCH_SEED_OFFSET");
    if(bench_seed) {
      char *end;
      unsigned long long offset=strtoull(bench_seed,&end,10);
      if(!*bench_seed || *end) return 1;
      /* Deterministic benchmark fixture; this harness also derives prover
         sampling from it and is not a production randomness design. */
      for(size_t i=0;i<16;i++) pp.nizk1_seed[i]=(uint8_t)i;
      unsigned int carry=0;
      for(size_t i=0;i<16;i++) {
        const unsigned int sum=pp.nizk1_seed[i]+(i<8 ? (offset>>(8*i))&255 : 0)+carry;
        pp.nizk1_seed[i]=(uint8_t)sum;
        carry=sum>>8;
      }
    }
    if(!strcmp(argv[1],"epoch-setup") && getrandom(&pp,sizeof(pp),0) != (ssize_t)sizeof(pp)) return 1;
    if(leopard_write_pp_file(argv[2],&pp) != 0)
      return 1;
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    printf("ok=1\n");
    printf("time_sec=%.6f\n",sec);
    printf("pp_bytes=%zu\n",sizeof(pp));
    return 0;
  }

  if(strcmp(argv[1],"preproc-client") == 0) {
    leopard_pp pp;
    leopard_preq_data preq = {0};
    leopard_pre_state_data pre_st = {0};
    size_t slot_count;

    if(argc != 6) {
      usage(argv[0]);
      return 1;
    }
    slot_count = (size_t)strtoull(argv[3],NULL,10);
    if(slot_count == 0 || leopard_read_pp_file(argv[2],&pp) != 0)
      return 1;
    clock_gettime(CLOCK_MONOTONIC,&t0);
    rc = leopard_nizk1_preproc_client_data(&preq,&pre_st,&pp,slot_count);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0 &&
       leopard_write_preq_file(argv[4],&preq) == 0 &&
       leopard_write_pre_state_file(argv[5],&pre_st) == 0) {
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("slot_count=%zu\n",slot_count);
      printf("preq_bytes=%zu\n",file_size_bytes(argv[4]));
      printf("pre_state_bytes=%zu\n",file_size_bytes(argv[5]));
      rc = 0;
    }
    leopard_preq_data_clear(&preq);
    leopard_pre_state_data_clear(&pre_st);
    free_comkey();
    return rc;
  }

  if(epoch_request || !strcmp(argv[1],"request")) {
    leopard_pp pp;
    leopard_pre_state_data pre_st = {0};
    leopard_request_data req = {0};
    leopard_state_data st = {0};
    leopard_nizk1_fast_public_profile fast_prof = {0};
    uint8_t tag[LEOPARD_TAG_BYTES];
    uint8_t xbytes[LEOPARD_X_BYTES];
    size_t nservers, threshold;
    uint64_t q_mask;

    if(argc != 11) {
      usage(argv[0]);
      return 1;
    }
    nservers = (size_t)strtoull(argv[4],NULL,10);
    threshold = (size_t)strtoull(argv[5],NULL,10);
    if(leopard_read_pp_file(argv[2],&pp) != 0 ||
       (!epoch_request && leopard_read_pre_state_file(argv[3],&pre_st) != 0) ||
       leopard_parse_hex(tag,sizeof(tag),argv[7]) != 0 ||
       leopard_parse_hex(xbytes,sizeof(xbytes),argv[8]) != 0 ||
       leopard_parse_q_spec(argv[6],nservers,threshold,&q_mask) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    if(epoch_request) {
      uint8_t seed[16];
      leopard_preq_data scratch = {0};
      if(leopard_parse_hex(seed,16,argv[3])) return 1;
      leopard_client_epoch(seed,1);
      if(leopard_nizk1_preproc_client_data(&scratch,&pre_st,&pp,1)) return 1;
      leopard_preq_data_clear(&scratch);
    }
    rc = leopard_nizk1_request_prove_preprocessed_data(&req,&st,&pre_st,&pp,
                                                       tag,sizeof(tag),xbytes);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    req.nservers = nservers;
    req.threshold = threshold;
    req.q_mask = q_mask;
    req.cr_external = !epoch_request && getenv("LEOPARD_REQUEST_CACHED_PREPROC") &&
                      strcmp(getenv("LEOPARD_REQUEST_CACHED_PREPROC"),"1") == 0;
    st.nservers = nservers;
    st.threshold = threshold;
    st.q_mask = q_mask;
    if(rc == 0 && !epoch_request &&
       leopard_nizk1_request_fast_public_preprocessed_data(&fast_prof,&pp,&pre_st,
                                                           req.slot_j,tag,sizeof(tag),
                                                           xbytes,&req) != 0) {
      leopard_request_data_clear(&req);
      leopard_state_data_clear(&st);
      leopard_pre_state_data_clear(&pre_st);
      free_comkey();
      return 1;
    }
    if(rc == 0 &&
       leopard_write_request_file(argv[9],&req) == 0 &&
       leopard_write_state_file(argv[10],&st) == 0 &&
       (epoch_request || leopard_write_pre_state_file(argv[3],&pre_st) == 0)) {
      const size_t req_file_bytes = file_size_bytes(argv[9]);
      const size_t req_public_bytes = req_file_bytes - req.proof_len;
      const size_t req_cx_bytes = packed_i64_wire_bytes(
          req.cx_coeffs,req.cx_rows * req.cx_cols * LEOPARD_RING_DEG);
      const leopard_nizk1_request_profile *prof = leopard_nizk1_last_request_profile();

      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("slot_j=%zu\n",req.slot_j);
      printf("q_mask=0x%llx\n",(unsigned long long)q_mask);
      printf("req_bytes=%zu\n",req_file_bytes);
      printf("req_cached_preprocessing=%d\n",req.cr_external);
      printf("req_cr_omitted_bytes=%zu\n",req.cr_external ?
          packed_i64_wire_bytes(req.cr_coeffs,req.cr_rows*LEOPARD_RING_DEG) : 0);
      printf("req_cr_wire_bytes=%zu\n",req.cr_external ? 0 :
          packed_i64_wire_bytes(req.cr_coeffs,req.cr_rows*LEOPARD_RING_DEG));
      printf("req_public_bytes=%zu\n",req_public_bytes);
      printf("req_public_cx_bytes=%zu\n",req_cx_bytes);
      printf("req_public_other_bytes=%zu\n",
             req_public_bytes >= req_cx_bytes ? req_public_bytes - req_cx_bytes : 0);
      printf("req_proof_bytes=%zu\n",req.proof_len);
      printf("req_protocol_input_sec=%.6f\n",prof->protocol_input_sec);
      printf("req_protocol_transcript_sec=%.6f\n",prof->protocol_transcript_sec);
      printf("req_protocol_selector_material_sec=%.6f\n",prof->protocol_selector_material_sec);
      printf("req_protocol_slot_load_sec=%.6f\n",prof->protocol_slot_load_sec);
      printf("req_protocol_hash_eval_sec=%.6f\n",prof->protocol_hash_eval_sec);
      printf("req_protocol_cx_sec=%.6f\n",prof->protocol_cx_sec);
      printf("req_protocol_state_sec=%.6f\n",prof->protocol_state_sec);
      printf("req_protocol_public_export_sec=%.6f\n",prof->protocol_public_export_sec);
      printf("req_fast_public_total_sec=%.6f\n",fast_prof.total_sec);
      printf("req_fast_public_selector_material_sec=%.6f\n",fast_prof.selector_material_sec);
      printf("req_fast_public_slot_load_sec=%.6f\n",fast_prof.slot_load_sec);
      printf("req_fast_public_chain_sec=%.6f\n",fast_prof.chain_sec);
      printf("req_fast_public_cx_sec=%.6f\n",fast_prof.cx_sec);
      printf("req_nizk_prepare_sec=%.6f\n",prof->nizk_prepare_sec);
      printf("req_nizk_inner_setup_sec=%.6f\n",prof->nizk_inner_setup_sec);
      printf("req_inner_instance_sec=%.6f\n",prof->inner_instance_sec);
      printf("req_inner_prepare_sec=%.6f\n",prof->inner_prepare_sec);
      printf("req_inner_prove_sec=%.6f\n",prof->inner_prove_sec);
      printf("req_inner_selfcheck_sec=%.6f\n",prof->inner_selfcheck_sec);
      printf("req_inner_rejection_rounds=%zu\n",prof->inner_rejection_rounds);
      printf("req_inner_batch_rounds=%zu\n",prof->inner_batch_rounds);
      printf("req_a1open_cache=%d\n",prof->a1open_cache);
      printf("req_outer_witness_sec=%.6f\n",prof->outer_witness_sec);
      printf("req_outer_statement_sec=%.6f\n",prof->outer_statement_sec);
      printf("req_outer_stmt_fold_sec=%.6f\n",prof->outer_stmt_fold_sec);
      printf("req_outer_stmt_quad_sec=%.6f\n",prof->outer_stmt_quad_sec);
      printf("req_outer_stmt_cnst1_sec=%.6f\n",prof->outer_stmt_cnst1_sec);
      printf("req_outer_stmt_cnst2_sec=%.6f\n",prof->outer_stmt_cnst2_sec);
      printf("req_outer_stmt_cnst3_sec=%.6f\n",prof->outer_stmt_cnst3_sec);
      printf("req_outer_stmt_hash_sec=%.6f\n",prof->outer_stmt_hash_sec);
      printf("req_outer_prove_sec=%.6f\n",prof->outer_prove_sec);
      printf("req_serialize_sec=%.6f\n",prof->serialize_sec);
      printf("req_export_sec=%.6f\n",prof->export_sec);
      rc = 0;
    }
    leopard_request_data_clear(&req);
    leopard_state_data_clear(&st);
    leopard_pre_state_data_clear(&pre_st);
    free_comkey();
    return rc;
  }

  if(epoch_verify || !strcmp(argv[1],"verify-request")) {
    leopard_pp pp;
    leopard_request_data req = {0};
    leopard_preq_data cached = {0};
    uint8_t tag[LEOPARD_TAG_BYTES];

    if(argc != 5 && argc != 6) {
      usage(argv[0]);
      return 1;
    }
    if(leopard_read_pp_file(argv[2],&pp) != 0 ||
       leopard_parse_hex(tag,sizeof(tag),argv[3]) != 0 ||
       leopard_read_request_file(argv[4],&req) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    /* Cache lookup is timed and never replaces the proof verification. */
    if(epoch_verify) {
      uint8_t seed[16];
      if(argc != 6 || req.nservers != 3 || req.threshold != 2 || req.cr_external || req.slot_j ||
         leopard_parse_hex(seed,16,argv[5]) || memcmp(req.ar_seed,seed,16)) return 1;
      leopard_client_epoch(seed,0);
    } else if((req.cr_external && argc != 6) ||
       (argc == 6 && (leopard_read_preq_file(argv[5],&cached) != 0 ||
                     leopard_request_attach_preprocessing(&req,&cached) != 0))) {
      leopard_preq_data_clear(&cached);
      leopard_request_data_clear(&req);
      fprintf(stderr,"missing or mismatched preprocessed commitment\n");
      return 1;
    }
    leopard_preq_data_clear(&cached);
    rc = leopard_nizk1_request_verify_data(&pp,tag,sizeof(tag),&req);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    {
      const leopard_nizk1_verify_profile *prof = leopard_nizk1_last_verify_profile();

      printf("ok=%d\n",rc == 0 ? 1 : 0);
      printf("time_sec=%.6f\n",sec);
      printf("verify_public_ctx_sec=%.6f\n",prof->public_ctx_sec);
      printf("verify_a1open_cache=%d\n",prof->a1open_cache);
      printf("verify_fs_seed_sec=%.6f\n",prof->fs_seed_sec);
      printf("verify_deserialize_sec=%.6f\n",prof->deserialize_sec);
      printf("verify_outer_statement_sec=%.6f\n",prof->outer_statement_sec);
      printf("verify_outer_stmt_fold_sec=%.6f\n",prof->outer_stmt_fold_sec);
      printf("verify_outer_stmt_quad_sec=%.6f\n",prof->outer_stmt_quad_sec);
      printf("verify_outer_stmt_cnst1_sec=%.6f\n",prof->outer_stmt_cnst1_sec);
      printf("verify_outer_stmt_cnst2_sec=%.6f\n",prof->outer_stmt_cnst2_sec);
      printf("verify_outer_stmt_cnst3_sec=%.6f\n",prof->outer_stmt_cnst3_sec);
      printf("verify_outer_stmt_hash_sec=%.6f\n",prof->outer_stmt_hash_sec);
      printf("verify_outer_verify_sec=%.6f\n",prof->outer_verify_sec);
    }
    leopard_request_data_clear(&req);
    free_comkey();
    return rc;
  }

  const int epoch_final=!strcmp(argv[1],"epoch-finalize");
  const int epoch_eval=!strcmp(argv[1],"epoch-eval");
  if(epoch_final || epoch_eval) {
    uint8_t pk[32], secret[16];
    if(argc!=(epoch_eval ? 8 : 6) || leopard_parse_hex(pk,32,argv[epoch_eval ? 6 : 5])) return 1;
    if(epoch_eval) {
      FILE *f=fopen(argv[7],"rb");
      if(!f || fread(secret,1,16,f)!=16 || fgetc(f)!=EOF) return 1;
      fclose(f);
    }
    leopard_output_context(pk,epoch_eval ? secret : NULL);
  }
  if(epoch_final || !strcmp(argv[1],"finalize")) {
    leopard_state_data st = {0};
    leopard_reply_data rep = {0};
    leopard_prep_data prep = {0};
    leopard_mask_list local = {0};
    uint8_t y[LEOPARD_OUTPUT_BYTES];
    char y_hex[2 * LEOPARD_OUTPUT_BYTES + 1];
    char rep_path[PATH_MAX];
    leopard_coeff *agg_ux = NULL;
    leopard_coeff *agg_vk = NULL;
    size_t agg_ux_rows = 0;
    size_t agg_vk_rows = 0;
    size_t total_rep_bytes = 0;
    size_t bundle_total = 0;
    size_t component_total = 0;
    size_t i, li;

    if(argc != (epoch_final ? 6 : 5)) {
      usage(argv[0]);
      return 1;
    }
    if(leopard_read_state_file(argv[2],&st) != 0 ||
       leopard_validate_q_mask(st.nservers,st.threshold,st.q_mask) != 0) {
      return 1;
    }

    rc = 0;
    for(i=0;i<st.nservers;i++) {
      if(((st.q_mask >> i) & 1u) == 0)
        continue;
      if(leopard_collect_local_masks(st.nservers,st.threshold,i,&local) != 0) {
        rc = 1;
        goto finalize_cleanup;
      }
      size_t bundle_indices[LEOPARD_BATCH_MAX];
      uint64_t bundle_masks[LEOPARD_BATCH_MAX];
      size_t bundle_count = 0;
      size_t bundle_idx = 0;

      for(li=0;li<local.count;li++) {
        if(leopard_component_owner(local.masks[li],st.q_mask) != (int)i)
          continue;
        bundle_indices[bundle_count] = li;
        bundle_masks[bundle_count] = local.masks[li];
        bundle_count++;
        component_total++;
        if(bundle_count != LEOPARD_BATCH_MAX)
          continue;

        if(leopard_load_finalize_bundle(&prep,argv[3],i,bundle_indices,
                                        bundle_masks,bundle_count) != 0 ||
           build_bundle_path(rep_path,sizeof(rep_path),argv[4],"repb",i,bundle_idx) != 0 ||
           leopard_read_reply_file(rep_path,&rep) != 0 ||
           rep.component_mask != leopard_bundle_id(bundle_masks,bundle_count) ||
           st.slot_j >= prep.slot_count) {
          rc = 1;
          goto finalize_cleanup;
        }
        total_rep_bytes += file_size_bytes(rep_path);
        if(agg_ux == NULL) {
          agg_ux_rows = rep.ux_rows;
          agg_vk_rows = prep.vk_rows;
          agg_ux = calloc(agg_ux_rows * LEOPARD_RING_DEG,sizeof(*agg_ux));
          agg_vk = calloc(agg_vk_rows * LEOPARD_RING_DEG,sizeof(*agg_vk));
          if(agg_ux == NULL || agg_vk == NULL) {
            rc = 1;
            goto finalize_cleanup;
          }
        }
        if(rep.ux_rows != agg_ux_rows || prep.vk_rows != agg_vk_rows) {
          rc = 1;
          goto finalize_cleanup;
        }
        for(size_t j=0;j<agg_ux_rows * LEOPARD_RING_DEG;j++)
          agg_ux[j] += rep.ux_coeffs[j];
        {
          const leopard_coeff *vk_slot = &prep.vk_coeffs[st.slot_j * prep.vk_rows * LEOPARD_RING_DEG];
          for(size_t j=0;j<agg_vk_rows * LEOPARD_RING_DEG;j++)
            agg_vk[j] += vk_slot[j];
        }
        leopard_reply_data_clear(&rep);
        leopard_prep_data_clear(&prep);
        bundle_idx++;
        bundle_total++;
        bundle_count = 0;
      }
      if(bundle_count != 0) {
        if(leopard_load_finalize_bundle(&prep,argv[3],i,bundle_indices,
                                        bundle_masks,bundle_count) != 0 ||
           build_bundle_path(rep_path,sizeof(rep_path),argv[4],"repb",i,bundle_idx) != 0 ||
           leopard_read_reply_file(rep_path,&rep) != 0 ||
           rep.component_mask != leopard_bundle_id(bundle_masks,bundle_count) ||
           st.slot_j >= prep.slot_count) {
          rc = 1;
          goto finalize_cleanup;
        }
        total_rep_bytes += file_size_bytes(rep_path);
        if(agg_ux == NULL) {
          agg_ux_rows = rep.ux_rows;
          agg_vk_rows = prep.vk_rows;
          agg_ux = calloc(agg_ux_rows * LEOPARD_RING_DEG,sizeof(*agg_ux));
          agg_vk = calloc(agg_vk_rows * LEOPARD_RING_DEG,sizeof(*agg_vk));
          if(agg_ux == NULL || agg_vk == NULL) {
            rc = 1;
            goto finalize_cleanup;
          }
        }
        if(rep.ux_rows != agg_ux_rows || prep.vk_rows != agg_vk_rows) {
          rc = 1;
          goto finalize_cleanup;
        }
        for(size_t j=0;j<agg_ux_rows * LEOPARD_RING_DEG;j++)
          agg_ux[j] += rep.ux_coeffs[j];
        {
          const leopard_coeff *vk_slot = &prep.vk_coeffs[st.slot_j * prep.vk_rows * LEOPARD_RING_DEG];
          for(size_t j=0;j<agg_vk_rows * LEOPARD_RING_DEG;j++)
            agg_vk[j] += vk_slot[j];
        }
        leopard_reply_data_clear(&rep);
        leopard_prep_data_clear(&prep);
        bundle_idx++;
        bundle_total++;
      }
      leopard_mask_list_clear(&local);
    }

    if(agg_ux == NULL || agg_vk == NULL) {
      rc = 1;
      goto finalize_cleanup;
    }

    clock_gettime(CLOCK_MONOTONIC,&t0);
    rc = leopard_nizk1_finalize_raw_output_data(y,&st,agg_ux,agg_ux_rows,
                                                agg_vk,agg_vk_rows);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0) {
      leopard_hex_encode(y_hex,y,sizeof(y));
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("rep_total_bytes=%zu\n",total_rep_bytes);
      printf("bundle_count=%zu\n",bundle_total);
      printf("component_count=%zu\n",component_total);
      printf("nservers=%zu\n",st.nservers);
      printf("threshold=%zu\n",st.threshold);
      printf("q_mask=0x%llx\n",(unsigned long long)st.q_mask);
      printf("y_hex=%s\n",y_hex);
    }

finalize_cleanup:
    free(agg_ux);
    free(agg_vk);
    leopard_reply_data_clear(&rep);
    leopard_prep_data_clear(&prep);
    leopard_mask_list_clear(&local);
    leopard_state_data_clear(&st);
    free_comkey();
    return rc;
  }


  if(epoch_eval || !strcmp(argv[1],"eval")) {
    leopard_pp pp;
    uint8_t tag[LEOPARD_TAG_BYTES];
    uint8_t xbytes[LEOPARD_X_BYTES];
    uint8_t y[LEOPARD_OUTPUT_BYTES];
    char y_hex[2 * LEOPARD_OUTPUT_BYTES + 1];
    leopard_coeff *k_coeffs = NULL;
    size_t k_rows = 0;

    if(argc != (epoch_eval ? 8 : 6)) {
      usage(argv[0]);
      return 1;
    }
    if(leopard_read_pp_file(argv[2],&pp) != 0 ||
       leopard_parse_hex(tag,sizeof(tag),argv[3]) != 0 ||
       leopard_parse_hex(xbytes,sizeof(xbytes),argv[4]) != 0 ||
       leopard_read_master_key_file(argv[5],&k_coeffs,&k_rows) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    rc = leopard_nizk1_eval(y,&pp,tag,sizeof(tag),xbytes,k_coeffs,k_rows);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0) {
      leopard_hex_encode(y_hex,y,sizeof(y));
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("y_hex=%s\n",y_hex);
    }
    free(k_coeffs);
    free_comkey();
    return rc;
  }

  usage(argv[0]);
  free_comkey();
  return 1;
}
