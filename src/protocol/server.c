#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../lazer.h"

#include "common.h"
#include "nizk2.h"
#include "wire.h"

#define LEOPARD_BATCH_MAX 2
#if NIZK2_SMALL_COMMITMENT && NIZK2_AGGREGATE_MAX != 2
#error "The integrated server requires aggregate-capable small-commitment parameters"
#endif
#define LEOPARD_BUNDLE_FLAG (UINT64_C(1) << 63)

#ifndef LEOPARD_PARALLEL_BLINDEVAL
#define LEOPARD_PARALLEL_BLINDEVAL 0
#endif

#ifndef LEOPARD_PARALLEL_VERIFY_REPLIES
#define LEOPARD_PARALLEL_VERIFY_REPLIES 1
#endif

typedef struct {
  uint64_t *masks;
  size_t count;
  size_t cap;
} leopard_mask_list;

typedef struct {
  size_t server_idx;
  size_t bundle_idx;
  size_t count;
  size_t local_indices[LEOPARD_BATCH_MAX];
  uint64_t masks[LEOPARD_BATCH_MAX];
  size_t rep_total_bytes;
  size_t rep_public_total_bytes;
  size_t rep_proof_total_bytes;
  double build_total;
  double protocol_eval_total;
  double nizk_build_total;
  double prove_total;
  double export_total;
  double verify_import_total;
  double verify_build_total;
  double verify_core_total;
  double verify_total;
  int rc;
} leopard_blindeval_job;

typedef struct {
  leopard_blindeval_job *jobs;
  size_t count;
  size_t cap;
} leopard_blindeval_job_list;

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

static int ensure_dir(const char *path) {
  if(mkdir(path,0777) == 0 || errno == EEXIST)
    return 0;
  return 1;
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

static int leopard_blindeval_job_list_push(leopard_blindeval_job_list *list,
                                           size_t server_idx,
                                           size_t bundle_idx,
                                           const size_t *local_indices,
                                           const uint64_t *masks,
                                           size_t count) {
  leopard_blindeval_job *tmp;
  size_t i;

  if(count > LEOPARD_BATCH_MAX)
    return 1;
  if(list->count == list->cap) {
    size_t new_cap = (list->cap == 0) ? 8 : 2 * list->cap;

    tmp = realloc(list->jobs,new_cap * sizeof(*tmp));
    if(tmp == NULL)
      return 1;
    list->jobs = tmp;
    list->cap = new_cap;
  }
  leopard_blindeval_job *job = &list->jobs[list->count++];
  memset(job,0,sizeof(*job));
  job->server_idx = server_idx;
  job->bundle_idx = bundle_idx;
  job->count = count;
  for(i=0;i<count;i++) {
    job->local_indices[i] = local_indices[i];
    job->masks[i] = masks[i];
  }
  return 0;
}

static void leopard_blindeval_job_list_clear(leopard_blindeval_job_list *list) {
  free(list->jobs);
  memset(list,0,sizeof(*list));
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

static int leopard_popcount64(uint64_t x) {
  int c = 0;

  while(x != 0) {
    x &= x - 1;
    c++;
  }
  return c;
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

static leopard_coeff leopard_center_q_i128(__int128 x) {
  const __int128 q = (__int128)LEOPARD_MOD_Q;
  const __int128 half = q / 2;
  __int128 y = x % q;

  if(y > half)
    y -= q;
  if(y < -half)
    y += q;
  return (leopard_coeff)y;
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
    (*out)[i] = leopard_center_q_i128((__int128)(*out)[i] + (__int128)in[i]);
  return 0;
}

static int leopard_add_pk(leopard_pk_data *out, const leopard_pk_data *in,
                          uint64_t bundle_id) {
  size_t count = in->ck_rows * LEOPARD_RING_DEG;
  size_t out_count = out->ck_rows * LEOPARD_RING_DEG;

  if(leopard_add_i64(&out->ck_coeffs,&out_count,in->ck_coeffs,count) != 0)
    return 1;
#if NIZK2_SMALL_COMMITMENT
  for(size_t i=0;i<count;i++) {
    leopard_coeff x = out->ck_coeffs[i] % NIZK2_COM_Q;
    if(x > NIZK2_COM_Q/2) x -= NIZK2_COM_Q;
    if(x < -NIZK2_COM_Q/2) x += NIZK2_COM_Q;
    out->ck_coeffs[i] = x;
  }
#endif
  out->ck_rows = in->ck_rows;
  out->component_mask = bundle_id;
  return 0;
}

static int leopard_add_sk_vec(leopard_coeff **out, size_t *out_rows,
                              const leopard_coeff *in, size_t in_rows) {
  size_t out_count = *out_rows * LEOPARD_RING_DEG;

  if(leopard_add_i64(out,&out_count,in,in_rows * LEOPARD_RING_DEG) != 0)
    return 1;
  *out_rows = in_rows;
  return 0;
}

static int leopard_add_sk(leopard_sk_data *out, const leopard_sk_data *in,
                          uint64_t bundle_id) {
  if(leopard_add_sk_vec(&out->r_coeffs,&out->r_rows,in->r_coeffs,in->r_rows) != 0 ||
     leopard_add_sk_vec(&out->e_coeffs,&out->e_rows,in->e_coeffs,in->e_rows) != 0 ||
     leopard_add_sk_vec(&out->ei_coeffs,&out->ei_rows,in->ei_coeffs,in->ei_rows) != 0 ||
     leopard_add_sk_vec(&out->k_coeffs,&out->k_rows,in->k_coeffs,in->k_rows) != 0) {
    return 1;
  }
  out->component_mask = bundle_id;
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

static int leopard_add_server_pre(leopard_server_pre_state_data *out,
                                  const leopard_server_pre_state_data *in,
                                  uint64_t bundle_id) {
  size_t vk_count = in->slot_count * in->vk_rows * LEOPARD_RING_DEG;
  size_t es_count = in->slot_count * in->es_rows * LEOPARD_RING_DEG;
  size_t out_vk_count = out->slot_count * out->vk_rows * LEOPARD_RING_DEG;
  size_t out_es_count = out->slot_count * out->es_rows * LEOPARD_RING_DEG;

  if(out->vk_coeffs != NULL &&
     (out->slot_count != in->slot_count || out->vk_rows != in->vk_rows ||
      out->es_rows != in->es_rows)) {
    return 1;
  }
  if(leopard_add_i64(&out->vk_coeffs,&out_vk_count,in->vk_coeffs,vk_count) != 0 ||
     leopard_add_i64(&out->es_coeffs,&out_es_count,in->es_coeffs,es_count) != 0) {
    return 1;
  }
  out->slot_count = in->slot_count;
  out->vk_rows = in->vk_rows;
  out->es_rows = in->es_rows;
  out->component_mask = bundle_id;
  return 0;
}

static int leopard_load_blindeval_bundle(leopard_pk_data *out_pk,
                                         leopard_sk_data *out_sk,
                                         leopard_server_pre_state_data *out_pre,
                                         const char *pk_dir,
                                         const char *sk_dir,
                                         const char *sps_dir,
                                         size_t server_idx,
                                         const size_t *local_indices,
                                         const uint64_t *masks,
                                         size_t count) {
  uint64_t bid = leopard_bundle_id(masks,count);
  size_t bi;
  char path[PATH_MAX];

  memset(out_pk,0,sizeof(*out_pk));
  memset(out_sk,0,sizeof(*out_sk));
  memset(out_pre,0,sizeof(*out_pre));
  for(bi=0;bi<count;bi++) {
    leopard_pk_data pk = {0};
    leopard_sk_data sk = {0};
    leopard_server_pre_state_data pre = {0};

    if(build_component_path(path,sizeof(path),pk_dir,"pk",server_idx,local_indices[bi]) != 0 ||
       leopard_read_pk_file(path,&pk) != 0 ||
       pk.component_mask != masks[bi] ||
       build_component_path(path,sizeof(path),sk_dir,"sk",server_idx,local_indices[bi]) != 0 ||
       leopard_read_sk_file(path,&sk) != 0 ||
       sk.component_mask != masks[bi] ||
       build_component_path(path,sizeof(path),sps_dir,"sps",server_idx,local_indices[bi]) != 0 ||
       leopard_read_server_pre_state_file(path,&pre) != 0 ||
       pre.component_mask != masks[bi] ||
       leopard_add_pk(out_pk,&pk,bid) != 0 ||
       leopard_add_sk(out_sk,&sk,bid) != 0 ||
       leopard_add_server_pre(out_pre,&pre,bid) != 0) {
      leopard_pk_data_clear(&pk);
      leopard_sk_data_clear(&sk);
      leopard_server_pre_state_data_clear(&pre);
      return 1;
    }
    leopard_pk_data_clear(&pk);
    leopard_sk_data_clear(&sk);
    leopard_server_pre_state_data_clear(&pre);
  }
  return 0;
}

static int leopard_load_verify_bundle(leopard_pk_data *out_pk,
                                      leopard_prep_data *out_prep,
                                      const char *pk_dir,
                                      const char *prep_dir,
                                      size_t server_idx,
                                      const size_t *local_indices,
                                      const uint64_t *masks,
                                      size_t count) {
  uint64_t bid = leopard_bundle_id(masks,count);
  size_t bi;
  char path[PATH_MAX];

  memset(out_pk,0,sizeof(*out_pk));
  memset(out_prep,0,sizeof(*out_prep));
  for(bi=0;bi<count;bi++) {
    leopard_pk_data pk = {0};
    leopard_prep_data prep = {0};

    if(build_component_path(path,sizeof(path),pk_dir,"pk",server_idx,local_indices[bi]) != 0 ||
       leopard_read_pk_file(path,&pk) != 0 ||
       pk.component_mask != masks[bi] ||
       build_component_path(path,sizeof(path),prep_dir,"prep",server_idx,local_indices[bi]) != 0 ||
       leopard_read_prep_file(path,&prep) != 0 ||
       prep.component_mask != masks[bi] ||
       leopard_add_pk(out_pk,&pk,bid) != 0 ||
       leopard_add_prep(out_prep,&prep,bid) != 0) {
      leopard_pk_data_clear(&pk);
      leopard_prep_data_clear(&prep);
      return 1;
    }
    leopard_pk_data_clear(&pk);
    leopard_prep_data_clear(&prep);
  }
  return 0;
}

static int leopard_run_blindeval_bundle(
    const leopard_pp *pp,
    const char *pk_dir,
    const char *sk_dir,
    const char *sps_dir,
    const char *rep_dir,
    const uint8_t *tag,
    size_t tag_len,
    const leopard_request_data *req,
    size_t server_idx,
    size_t bundle_idx,
    const size_t *local_indices,
    const uint64_t *masks,
    size_t count,
    size_t *rep_total_bytes,
    size_t *rep_public_total_bytes,
    size_t *rep_proof_total_bytes,
    double *build_total,
    double *protocol_eval_total,
    double *nizk_build_total,
    double *prove_total,
    double *export_total) {
  leopard_pk_data pk = {0};
  leopard_sk_data sk = {0};
  leopard_server_pre_state_data pre_st = {0};
  leopard_reply_data rep = {0};
  const leopard_nizk2_blindeval_profile *prof;
  char rep_path[PATH_MAX];
  size_t rep_bytes;
  int ret = 1;

  if(leopard_load_blindeval_bundle(&pk,&sk,&pre_st,pk_dir,sk_dir,sps_dir,
                                   server_idx,local_indices,masks,count) != 0) {
    goto end;
  }
  if(build_bundle_path(rep_path,sizeof(rep_path),rep_dir,"repb",server_idx,bundle_idx) != 0) {
    goto end;
  }
  if(leopard_nizk2_blindeval_preprocessed_data(&rep,pp,&pk,&sk,&pre_st,
                                               tag,tag_len,req) != 0) {
    goto end;
  }
  if(rep.component_mask != leopard_bundle_id(masks,count)) {
    goto end;
  }
  if(leopard_write_reply_file(rep_path,&rep) != 0) {
    goto end;
  }

  rep_bytes = file_size_bytes(rep_path);
  *rep_total_bytes += rep_bytes;
  *rep_proof_total_bytes += rep.proof_len;
  *rep_public_total_bytes += rep_bytes - rep.proof_len;
  prof = leopard_nizk2_last_blindeval_profile();
  *protocol_eval_total += prof->protocol_eval_sec;
  *nizk_build_total += prof->nizk_build_sec;
  *build_total += prof->build_instance_sec;
  *prove_total += prof->prover_sec;
  *export_total += prof->export_sec;
  ret = 0;

end:
  leopard_pk_data_clear(&pk);
  leopard_sk_data_clear(&sk);
  leopard_server_pre_state_data_clear(&pre_st);
  leopard_reply_data_clear(&rep);
  return ret;
}

static int leopard_verify_reply_bundle(const leopard_pp *pp,
                                       const char *pk_dir,
                                       const char *prep_dir,
                                       const char *rep_dir,
                                       const uint8_t *tag,
                                       size_t tag_len,
                                       const leopard_request_data *req,
                                       size_t server_idx,
                                       size_t bundle_idx,
                                       const size_t *local_indices,
                                       const uint64_t *masks,
                                       size_t count,
                                       double *verify_import_total,
                                       double *verify_build_total,
                                       double *verify_core_total,
                                       double *verify_total) {
  leopard_pk_data pk = {0};
  leopard_prep_data prep = {0};
  leopard_reply_data rep = {0};
  const leopard_nizk2_verify_profile *prof;
  char rep_path[PATH_MAX];
  int ret = 1;

  if(leopard_load_verify_bundle(&pk,&prep,pk_dir,prep_dir,server_idx,
                                local_indices,masks,count) != 0 ||
     build_bundle_path(rep_path,sizeof(rep_path),rep_dir,"repb",server_idx,bundle_idx) != 0 ||
     leopard_read_reply_file(rep_path,&rep) != 0 ||
     rep.component_mask != leopard_bundle_id(masks,count)) {
    goto end;
  }
  if(leopard_nizk2_reply_verify_data(pp,&pk,&prep,tag,tag_len,req,&rep) != 0) {
    goto end;
  }
  prof = leopard_nizk2_last_verify_profile();
  *verify_import_total += prof->import_sec;
  *verify_build_total += prof->build_instance_sec;
  *verify_core_total += prof->verifier_sec;
  *verify_total += prof->total_sec;
  ret = 0;

end:
  leopard_pk_data_clear(&pk);
  leopard_prep_data_clear(&prep);
  leopard_reply_data_clear(&rep);
  return ret;
}

static void usage(const char *prog) {
  fprintf(stderr,
          "usage:\n"
          "  %s keygen <pp.bin> <nservers> <threshold> <pk_dir> <sk_dir> <master.bin>\n"
          "  %s preproc-server <pp.bin> <sk_dir> <nservers> <threshold> <preq.bin> <prep_dir> <server_pre_dir>\n"
          "  %s blindeval <pp.bin> <pk_dir> <sk_dir> <server_pre_dir> <tag_hex> <req.bin> <rep_dir>\n"
          "  %s verify-replies <pp.bin> <pk_dir> <prep_dir> <tag_hex> <req.bin> <rep_dir>\n",
          prog,prog,prog,prog);
}

int main(int argc, char **argv) {
  struct timespec t0, t1;
  double sec;
  int rc = 1;

  if(argc < 2) {
    usage(argv[0]);
    return 1;
  }

  lazer_init();

  if(strcmp(argv[1],"keygen") == 0) {
    leopard_pp pp;
    leopard_pk_data pk = {0}, sample_pk = {0};
    leopard_sk_data sk = {0}, sample_sk = {0};
    leopard_mask_list global = {0}, local = {0};
    leopard_coeff *master_k = NULL;
    size_t master_rows = 0;
    size_t nservers, threshold, gi, si, li;
    size_t pk_total_bytes = 0;
    double total_party_sec = 0.0;
    uint64_t nonce_base = 10000;
    char pathbuf[PATH_MAX];

    if(argc != 8) {
      usage(argv[0]);
      return 1;
    }
    nservers = (size_t)strtoull(argv[3],NULL,10);
    threshold = (size_t)strtoull(argv[4],NULL,10);
    if(nservers == 0 || threshold == 0 || threshold > nservers || nservers > 63 ||
       leopard_read_pp_file(argv[2],&pp) != 0 ||
       ensure_dir(argv[5]) != 0 ||
       ensure_dir(argv[6]) != 0 ||
       leopard_collect_global_masks(nservers,threshold,&global) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    for(gi=0;gi<global.count;gi++) {
      struct timespec component_start,component_end;
      clock_gettime(CLOCK_MONOTONIC,&component_start);
      if(leopard_nizk2_keygen_data(&sample_pk,&sample_sk,&pp) != 0) {
        rc = 1;
        break;
      }
      clock_gettime(CLOCK_MONOTONIC,&component_end);
      if(NIZK2_DOUBLE_COMMITMENT) {
        sample_pk.component_mask=global.masks[gi];
        sample_sk.component_mask=global.masks[gi];
        total_party_sec+=elapsed_sec_ts(&component_start,&component_end);
      }
      if(master_k == NULL) {
        master_rows = sample_sk.k_rows;
        master_k = calloc(master_rows * LEOPARD_RING_DEG,sizeof(*master_k));
        if(master_k == NULL) {
          rc = 1;
          break;
        }
      }
      for(li=0;li<sample_sk.k_rows * LEOPARD_RING_DEG;li++)
        master_k[li] = leopard_center_q_i128((__int128)master_k[li] + sample_sk.k_coeffs[li]);

      for(si=0;si<nservers;si++) {
        struct timespec p0, p1;

        if(((global.masks[gi] >> si) & 1u) == 0)
          continue;
        if(leopard_collect_local_masks(nservers,threshold,si,&local) != 0) {
          rc = 1;
          break;
        }
        for(li=0;li<local.count;li++) {
          if(local.masks[li] == global.masks[gi])
            break;
        }
        if(li == local.count) {
          leopard_mask_list_clear(&local);
          rc = 1;
          break;
        }
        clock_gettime(CLOCK_MONOTONIC,&p0);
        rc = NIZK2_DOUBLE_COMMITMENT ? 0 :
          leopard_nizk2_keygen_given_key_data(&pk,&sk,&pp,sample_sk.k_coeffs,
                                             global.masks[gi],nonce_base);
        clock_gettime(CLOCK_MONOTONIC,&p1);
        nonce_base += 16;
        leopard_mask_list_clear(&local);
        if(rc != 0)
          break;
        if(!NIZK2_DOUBLE_COMMITMENT) total_party_sec += elapsed_sec_ts(&p0,&p1);
        const leopard_pk_data *stored_pk=NIZK2_DOUBLE_COMMITMENT ? &sample_pk : &pk;
        const leopard_sk_data *stored_sk=NIZK2_DOUBLE_COMMITMENT ? &sample_sk : &sk;
        if(build_component_path(pathbuf,sizeof(pathbuf),argv[5],"pk",si,li) != 0 ||
           leopard_write_pk_file(pathbuf,stored_pk) != 0) {
          rc = 1;
          break;
        }
        /* Holder files replicate the same public component, sent only once. */
        if(!NIZK2_DOUBLE_COMMITMENT || si==(size_t)__builtin_ctzll(global.masks[gi]))
          pk_total_bytes += file_size_bytes(pathbuf);
        if(build_component_path(pathbuf,sizeof(pathbuf),argv[6],"sk",si,li) != 0 ||
           leopard_write_sk_file(pathbuf,stored_sk) != 0) {
          rc = 1;
          break;
        }
        leopard_pk_data_clear(&pk);
        leopard_sk_data_clear(&sk);
      }
      leopard_pk_data_clear(&sample_pk);
      leopard_sk_data_clear(&sample_sk);
      if(rc != 0)
        break;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0 &&
       leopard_write_master_key_file(argv[7],master_k,master_rows) == 0) {
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("nservers=%zu\n",nservers);
      printf("threshold=%zu\n",threshold);
      printf("pk_total_bytes=%zu\n",pk_total_bytes);
      printf("pk_avg_bytes=%zu\n",pk_total_bytes / nservers);
#if NIZK2_DOUBLE_COMMITMENT
      printf("commitment_mode=double\nkey_distribution=uniform_mod_q\n");
#elif NIZK2_SMALL_COMMITMENT
      printf("commitment_mode=small\ncommitment_modulus=%lld\ncommitment_rank=%d\n",
             (long long)NIZK2_COM_Q,NIZK2_COM_RANK);
#else
      printf("commitment_mode=large\n");
#endif
      printf("keygen_party_total_sec=%.6f\n",total_party_sec);
      printf("keygen_party_avg_sec=%.6f\n",total_party_sec / (double)nservers);
      rc = 0;
    }
    free(master_k);
    leopard_pk_data_clear(&pk);
    leopard_sk_data_clear(&sk);
    leopard_pk_data_clear(&sample_pk);
    leopard_sk_data_clear(&sample_sk);
    leopard_mask_list_clear(&global);
    leopard_mask_list_clear(&local);
    return rc;
  }

  if(strcmp(argv[1],"preproc-server") == 0) {
    leopard_pp pp;
    leopard_sk_data sk = {0};
    leopard_preq_data preq = {0};
    leopard_prep_data prep = {0};
    leopard_server_pre_state_data pre_st = {0};
    leopard_mask_list local = {0};
    size_t nservers, threshold, i, li;
    size_t prep_total_bytes = 0;
    double total_party_sec = 0.0;
    char sk_path[PATH_MAX];
    char prep_path[PATH_MAX];
    char sps_path[PATH_MAX];

    if(argc != 9) {
      usage(argv[0]);
      return 1;
    }
    nservers = (size_t)strtoull(argv[4],NULL,10);
    threshold = (size_t)strtoull(argv[5],NULL,10);
    if(nservers == 0 || threshold == 0 || threshold > nservers || nservers > 63 ||
       leopard_read_pp_file(argv[2],&pp) != 0 ||
       leopard_read_preq_file(argv[6],&preq) != 0 ||
       ensure_dir(argv[7]) != 0 ||
       ensure_dir(argv[8]) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    for(i=0;i<nservers;i++) {
      if(leopard_collect_local_masks(nservers,threshold,i,&local) != 0) {
        rc = 1;
        break;
      }
      for(li=0;li<local.count;li++) {
        struct timespec p0, p1;

        if(build_component_path(sk_path,sizeof(sk_path),argv[3],"sk",i,li) != 0 ||
           build_component_path(prep_path,sizeof(prep_path),argv[7],"prep",i,li) != 0 ||
           build_component_path(sps_path,sizeof(sps_path),argv[8],"sps",i,li) != 0 ||
           leopard_read_sk_file(sk_path,&sk) != 0) {
          rc = 1;
          break;
        }
        clock_gettime(CLOCK_MONOTONIC,&p0);
        rc = leopard_nizk2_preproc_server_data(&prep,&pre_st,&pp,&sk,&preq);
        clock_gettime(CLOCK_MONOTONIC,&p1);
        if(rc != 0 ||
           leopard_write_prep_file(prep_path,&prep) != 0 ||
           leopard_write_server_pre_state_file(sps_path,&pre_st) != 0) {
          rc = 1;
          break;
        }
        total_party_sec += elapsed_sec_ts(&p0,&p1);
        prep_total_bytes += file_size_bytes(prep_path);
        leopard_sk_data_clear(&sk);
        leopard_prep_data_clear(&prep);
        leopard_server_pre_state_data_clear(&pre_st);
      }
      leopard_mask_list_clear(&local);
      if(rc != 0)
        break;
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0) {
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("nservers=%zu\n",nservers);
      printf("threshold=%zu\n",threshold);
      printf("slot_count=%zu\n",preq.slot_count);
      printf("prep_total_bytes=%zu\n",prep_total_bytes);
      printf("prep_avg_bytes=%zu\n",prep_total_bytes / nservers);
      printf("preproc_party_total_sec=%.6f\n",total_party_sec);
      printf("preproc_party_avg_sec=%.6f\n",total_party_sec / (double)nservers);
      rc = 0;
    }
    leopard_sk_data_clear(&sk);
    leopard_prep_data_clear(&prep);
    leopard_server_pre_state_data_clear(&pre_st);
    leopard_preq_data_clear(&preq);
    leopard_mask_list_clear(&local);
    return rc;
  }

  if(strcmp(argv[1],"blindeval") == 0) {
    leopard_pp pp;
    leopard_request_data req = {0};
    leopard_mask_list local = {0};
    leopard_blindeval_job_list jobs = {0};
    uint8_t tag[LEOPARD_TAG_BYTES];
    size_t i, li, ji;
    size_t bundle_total = 0;
    size_t component_total = 0;
    size_t rep_total_bytes = 0;
    size_t rep_public_total_bytes = 0;
    size_t rep_proof_total_bytes = 0;
    double build_total = 0.0;
    double protocol_eval_total = 0.0;
    double nizk_build_total = 0.0;
    double prove_total = 0.0;
    double export_total = 0.0;

    if(argc != 9) {
      usage(argv[0]);
      return 1;
    }
    if(leopard_read_pp_file(argv[2],&pp) != 0 ||
       leopard_parse_hex(tag,sizeof(tag),argv[6]) != 0 ||
       leopard_read_request_file(argv[7],&req) != 0 ||
       leopard_validate_q_mask(req.nservers,req.threshold,req.q_mask) != 0 ||
       ensure_dir(argv[8]) != 0) {
      return 1;
    }
    rc = 0;
    clock_gettime(CLOCK_MONOTONIC,&t0);
    for(i=0;i<req.nservers;i++) {
      if(((req.q_mask >> i) & 1u) == 0)
        continue;
      if(leopard_collect_local_masks(req.nservers,req.threshold,i,&local) != 0) {
        rc = 1;
        break;
      }
      size_t bundle_indices[LEOPARD_BATCH_MAX];
      uint64_t bundle_masks[LEOPARD_BATCH_MAX];
      size_t bundle_count = 0;
      size_t bundle_idx = 0;

      for(li=0;li<local.count;li++) {
        if(leopard_component_owner(local.masks[li],req.q_mask) != (int)i)
          continue;
        bundle_indices[bundle_count] = li;
        bundle_masks[bundle_count] = local.masks[li];
        bundle_count++;
        component_total++;
        if(bundle_count == LEOPARD_BATCH_MAX) {
          if(leopard_blindeval_job_list_push(&jobs,i,bundle_idx,
                                             bundle_indices,bundle_masks,bundle_count) != 0) {
            rc = 1;
            break;
          }
          bundle_idx++;
          bundle_count = 0;
        }
      }
      if(rc == 0 && bundle_count != 0) {
        if(leopard_blindeval_job_list_push(&jobs,i,bundle_idx,
                                           bundle_indices,bundle_masks,bundle_count) != 0) {
          rc = 1;
        } else {
          bundle_idx++;
        }
      }
      leopard_mask_list_clear(&local);
      if(rc != 0)
        break;
    }

    if(rc == 0) {
#if defined(_OPENMP) && LEOPARD_PARALLEL_BLINDEVAL
#pragma omp parallel for schedule(dynamic) if(jobs.count > 1)
#endif
      for(ji=0;ji<jobs.count;ji++) {
        leopard_blindeval_job *job = &jobs.jobs[ji];

        job->rc = leopard_run_blindeval_bundle(&pp,argv[3],argv[4],argv[5],argv[8],
                                               tag,sizeof(tag),&req,job->server_idx,
                                               job->bundle_idx,job->local_indices,
                                               job->masks,job->count,
                                               &job->rep_total_bytes,
                                               &job->rep_public_total_bytes,
                                               &job->rep_proof_total_bytes,
                                               &job->build_total,
                                               &job->protocol_eval_total,
                                               &job->nizk_build_total,
                                               &job->prove_total,
                                               &job->export_total);
      }

      for(ji=0;ji<jobs.count;ji++) {
        leopard_blindeval_job *job = &jobs.jobs[ji];

        if(job->rc != 0) {
          rc = 1;
          break;
        }
        rep_total_bytes += job->rep_total_bytes;
        rep_public_total_bytes += job->rep_public_total_bytes;
        rep_proof_total_bytes += job->rep_proof_total_bytes;
        build_total += job->build_total;
        protocol_eval_total += job->protocol_eval_total;
        nizk_build_total += job->nizk_build_total;
        prove_total += job->prove_total;
        export_total += job->export_total;
      }
    }
    bundle_total = jobs.count;
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    if(rc == 0) {
      printf("ok=1\n");
      printf("time_sec=%.6f\n",sec);
      printf("nservers=%zu\n",req.nservers);
      printf("threshold=%zu\n",req.threshold);
      printf("q_mask=0x%llx\n",(unsigned long long)req.q_mask);
      printf("bundle_count=%zu\n",bundle_total);
      printf("component_count=%zu\n",component_total);
      printf("rep_total_bytes=%zu\n",rep_total_bytes);
      printf("rep_public_total_bytes=%zu\n",rep_public_total_bytes);
      printf("rep_proof_total_bytes=%zu\n",rep_proof_total_bytes);
      printf("rep_avg_bytes=%.6f\n",(double)rep_total_bytes / req.threshold);
      printf("rep_public_avg_bytes=%.6f\n",(double)rep_public_total_bytes / req.threshold);
      printf("rep_proof_avg_bytes=%.6f\n",(double)rep_proof_total_bytes / req.threshold);
      printf("blindeval_protocol_eval_total_sec=%.6f\n",protocol_eval_total);
      printf("blindeval_nizk_build_total_sec=%.6f\n",nizk_build_total);
      printf("blindeval_build_instance_total_sec=%.6f\n",build_total);
      printf("blindeval_prover_total_sec=%.6f\n",prove_total);
      printf("blindeval_export_total_sec=%.6f\n",export_total);
      printf("blindeval_protocol_eval_avg_sec=%.6f\n",protocol_eval_total / (double)req.threshold);
      printf("blindeval_nizk_build_avg_sec=%.6f\n",nizk_build_total / (double)req.threshold);
      printf("blindeval_build_instance_avg_sec=%.6f\n",build_total / (double)req.threshold);
      printf("blindeval_prover_avg_sec=%.6f\n",prove_total / (double)req.threshold);
      printf("blindeval_export_avg_sec=%.6f\n",export_total / (double)req.threshold);
      rc = 0;
    }
    leopard_request_data_clear(&req);
    leopard_mask_list_clear(&local);
    leopard_blindeval_job_list_clear(&jobs);
    return rc;
  }

  if(strcmp(argv[1],"verify-replies") == 0) {
    leopard_pp pp;
    leopard_request_data req = {0};
    leopard_mask_list local = {0};
    leopard_blindeval_job_list jobs = {0};
    uint8_t tag[LEOPARD_TAG_BYTES];
    size_t i, li, ji;
    size_t accepted = 0;
    size_t component_total = 0;
    double verify_import_total = 0.0;
    double verify_build_total = 0.0;
    double verify_core_total = 0.0;
    double verify_profile_total = 0.0;

    if(argc != 8) {
      usage(argv[0]);
      return 1;
    }
    if(leopard_read_pp_file(argv[2],&pp) != 0 ||
       leopard_parse_hex(tag,sizeof(tag),argv[5]) != 0 ||
       leopard_read_request_file(argv[6],&req) != 0 ||
       leopard_validate_q_mask(req.nservers,req.threshold,req.q_mask) != 0) {
      return 1;
    }
    clock_gettime(CLOCK_MONOTONIC,&t0);
    rc = 0;
    for(i=0;i<req.nservers;i++) {
      if(((req.q_mask >> i) & 1u) == 0)
        continue;
      if(leopard_collect_local_masks(req.nservers,req.threshold,i,&local) != 0) {
        rc = 1;
        break;
      }
      size_t bundle_indices[LEOPARD_BATCH_MAX];
      uint64_t bundle_masks[LEOPARD_BATCH_MAX];
      size_t bundle_count = 0;
      size_t bundle_idx = 0;

      for(li=0;li<local.count;li++) {
        if(leopard_component_owner(local.masks[li],req.q_mask) != (int)i)
          continue;
        bundle_indices[bundle_count] = li;
        bundle_masks[bundle_count] = local.masks[li];
        bundle_count++;
        component_total++;
        if(bundle_count == LEOPARD_BATCH_MAX) {
          if(leopard_blindeval_job_list_push(&jobs,i,bundle_idx,
                                             bundle_indices,bundle_masks,bundle_count) != 0) {
            rc = 1;
            break;
          }
          bundle_idx++;
          bundle_count = 0;
        }
      }
      if(rc == 0 && bundle_count != 0) {
        if(leopard_blindeval_job_list_push(&jobs,i,bundle_idx,
                                           bundle_indices,bundle_masks,bundle_count) != 0) {
          rc = 1;
        } else {
          bundle_idx++;
        }
      }
      leopard_mask_list_clear(&local);
      if(rc != 0)
        break;
    }

    if(rc == 0) {
#if defined(_OPENMP) && LEOPARD_PARALLEL_VERIFY_REPLIES
#pragma omp parallel for schedule(dynamic) if(jobs.count > 1)
#endif
      for(ji=0;ji<jobs.count;ji++) {
        leopard_blindeval_job *job = &jobs.jobs[ji];

        job->rc = leopard_verify_reply_bundle(&pp,argv[3],argv[4],argv[7],
                                              tag,sizeof(tag),&req,job->server_idx,
                                              job->bundle_idx,job->local_indices,
                                              job->masks,job->count,
                                              &job->verify_import_total,
                                              &job->verify_build_total,
                                              &job->verify_core_total,
                                              &job->verify_total);
      }
      for(ji=0;ji<jobs.count;ji++) {
        if(jobs.jobs[ji].rc != 0) {
          rc = 1;
          break;
        }
        verify_import_total += jobs.jobs[ji].verify_import_total;
        verify_build_total += jobs.jobs[ji].verify_build_total;
        verify_core_total += jobs.jobs[ji].verify_core_total;
        verify_profile_total += jobs.jobs[ji].verify_total;
        accepted++;
      }
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    sec = elapsed_sec_ts(&t0,&t1);
    printf("ok=%d\n",rc == 0 ? 1 : 0);
    printf("time_sec=%.6f\n",sec);
    printf("accepted=%zu\n",accepted);
    printf("component_count=%zu\n",component_total);
    printf("nservers=%zu\n",req.nservers);
    printf("threshold=%zu\n",req.threshold);
    printf("q_mask=0x%llx\n",(unsigned long long)req.q_mask);
    {
      const double denom = accepted == 0 ? 1.0 : (double)accepted;
      printf("parallel_verify_replies=%d\n",LEOPARD_PARALLEL_VERIFY_REPLIES);
      printf("verify_reply_import_total_sec=%.6f\n",verify_import_total);
      printf("verify_reply_build_instance_total_sec=%.6f\n",verify_build_total);
      printf("verify_reply_verifier_total_sec=%.6f\n",verify_core_total);
      printf("verify_reply_profile_total_sec=%.6f\n",verify_profile_total);
      printf("verify_reply_import_avg_sec=%.6f\n",verify_import_total / denom);
      printf("verify_reply_build_instance_avg_sec=%.6f\n",verify_build_total / denom);
      printf("verify_reply_verifier_avg_sec=%.6f\n",verify_core_total / denom);
      printf("verify_reply_profile_avg_sec=%.6f\n",verify_profile_total / denom);
    }
    leopard_request_data_clear(&req);
    leopard_mask_list_clear(&local);
    leopard_blindeval_job_list_clear(&jobs);
    return rc;
  }

  usage(argv[0]);
  return 1;
}
