#ifndef LEOPARD_NIZK1_H
#define LEOPARD_NIZK1_H

#include <stddef.h>
#include <stdint.h>

#include "common.h"

typedef struct leopard_nizk1_request leopard_nizk1_request;
typedef struct leopard_nizk1_state leopard_nizk1_state;

typedef struct {
  double total_sec;
  double selector_material_sec;
  double slot_load_sec;
  double chain_sec;
  double cx_sec;
} leopard_nizk1_fast_public_profile;

typedef struct {
  double protocol_input_sec;
  double protocol_transcript_sec;
  double protocol_selector_material_sec;
  double protocol_slot_load_sec;
  double protocol_hash_eval_sec;
  double protocol_cx_sec;
  double protocol_state_sec;
  double protocol_public_export_sec;
  double nizk_prepare_sec;
  double nizk_inner_setup_sec;
  double inner_instance_sec;
  double inner_prepare_sec;
  double inner_prove_sec;
  double inner_selfcheck_sec;
  size_t inner_rejection_rounds;
  size_t inner_batch_rounds;
  int a1open_cache;
  double outer_witness_sec;
  double outer_statement_sec;
  double outer_stmt_fold_sec;
  double outer_stmt_quad_sec;
  double outer_stmt_cnst1_sec;
  double outer_stmt_cnst2_sec;
  double outer_stmt_cnst3_sec;
  double outer_stmt_hash_sec;
  double outer_prove_sec;
  double serialize_sec;
  double export_sec;
} leopard_nizk1_request_profile;

typedef struct {
  double public_ctx_sec;
  int a1open_cache;
  double fs_seed_sec;
  double deserialize_sec;
  double outer_statement_sec;
  double outer_stmt_fold_sec;
  double outer_stmt_quad_sec;
  double outer_stmt_cnst1_sec;
  double outer_stmt_cnst2_sec;
  double outer_stmt_cnst3_sec;
  double outer_stmt_hash_sec;
  double outer_verify_sec;
} leopard_nizk1_verify_profile;

int leopard_nizk1_request_prove_data(leopard_request_data *out_req,
                                     leopard_state_data *out_st,
                                     const leopard_pp *pp,
                                     const uint8_t *tag, size_t tag_len,
                                     const uint8_t xbytes[LEOPARD_X_BYTES]);

int leopard_nizk1_preproc_client_data(leopard_preq_data *out_preq,
                                      leopard_pre_state_data *out_st,
                                      const leopard_pp *pp,
                                      size_t slot_count);

int leopard_nizk1_request_prove_preprocessed_data(leopard_request_data *out_req,
                                                  leopard_state_data *out_st,
                                                  leopard_pre_state_data *io_pre_st,
                                                  const leopard_pp *pp,
                                                  const uint8_t *tag, size_t tag_len,
                                                  const uint8_t xbytes[LEOPARD_X_BYTES]);

int leopard_nizk1_request_fast_public_preprocessed_data(
    leopard_nizk1_fast_public_profile *out_prof,
    const leopard_pp *pp,
    const leopard_pre_state_data *pre,
    size_t slot_j,
    const uint8_t *tag, size_t tag_len,
    const uint8_t xbytes[LEOPARD_X_BYTES],
    const leopard_request_data *expected_req);

int leopard_nizk1_request_verify_data(const leopard_pp *pp,
                                      const uint8_t *tag, size_t tag_len,
                                      const leopard_request_data *req);

size_t leopard_nizk1_request_data_bytes(const leopard_request_data *req);

int leopard_nizk1_finalize_output_data(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                       const leopard_state_data *st,
                                       const leopard_reply_data *rep);

int leopard_nizk1_finalize_raw_output_data(uint8_t y[LEOPARD_OUTPUT_BYTES],
                                           const leopard_state_data *st,
                                           const leopard_coeff *ux_coeffs,
                                           size_t ux_rows,
                                           const leopard_coeff *vk_coeffs,
                                           size_t vk_rows);

int leopard_nizk1_eval(uint8_t y[LEOPARD_OUTPUT_BYTES],
                       const leopard_pp *pp,
                       const uint8_t *tag, size_t tag_len,
                       const uint8_t xbytes[LEOPARD_X_BYTES],
                       const leopard_coeff *k_coeffs, size_t k_rows);

const leopard_nizk1_request_profile *leopard_nizk1_last_request_profile(void);
const leopard_nizk1_verify_profile *leopard_nizk1_last_verify_profile(void);

#endif
