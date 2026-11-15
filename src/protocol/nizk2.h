#ifndef LEOPARD_NIZK2_H
#define LEOPARD_NIZK2_H

#include <stddef.h>
#include <stdint.h>

#include "common.h"

typedef struct {
  double protocol_eval_sec;
  double nizk_build_sec;
  double build_instance_sec;
  double prover_sec;
  double export_sec;
} leopard_nizk2_blindeval_profile;

typedef struct {
  double import_sec;
  double build_instance_sec;
  double verifier_sec;
  double total_sec;
} leopard_nizk2_verify_profile;

int leopard_nizk2_keygen_data(leopard_pk_data *out_pk,
                              leopard_sk_data *out_sk,
                              const leopard_pp *pp);

int leopard_nizk2_keygen_given_key_data(leopard_pk_data *out_pk,
                                        leopard_sk_data *out_sk,
                                        const leopard_pp *pp,
                                        const leopard_coeff *k_coeffs,
                                        uint64_t component_mask,
                                        uint64_t nonce_base);

size_t leopard_nizk2_pk_data_bytes(const leopard_pk_data *pk);

int leopard_nizk2_preproc_server_data(leopard_prep_data *out_prep,
                                      leopard_server_pre_state_data *out_st,
                                      const leopard_pp *pp,
                                      const leopard_sk_data *sk,
                                      const leopard_preq_data *preq);

int leopard_nizk2_blindeval_data(leopard_reply_data *out_rep,
                                 const leopard_pp *pp,
                                 const leopard_pk_data *pk,
                                 const leopard_sk_data *sk,
                                 const uint8_t *tag, size_t tag_len,
                                 const leopard_request_data *req);

int leopard_nizk2_blindeval_preprocessed_data(
    leopard_reply_data *out_rep,
    const leopard_pp *pp,
    const leopard_pk_data *pk,
    const leopard_sk_data *sk,
    const leopard_server_pre_state_data *pre_st,
    const uint8_t *tag, size_t tag_len,
    const leopard_request_data *req);

int leopard_nizk2_reply_verify_data(const leopard_pp *pp,
                                    const leopard_pk_data *pk,
                                    const leopard_prep_data *prep,
                                    const uint8_t *tag, size_t tag_len,
                                    const leopard_request_data *req,
                                    const leopard_reply_data *rep);

size_t leopard_nizk2_reply_data_bytes(const leopard_reply_data *rep);

const leopard_nizk2_blindeval_profile *leopard_nizk2_last_blindeval_profile(void);
const leopard_nizk2_verify_profile *leopard_nizk2_last_verify_profile(void);

#endif
