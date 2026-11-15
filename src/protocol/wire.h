#ifndef LEOPARD_POPRF_WIRE_H
#define LEOPARD_POPRF_WIRE_H

#include "common.h"

void leopard_request_data_clear(leopard_request_data *req);
void leopard_state_data_clear(leopard_state_data *st);
void leopard_pk_data_clear(leopard_pk_data *pk);
void leopard_sk_data_clear(leopard_sk_data *sk);
void leopard_reply_data_clear(leopard_reply_data *rep);
void leopard_preq_data_clear(leopard_preq_data *preq);
void leopard_prep_data_clear(leopard_prep_data *prep);
void leopard_pre_state_data_clear(leopard_pre_state_data *st);
void leopard_server_pre_state_data_clear(leopard_server_pre_state_data *st);

int leopard_write_pp_file(const char *path, const leopard_pp *pp);
int leopard_read_pp_file(const char *path, leopard_pp *pp);

int leopard_write_preq_file(const char *path, const leopard_preq_data *preq);
int leopard_read_preq_file(const char *path, leopard_preq_data *preq);

int leopard_write_prep_file(const char *path, const leopard_prep_data *prep);
int leopard_read_prep_file(const char *path, leopard_prep_data *prep);

int leopard_write_pre_state_file(const char *path, const leopard_pre_state_data *st);
int leopard_read_pre_state_file(const char *path, leopard_pre_state_data *st);

int leopard_write_server_pre_state_file(const char *path, const leopard_server_pre_state_data *st);
int leopard_read_server_pre_state_file(const char *path, leopard_server_pre_state_data *st);

int leopard_write_request_file(const char *path, const leopard_request_data *req);
int leopard_read_request_file(const char *path, leopard_request_data *req);

int leopard_write_state_file(const char *path, const leopard_state_data *st);
int leopard_read_state_file(const char *path, leopard_state_data *st);

int leopard_write_pk_file(const char *path, const leopard_pk_data *pk);
int leopard_read_pk_file(const char *path, leopard_pk_data *pk);

int leopard_write_sk_file(const char *path, const leopard_sk_data *sk);
int leopard_read_sk_file(const char *path, leopard_sk_data *sk);

int leopard_write_reply_file(const char *path, const leopard_reply_data *rep);
int leopard_read_reply_file(const char *path, leopard_reply_data *rep);

int leopard_write_master_key_file(const char *path,
                                  const leopard_coeff *k_coeffs,
                                  size_t k_rows);
int leopard_read_master_key_file(const char *path,
                                 leopard_coeff **k_coeffs,
                                 size_t *k_rows);

int leopard_parse_hex(uint8_t *out, size_t out_len, const char *hex);
void leopard_hex_encode(char *out, const uint8_t *in, size_t in_len);

#endif
