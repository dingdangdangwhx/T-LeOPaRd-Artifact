#ifndef LEOPARD_POPRF_COMMON_H
#define LEOPARD_POPRF_COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "../../third_party/nizk2/commitment.h"

#ifndef LEOPARD_LOGQ
#ifdef LOGQ
#define LEOPARD_LOGQ LOGQ
#else
#define LEOPARD_LOGQ 48
#endif
#endif

#if NIZK2_SMALL_COMMITMENT && LEOPARD_LOGQ != 56 && LEOPARD_LOGQ != 64
#error "small server commitment parameters are available only for LOGQ56 and LOGQ64"
#endif

#if LEOPARD_LOGQ == 48
#define LEOPARD_MOD_Q 281474976710597LL
#elif LEOPARD_LOGQ == 56
#define LEOPARD_MOD_Q 72057594037927909LL
#elif LEOPARD_LOGQ == 64
#define LEOPARD_MOD_Q 18446744073709551557ULL
#elif LEOPARD_LOGQ == 72
#define LEOPARD_MOD_Q (((unsigned __int128)1 << 72) - 107)
#else
#error "unsupported LEOPARD_LOGQ; expected 48, 56, 64 or 72"
#endif

#if LEOPARD_LOGQ > 64
typedef __int128 leopard_coeff;
typedef unsigned __int128 leopard_ucoeff;
#define LEOPARD_COEFF_WIRE_BYTES 9
#else
typedef int64_t leopard_coeff;
typedef uint64_t leopard_ucoeff;
#define LEOPARD_COEFF_WIRE_BYTES 8
#endif

enum {
  LEOPARD_TAG_BYTES = 16,
  LEOPARD_OUTPUT_BYTES = 32,
  LEOPARD_X_BITS = 80,
  LEOPARD_X_BYTES = LEOPARD_X_BITS / 8,
  LEOPARD_RING_DEG = 64,
  LEOPARD_Q_BYTES = (LEOPARD_LOGQ + 7) / 8,
  LEOPARD_AR_LIMBS = (LEOPARD_LOGQ + 13) / 14
};

#ifndef NIZK1_SEED_OFFSET
#define NIZK1_SEED_OFFSET 9
#endif

typedef struct {
  uint8_t nizk1_seed[16];
  uint8_t nizk2_pp[32];
} leopard_pp;

typedef struct {
  size_t slot_count;
  size_t cr_rows;
  uint8_t *ar_seeds;
  leopard_coeff *cr_coeffs;
} leopard_preq_data;

typedef struct {
  size_t slot_count;
  uint64_t component_mask;
  size_t vk_rows;
  leopard_coeff *vk_coeffs;
} leopard_prep_data;

typedef struct {
  size_t slot_count;
  size_t next_j;
  size_t cr_rows;
  size_t rhat_rows;
  size_t eua_rows;
  size_t eub_rows;
  size_t r_rows;
  size_t va_rows;
  size_t vb_rows;
  uint8_t *ar_seeds;
  leopard_coeff *cr_coeffs;
  leopard_coeff *rhat_coeffs;
  leopard_coeff *eua_coeffs;
  leopard_coeff *eub_coeffs;
  leopard_coeff *r_coeffs;
  leopard_coeff *va_coeffs;
  leopard_coeff *vb_coeffs;
} leopard_pre_state_data;

typedef struct {
  size_t slot_count;
  uint64_t component_mask;
  size_t vk_rows;
  size_t es_rows;
  leopard_coeff *vk_coeffs;
  leopard_coeff *es_coeffs;
} leopard_server_pre_state_data;

typedef struct {
  size_t nservers;
  size_t threshold;
  uint64_t q_mask;
  size_t slot_j;
  size_t cr_rows;
  size_t cx_rows;
  size_t cx_cols;
  size_t kmsis;
  size_t rbit_const_len;
  size_t proof_len;
  uint64_t outer_betasq;
  uint8_t ar_seed[16];
  leopard_coeff *cr_coeffs;
  leopard_coeff *cx_coeffs;
  leopard_coeff *tA_coeffs;
  leopard_coeff *w_coeffs;
  leopard_coeff *rbit_const_tg_coeffs;
  leopard_coeff *rbit_const_h_coeffs;
  leopard_coeff *vlin_coeffs;
  leopard_coeff *tq_coeffs;
  leopard_coeff *vq_coeffs;
  uint8_t *proof;
} leopard_request_data;

typedef struct {
  size_t nservers;
  size_t threshold;
  uint64_t q_mask;
  size_t slot_j;
  size_t rhat_rows;
  leopard_coeff *rhat_coeffs;
  uint8_t tag_hash[LEOPARD_TAG_BYTES];
  uint8_t xbits[LEOPARD_X_BITS];
} leopard_state_data;

typedef struct {
  uint64_t component_mask;
  size_t ck_rows;
  leopard_coeff *ck_coeffs;
} leopard_pk_data;

typedef struct {
  uint64_t component_mask;
  size_t r_rows;
  size_t e_rows;
  size_t ei_rows;
  size_t k_rows;
  leopard_coeff *r_coeffs;
  leopard_coeff *e_coeffs;
  leopard_coeff *ei_coeffs;
  leopard_coeff *k_coeffs;
} leopard_sk_data;

typedef struct {
  uint64_t component_mask;
  size_t ux_rows;
  size_t proof_len;
  leopard_coeff *ux_coeffs;
  uint8_t *proof;
} leopard_reply_data;

static inline void leopard_pp_init(leopard_pp *pp) {
  static const uint8_t nizk1_base_seed[16] = {
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
    0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f
  };
  static const uint8_t nizk2_pp[32] = {
    0x6c, 0x65, 0x6f, 0x70, 0x61, 0x72, 0x64, 0x2d,
    0x70, 0x70, 0x2d, 0x73, 0x65, 0x65, 0x64, 0x00,
    0x10, 0x32, 0x54, 0x76, 0x98, 0xba, 0xdc, 0xfe,
    0x13, 0x57, 0x9b, 0xdf, 0x24, 0x68, 0xac, 0xe0
  };

  memcpy(pp->nizk1_seed,nizk1_base_seed,sizeof(pp->nizk1_seed));
  {
    uint64_t offset = (uint64_t)NIZK1_SEED_OFFSET;
    int carry = 0;
    size_t i;

    for(i=0;i<8;i++) {
      const uint16_t sum = (uint16_t)pp->nizk1_seed[i] +
                           (uint16_t)(offset & 0xffu) + (uint16_t)carry;
      pp->nizk1_seed[i] = (uint8_t)sum;
      carry = (sum >> 8) & 1;
      offset >>= 8;
    }
    for(i=8;i<16 && carry;i++) {
      const uint16_t sum = (uint16_t)pp->nizk1_seed[i] + (uint16_t)carry;
      pp->nizk1_seed[i] = (uint8_t)sum;
      carry = (sum >> 8) & 1;
    }
  }
  memcpy(pp->nizk2_pp,nizk2_pp,sizeof(pp->nizk2_pp));
}

#endif
