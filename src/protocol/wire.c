#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wire.h"

enum {
  LEOPARD_MAGIC_PP = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x50505231U)),
  LEOPARD_MAGIC_PREQ = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x50525131U)),
  LEOPARD_MAGIC_PREP = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x50525031U)),
  LEOPARD_MAGIC_PST = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x50535431U)),
  LEOPARD_MAGIC_SPS = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x53505331U)),
  LEOPARD_MAGIC_REQ = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x52455131U)),
  LEOPARD_MAGIC_REQ_CACHED = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x52455132U)),
  LEOPARD_MAGIC_ST = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x53544131U)),
  LEOPARD_MAGIC_PK = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (NIZK2_DOUBLE_COMMITMENT ? 0x504b5933U : (NIZK2_SMALL_COMMITMENT ? 0x504b5932U : 0x504b5931U))),
  LEOPARD_MAGIC_SK = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (NIZK2_DOUBLE_COMMITMENT ? 0x534b5932U : 0x534b5931U)),
  LEOPARD_MAGIC_REP = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x52455031U)),
  LEOPARD_MAGIC_MSK = ((LEOPARD_LOGQ > 64 ? 0x00800000U : 0) ^ (0x4d534b31U))
};

static leopard_coeff *alloc_i64(size_t count) {
  leopard_coeff *buf = malloc(count * sizeof(*buf));

  if(buf == NULL)
    return NULL;
  memset(buf,0,count * sizeof(*buf));
  return buf;
}

static uint8_t *alloc_u8(size_t count) {
  uint8_t *buf = malloc(count);

  if(buf == NULL)
    return NULL;
  memset(buf,0,count);
  return buf;
}

static int write_bytes(FILE *fp, const void *buf, size_t len) {
  return fwrite(buf,1,len,fp) == len ? 0 : 1;
}

static int read_bytes(FILE *fp, void *buf, size_t len) {
  return fread(buf,1,len,fp) == len ? 0 : 1;
}

static int write_u64(FILE *fp, uint64_t v) {
  return write_bytes(fp,&v,sizeof(v));
}

static int read_u64(FILE *fp, uint64_t *v) {
  return read_bytes(fp,v,sizeof(*v));
}

static int write_i64_array(FILE *fp, const leopard_coeff *v, size_t count) {
#if LEOPARD_LOGQ > 64
  for(size_t i=0;i<count;i++) {
    uint8_t bytes[LEOPARD_COEFF_WIRE_BYTES];
    if(v[i] < -((leopard_coeff)LEOPARD_MOD_Q/2) || v[i] > (leopard_coeff)LEOPARD_MOD_Q/2) return 1;
    for(size_t j=0;j<sizeof(bytes);j++) bytes[j]=(uint8_t)((leopard_ucoeff)v[i]>>(8*j));
    if(write_bytes(fp,bytes,sizeof(bytes))) return 1;
  }
  return 0;
#else
  return write_bytes(fp,v,count * sizeof(*v));
#endif
}

static int read_i64_array(FILE *fp, leopard_coeff *v, size_t count) {
#if LEOPARD_LOGQ > 64
  for(size_t i=0;i<count;i++) {
    uint8_t bytes[LEOPARD_COEFF_WIRE_BYTES];
    leopard_ucoeff x=0;
    if(read_bytes(fp,bytes,sizeof(bytes))) return 1;
    for(size_t j=0;j<sizeof(bytes);j++) x|=(leopard_ucoeff)bytes[j]<<(8*j);
    v[i]=(bytes[sizeof(bytes)-1]&128) ? (leopard_coeff)x-((leopard_coeff)1<<72) : (leopard_coeff)x;
    if(v[i] < -((leopard_coeff)LEOPARD_MOD_Q/2) || v[i] > (leopard_coeff)LEOPARD_MOD_Q/2) return 1;
  }
  return 0;
#else
  return read_bytes(fp,v,count * sizeof(*v));
#endif
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

static void pack_i64_offsets(uint8_t *dst, const leopard_coeff *src, size_t count,
                             leopard_coeff minv, unsigned int bits) {
  size_t i;
  size_t bitpos = 0;
  const size_t packed_bytes = packed_payload_bytes(count,bits);

  if(packed_bytes == 0)
    return;
  memset(dst,0,packed_bytes);
  for(i=0;i<count;i++) {
    leopard_ucoeff x = (leopard_ucoeff)src[i] - (leopard_ucoeff)minv;
    unsigned int written = 0;

    while(written < bits) {
      const size_t byte_idx = bitpos >> 3;
      const unsigned int bit_idx = (unsigned int)(bitpos & 7);
      const unsigned int room = 8 - bit_idx;
      const unsigned int take = (bits - written < room) ? (bits - written) : room;
      const uint8_t mask = (take == 8) ? 0xffu : (uint8_t)((1u << take) - 1u);

      dst[byte_idx] |= (uint8_t)(((x >> written) & mask) << bit_idx);
      written += take;
      bitpos += take;
    }
  }
}

static void unpack_i64_offsets(leopard_coeff *dst, size_t count, leopard_coeff minv,
                               unsigned int bits, const uint8_t *src) {
  size_t i;
  size_t bitpos = 0;

  if(bits == 0) {
    for(i=0;i<count;i++)
      dst[i] = minv;
    return;
  }
  for(i=0;i<count;i++) {
    leopard_ucoeff x = 0;
    unsigned int read = 0;

    while(read < bits) {
      const size_t byte_idx = bitpos >> 3;
      const unsigned int bit_idx = (unsigned int)(bitpos & 7);
      const unsigned int room = 8 - bit_idx;
      const unsigned int take = (bits - read < room) ? (bits - read) : room;
      const uint8_t mask = (take == 8) ? 0xffu : (uint8_t)((1u << take) - 1u);

      x |= (leopard_ucoeff)((src[byte_idx] >> bit_idx) & mask) << read;
      read += take;
      bitpos += take;
    }
    dst[i] = minv + (leopard_coeff)x;
  }
}

static int write_packed_i64_array(FILE *fp, const leopard_coeff *v, size_t count) {
  leopard_coeff minv, maxv;
  uint8_t bits8;
  uint8_t *buf = NULL;
  size_t i;
  size_t packed_bytes;

  if(count == 0)
    return 1;
  minv = v[0];
  maxv = v[0];
  for(i=1;i<count;i++) {
    if(v[i] < minv)
      minv = v[i];
    if(v[i] > maxv)
      maxv = v[i];
  }
  bits8 = (uint8_t)packed_bits_u64((leopard_ucoeff)maxv - (leopard_ucoeff)minv);
  if(bits8 > LEOPARD_LOGQ || count > (SIZE_MAX-7)/(bits8 ? bits8 : 1)) return 1;
#if LEOPARD_LOGQ > 64
  if(maxv > (leopard_coeff)LEOPARD_MOD_Q/2) return 1;
#endif
  packed_bytes = packed_payload_bytes(count,bits8);
  if(write_i64_array(fp,&minv,1) != 0 ||
     write_bytes(fp,&bits8,sizeof(bits8)) != 0)
    return 1;
  if(packed_bytes == 0)
    return 0;
  buf = alloc_u8(packed_bytes);
  if(buf == NULL)
    return 1;
  pack_i64_offsets(buf,v,count,minv,bits8);
  i = write_bytes(fp,buf,packed_bytes);
  free(buf);
  return i;
}

static int read_packed_i64_array(FILE *fp, leopard_coeff *v, size_t count) {
  leopard_coeff minv;
  uint8_t bits8;
  uint8_t *buf = NULL;
  size_t packed_bytes;
  int ret = 1;

  if(count == 0)
    return 1;
  if(read_i64_array(fp,&minv,1) != 0 ||
     read_bytes(fp,&bits8,sizeof(bits8)) != 0)
    return 1;
  if(bits8 > LEOPARD_LOGQ || count > (SIZE_MAX-7)/(bits8 ? bits8 : 1)) return 1;
  packed_bytes = packed_payload_bytes(count,bits8);
  if(packed_bytes == 0) {
    unpack_i64_offsets(v,count,minv,bits8,NULL);
    return 0;
  }
  buf = alloc_u8(packed_bytes);
  if(buf == NULL)
    return 1;
  if(read_bytes(fp,buf,packed_bytes) != 0)
    goto end;
  unpack_i64_offsets(v,count,minv,bits8,buf);
#if LEOPARD_LOGQ > 64
  for(size_t k=0;k<count;k++)
    if(v[k] < -((leopard_coeff)LEOPARD_MOD_Q/2) || v[k] > (leopard_coeff)LEOPARD_MOD_Q/2) goto end;
#endif
  ret = 0;

end:
  free(buf);
  return ret;
}

void leopard_request_data_clear(leopard_request_data *req) {
  free(req->cr_coeffs);
  free(req->cx_coeffs);
  free(req->tA_coeffs);
  free(req->w_coeffs);
  free(req->rbit_const_tg_coeffs);
  free(req->rbit_const_h_coeffs);
  free(req->vlin_coeffs);
  free(req->tq_coeffs);
  free(req->vq_coeffs);
  free(req->proof);
  memset(req,0,sizeof(*req));
}

void leopard_state_data_clear(leopard_state_data *st) {
  free(st->rhat_coeffs);
  memset(st,0,sizeof(*st));
}

void leopard_pk_data_clear(leopard_pk_data *pk) {
  free(pk->ck_coeffs);
  memset(pk,0,sizeof(*pk));
}

void leopard_sk_data_clear(leopard_sk_data *sk) {
  free(sk->r_coeffs);
  free(sk->e_coeffs);
  free(sk->ei_coeffs);
  free(sk->k_coeffs);
  memset(sk,0,sizeof(*sk));
}

void leopard_reply_data_clear(leopard_reply_data *rep) {
  free(rep->ux_coeffs);
  free(rep->proof);
  memset(rep,0,sizeof(*rep));
}

void leopard_preq_data_clear(leopard_preq_data *preq) {
  free(preq->ar_seeds);
  free(preq->cr_coeffs);
  memset(preq,0,sizeof(*preq));
}

void leopard_prep_data_clear(leopard_prep_data *prep) {
  free(prep->vk_coeffs);
  memset(prep,0,sizeof(*prep));
}

void leopard_pre_state_data_clear(leopard_pre_state_data *st) {
  free(st->ar_seeds);
  free(st->cr_coeffs);
  free(st->rhat_coeffs);
  free(st->eua_coeffs);
  free(st->eub_coeffs);
  free(st->r_coeffs);
  free(st->va_coeffs);
  free(st->vb_coeffs);
  memset(st,0,sizeof(*st));
}

void leopard_server_pre_state_data_clear(leopard_server_pre_state_data *st) {
  free(st->vk_coeffs);
  free(st->es_coeffs);
  memset(st,0,sizeof(*st));
}

static int open_write(FILE **out, const char *path) {
  *out = fopen(path,"wb");
  return *out == NULL ? 1 : 0;
}

static int open_read(FILE **out, const char *path) {
  *out = fopen(path,"rb");
  return *out == NULL ? 1 : 0;
}

int leopard_write_pp_file(const char *path, const leopard_pp *pp) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_PP;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_bytes(fp,pp,sizeof(*pp)) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_pp_file(const char *path, leopard_pp *pp) {
  FILE *fp;
  uint32_t magic;

  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_PP ||
     read_bytes(fp,pp,sizeof(*pp)) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_write_preq_file(const char *path, const leopard_preq_data *preq) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_PREQ;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,preq->slot_count) != 0 ||
     write_u64(fp,preq->cr_rows) != 0 ||
     write_bytes(fp,preq->ar_seeds,16 * preq->slot_count) != 0 ||
     write_packed_i64_array(fp,preq->cr_coeffs,
                            preq->slot_count * preq->cr_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_preq_file(const char *path, leopard_preq_data *preq) {
  FILE *fp;
  uint32_t magic;
  uint64_t tmp64;
  size_t coeff_count;

  memset(preq,0,sizeof(*preq));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_PREQ ||
     read_u64(fp,&tmp64) != 0) {
    fclose(fp);
    return 1;
  }
  preq->slot_count = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0)
    goto fail;
  preq->cr_rows = (size_t)tmp64;
  coeff_count = preq->slot_count * preq->cr_rows * LEOPARD_RING_DEG;
  preq->ar_seeds = alloc_u8(16 * preq->slot_count);
  preq->cr_coeffs = alloc_i64(coeff_count);
  if(preq->ar_seeds == NULL || preq->cr_coeffs == NULL)
    goto fail;
  if(read_bytes(fp,preq->ar_seeds,16 * preq->slot_count) != 0 ||
     read_packed_i64_array(fp,preq->cr_coeffs,coeff_count) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_preq_data_clear(preq);
  return 1;
}

int leopard_write_prep_file(const char *path, const leopard_prep_data *prep) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_PREP;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,prep->slot_count) != 0 ||
     write_u64(fp,prep->component_mask) != 0 ||
     write_u64(fp,prep->vk_rows) != 0 ||
     write_packed_i64_array(fp,prep->vk_coeffs,
                            prep->slot_count * prep->vk_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_prep_file(const char *path, leopard_prep_data *prep) {
  FILE *fp;
  uint32_t magic;
  uint64_t tmp64;
  size_t coeff_count;

  memset(prep,0,sizeof(*prep));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_PREP ||
     read_u64(fp,&tmp64) != 0) {
    fclose(fp);
    return 1;
  }
  prep->slot_count = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0)
    goto fail;
  prep->component_mask = tmp64;
  if(read_u64(fp,&tmp64) != 0)
    goto fail;
  prep->vk_rows = (size_t)tmp64;
  coeff_count = prep->slot_count * prep->vk_rows * LEOPARD_RING_DEG;
  prep->vk_coeffs = alloc_i64(coeff_count);
  if(prep->vk_coeffs == NULL)
    goto fail;
  if(read_packed_i64_array(fp,prep->vk_coeffs,coeff_count) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_prep_data_clear(prep);
  return 1;
}

int leopard_write_pre_state_file(const char *path, const leopard_pre_state_data *st) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_PST;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,st->slot_count) != 0 ||
     write_u64(fp,st->next_j) != 0 ||
     write_u64(fp,st->cr_rows) != 0 ||
     write_u64(fp,st->rhat_rows) != 0 ||
     write_u64(fp,st->eua_rows) != 0 ||
     write_u64(fp,st->eub_rows) != 0 ||
     write_u64(fp,st->r_rows) != 0 ||
     write_u64(fp,st->va_rows) != 0 ||
     write_u64(fp,st->vb_rows) != 0 ||
     write_bytes(fp,st->ar_seeds,16 * st->slot_count) != 0 ||
     write_packed_i64_array(fp,st->cr_coeffs,st->slot_count * st->cr_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->rhat_coeffs,st->slot_count * st->rhat_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->eua_coeffs,st->slot_count * st->eua_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->eub_coeffs,st->slot_count * st->eub_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->r_coeffs,st->slot_count * st->r_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->va_coeffs,st->slot_count * st->va_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->vb_coeffs,st->slot_count * st->vb_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_pre_state_file(const char *path, leopard_pre_state_data *st) {
  FILE *fp;
  uint32_t magic;
  uint64_t tmp64;

  memset(st,0,sizeof(*st));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_PST ||
     read_u64(fp,&tmp64) != 0) {
    fclose(fp);
    return 1;
  }
  st->slot_count = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->next_j = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->cr_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->rhat_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->eua_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->eub_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->r_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->va_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->vb_rows = (size_t)tmp64;
  st->ar_seeds = alloc_u8(16 * st->slot_count);
  st->cr_coeffs = alloc_i64(st->slot_count * st->cr_rows * LEOPARD_RING_DEG);
  st->rhat_coeffs = alloc_i64(st->slot_count * st->rhat_rows * LEOPARD_RING_DEG);
  st->eua_coeffs = alloc_i64(st->slot_count * st->eua_rows * LEOPARD_RING_DEG);
  st->eub_coeffs = alloc_i64(st->slot_count * st->eub_rows * LEOPARD_RING_DEG);
  st->r_coeffs = alloc_i64(st->slot_count * st->r_rows * LEOPARD_RING_DEG);
  st->va_coeffs = alloc_i64(st->slot_count * st->va_rows * LEOPARD_RING_DEG);
  st->vb_coeffs = alloc_i64(st->slot_count * st->vb_rows * LEOPARD_RING_DEG);
  if(st->ar_seeds == NULL || st->cr_coeffs == NULL || st->rhat_coeffs == NULL ||
     st->eua_coeffs == NULL || st->eub_coeffs == NULL || st->r_coeffs == NULL ||
     st->va_coeffs == NULL || st->vb_coeffs == NULL)
    goto fail;
  if(read_bytes(fp,st->ar_seeds,16 * st->slot_count) != 0 ||
     read_packed_i64_array(fp,st->cr_coeffs,st->slot_count * st->cr_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->rhat_coeffs,st->slot_count * st->rhat_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->eua_coeffs,st->slot_count * st->eua_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->eub_coeffs,st->slot_count * st->eub_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->r_coeffs,st->slot_count * st->r_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->va_coeffs,st->slot_count * st->va_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->vb_coeffs,st->slot_count * st->vb_rows * LEOPARD_RING_DEG) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_pre_state_data_clear(st);
  return 1;
}

int leopard_write_server_pre_state_file(const char *path, const leopard_server_pre_state_data *st) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_SPS;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,st->slot_count) != 0 ||
     write_u64(fp,st->component_mask) != 0 ||
     write_u64(fp,st->vk_rows) != 0 ||
     write_u64(fp,st->es_rows) != 0 ||
     write_packed_i64_array(fp,st->vk_coeffs,st->slot_count * st->vk_rows * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,st->es_coeffs,st->slot_count * st->es_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_server_pre_state_file(const char *path, leopard_server_pre_state_data *st) {
  FILE *fp;
  uint32_t magic;
  uint64_t tmp64;

  memset(st,0,sizeof(*st));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_SPS ||
     read_u64(fp,&tmp64) != 0) {
    fclose(fp);
    return 1;
  }
  st->slot_count = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->component_mask = tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->vk_rows = (size_t)tmp64;
  if(read_u64(fp,&tmp64) != 0) goto fail;
  st->es_rows = (size_t)tmp64;
  st->vk_coeffs = alloc_i64(st->slot_count * st->vk_rows * LEOPARD_RING_DEG);
  st->es_coeffs = alloc_i64(st->slot_count * st->es_rows * LEOPARD_RING_DEG);
  if(st->vk_coeffs == NULL || st->es_coeffs == NULL)
    goto fail;
  if(read_packed_i64_array(fp,st->vk_coeffs,st->slot_count * st->vk_rows * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,st->es_coeffs,st->slot_count * st->es_rows * LEOPARD_RING_DEG) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_server_pre_state_data_clear(st);
  return 1;
}

int leopard_write_request_file(const char *path, const leopard_request_data *req) {
  FILE *fp;
  uint32_t magic = req->cr_external ? LEOPARD_MAGIC_REQ_CACHED : LEOPARD_MAGIC_REQ;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,req->nservers) != 0 ||
     write_u64(fp,req->threshold) != 0 ||
     write_u64(fp,req->q_mask) != 0 ||
     write_u64(fp,req->slot_j) != 0 ||
     write_u64(fp,req->cr_rows) != 0 ||
     write_u64(fp,req->cx_rows) != 0 ||
     write_u64(fp,req->cx_cols) != 0 ||
     write_u64(fp,req->kmsis) != 0 ||
     write_u64(fp,req->rbit_const_len) != 0 ||
     write_u64(fp,req->proof_len) != 0 ||
     write_u64(fp,req->outer_betasq) != 0 ||
     write_bytes(fp,req->ar_seed,sizeof(req->ar_seed)) != 0 ||
     (!req->cr_external && write_packed_i64_array(fp,req->cr_coeffs,req->cr_rows * LEOPARD_RING_DEG) != 0) ||
     write_packed_i64_array(fp,req->cx_coeffs,req->cx_rows * req->cx_cols * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->tA_coeffs,req->kmsis * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->w_coeffs,req->kmsis * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->rbit_const_tg_coeffs,req->rbit_const_len * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->rbit_const_h_coeffs,req->rbit_const_len * LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->vlin_coeffs,LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->tq_coeffs,LEOPARD_RING_DEG) != 0 ||
     write_packed_i64_array(fp,req->vq_coeffs,LEOPARD_RING_DEG) != 0 ||
     write_bytes(fp,req->proof,req->proof_len) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_request_file(const char *path, leopard_request_data *req) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;

  memset(req,0,sizeof(*req));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 ||
     (magic != LEOPARD_MAGIC_REQ && magic != LEOPARD_MAGIC_REQ_CACHED) ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  req->nservers = (size_t)rows64;
  req->cr_external = magic == LEOPARD_MAGIC_REQ_CACHED;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->threshold = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->q_mask = rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->slot_j = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->cr_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->cx_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->cx_cols = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->kmsis = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->rbit_const_len = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->proof_len = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  req->outer_betasq = rows64;
  if(read_bytes(fp,req->ar_seed,sizeof(req->ar_seed)) != 0) goto fail;

  if(!req->cr_external)
    req->cr_coeffs = alloc_i64(req->cr_rows * LEOPARD_RING_DEG);
  req->cx_coeffs = alloc_i64(req->cx_rows * req->cx_cols * LEOPARD_RING_DEG);
  req->tA_coeffs = alloc_i64(req->kmsis * LEOPARD_RING_DEG);
  req->w_coeffs = alloc_i64(req->kmsis * LEOPARD_RING_DEG);
  req->rbit_const_tg_coeffs = alloc_i64(req->rbit_const_len * LEOPARD_RING_DEG);
  req->rbit_const_h_coeffs = alloc_i64(req->rbit_const_len * LEOPARD_RING_DEG);
  req->vlin_coeffs = alloc_i64(LEOPARD_RING_DEG);
  req->tq_coeffs = alloc_i64(LEOPARD_RING_DEG);
  req->vq_coeffs = alloc_i64(LEOPARD_RING_DEG);
  req->proof = alloc_u8(req->proof_len);
  if((!req->cr_external && req->cr_coeffs == NULL) || req->cx_coeffs == NULL || req->tA_coeffs == NULL ||
     req->w_coeffs == NULL || req->rbit_const_tg_coeffs == NULL ||
     req->rbit_const_h_coeffs == NULL || req->vlin_coeffs == NULL || req->tq_coeffs == NULL ||
     req->vq_coeffs == NULL || req->proof == NULL)
    goto fail;
  if((!req->cr_external && read_packed_i64_array(fp,req->cr_coeffs,req->cr_rows * LEOPARD_RING_DEG) != 0) ||
     read_packed_i64_array(fp,req->cx_coeffs,req->cx_rows * req->cx_cols * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->tA_coeffs,req->kmsis * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->w_coeffs,req->kmsis * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->rbit_const_tg_coeffs,req->rbit_const_len * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->rbit_const_h_coeffs,req->rbit_const_len * LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->vlin_coeffs,LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->tq_coeffs,LEOPARD_RING_DEG) != 0 ||
     read_packed_i64_array(fp,req->vq_coeffs,LEOPARD_RING_DEG) != 0 ||
     read_bytes(fp,req->proof,req->proof_len) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_request_data_clear(req);
  return 1;
}

int leopard_request_attach_preprocessing(leopard_request_data *req, const leopard_preq_data *preq) {
  if(!preq || !preq->cr_coeffs || !preq->ar_seeds ||
     req->slot_j >= preq->slot_count || req->cr_rows != preq->cr_rows ||
     req->cr_rows > SIZE_MAX/(LEOPARD_RING_DEG*sizeof(leopard_coeff)) ||
     memcmp(req->ar_seed,preq->ar_seeds+16*req->slot_j,16))
    return 1;
  const size_t count=req->cr_rows*LEOPARD_RING_DEG;
  const leopard_coeff *expected=preq->cr_coeffs+req->slot_j*count;
  if(req->cr_coeffs)
    return memcmp(req->cr_coeffs,expected,count*sizeof(*expected)) != 0;
  if(!req->cr_external) return 1;
  req->cr_coeffs=alloc_i64(count);
  if(!req->cr_coeffs) return 1;
  memcpy(req->cr_coeffs,expected,count*sizeof(*expected));
  return 0;
}

int leopard_write_state_file(const char *path, const leopard_state_data *st) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_ST;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,st->nservers) != 0 ||
     write_u64(fp,st->threshold) != 0 ||
     write_u64(fp,st->q_mask) != 0 ||
     write_u64(fp,st->slot_j) != 0 ||
     write_u64(fp,st->rhat_rows) != 0 ||
     write_bytes(fp,st->tag_hash,sizeof(st->tag_hash)) != 0 ||
     write_bytes(fp,st->xbits,sizeof(st->xbits)) != 0 ||
     write_i64_array(fp,st->rhat_coeffs,st->rhat_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_state_file(const char *path, leopard_state_data *st) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;

  memset(st,0,sizeof(*st));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_ST ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  st->nservers = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  st->threshold = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  st->q_mask = rows64;
  if(read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  st->slot_j = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  st->rhat_rows = (size_t)rows64;
  st->rhat_coeffs = alloc_i64(st->rhat_rows * LEOPARD_RING_DEG);
  if(st->rhat_coeffs == NULL ||
     read_bytes(fp,st->tag_hash,sizeof(st->tag_hash)) != 0 ||
     read_bytes(fp,st->xbits,sizeof(st->xbits)) != 0 ||
     read_i64_array(fp,st->rhat_coeffs,st->rhat_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    leopard_state_data_clear(st);
    return 1;
  }
  return fclose(fp);
}

int leopard_write_pk_file(const char *path, const leopard_pk_data *pk) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_PK;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,pk->component_mask) != 0 ||
     write_u64(fp,pk->ck_rows) != 0) {
    fclose(fp);
    return 1;
  }
#if NIZK2_SMALL_COMMITMENT
  const size_t count = pk->ck_rows * LEOPARD_RING_DEG;
  const size_t bytes = packed_payload_bytes(count,NIZK2_COM_BITS);
  uint8_t *buf = alloc_u8(bytes);
  if(buf == NULL) { fclose(fp); return 1; }
  for(size_t i=0;i<count;i++)
    if(pk->ck_coeffs[i] < -NIZK2_COM_Q/2 || pk->ck_coeffs[i] > NIZK2_COM_Q/2) {
      free(buf); fclose(fp); return 1;
    }
  pack_i64_offsets(buf,pk->ck_coeffs,count,-NIZK2_COM_Q/2,NIZK2_COM_BITS);
  int failed = write_u64(fp,NIZK2_COM_Q) || write_bytes(fp,buf,bytes);
  free(buf);
  if(failed) { fclose(fp); return 1; }
#else
  if(write_i64_array(fp,pk->ck_coeffs,pk->ck_rows * LEOPARD_RING_DEG)) {
    fclose(fp); return 1;
  }
#endif
  return fclose(fp);
}

int leopard_read_pk_file(const char *path, leopard_pk_data *pk) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;

  memset(pk,0,sizeof(*pk));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_PK ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  pk->component_mask = rows64;
  if(read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  pk->ck_rows = (size_t)rows64;
#if NIZK2_SMALL_COMMITMENT
  const size_t expected_rows = NIZK2_COM_RANK + (LEOPARD_LOGQ == 56 ? 34 : 37);
  uint64_t qc;
  if(pk->ck_rows != expected_rows || read_u64(fp,&qc) || qc != NIZK2_COM_Q) {
    fclose(fp); return 1;
  }
#endif
  pk->ck_coeffs = alloc_i64(pk->ck_rows * LEOPARD_RING_DEG);
#if NIZK2_SMALL_COMMITMENT
  const size_t count = pk->ck_rows * LEOPARD_RING_DEG;
  const size_t bytes = packed_payload_bytes(count,NIZK2_COM_BITS);
  uint8_t *buf = alloc_u8(bytes);
  int failed = !buf || !pk->ck_coeffs;
  if(!failed) failed = read_bytes(fp,buf,bytes) || fgetc(fp) != EOF;
  if(!failed) {
    unpack_i64_offsets(pk->ck_coeffs,count,-NIZK2_COM_Q/2,NIZK2_COM_BITS,buf);
    for(size_t i=0;i<count;i++)
      if(pk->ck_coeffs[i] > NIZK2_COM_Q/2) failed = 1;
  }
  free(buf);
  if(failed) {
    fclose(fp); leopard_pk_data_clear(pk); return 1;
  }
#else
  if(pk->ck_coeffs == NULL ||
     read_i64_array(fp,pk->ck_coeffs,pk->ck_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    leopard_pk_data_clear(pk);
    return 1;
  }
#endif
  return fclose(fp);
}

int leopard_write_sk_file(const char *path, const leopard_sk_data *sk) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_SK;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,sk->component_mask) != 0 ||
     write_u64(fp,sk->r_rows) != 0 ||
     write_u64(fp,sk->e_rows) != 0 ||
     write_u64(fp,sk->ei_rows) != 0 ||
     write_u64(fp,sk->k_rows) != 0 ||
     write_i64_array(fp,sk->r_coeffs,sk->r_rows * LEOPARD_RING_DEG) != 0 ||
     write_i64_array(fp,sk->e_coeffs,sk->e_rows * LEOPARD_RING_DEG) != 0 ||
     write_i64_array(fp,sk->ei_coeffs,sk->ei_rows * LEOPARD_RING_DEG) != 0 ||
     write_i64_array(fp,sk->k_coeffs,sk->k_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_sk_file(const char *path, leopard_sk_data *sk) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;

  memset(sk,0,sizeof(*sk));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_SK ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  sk->component_mask = rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  sk->r_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  sk->e_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  sk->ei_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  sk->k_rows = (size_t)rows64;
  sk->r_coeffs = alloc_i64(sk->r_rows * LEOPARD_RING_DEG);
  sk->e_coeffs = alloc_i64(sk->e_rows * LEOPARD_RING_DEG);
  sk->ei_coeffs = alloc_i64(sk->ei_rows * LEOPARD_RING_DEG);
  sk->k_coeffs = alloc_i64(sk->k_rows * LEOPARD_RING_DEG);
  if(sk->r_coeffs == NULL || sk->e_coeffs == NULL || sk->ei_coeffs == NULL || sk->k_coeffs == NULL)
    goto fail;
  if(read_i64_array(fp,sk->r_coeffs,sk->r_rows * LEOPARD_RING_DEG) != 0 ||
     read_i64_array(fp,sk->e_coeffs,sk->e_rows * LEOPARD_RING_DEG) != 0 ||
     read_i64_array(fp,sk->ei_coeffs,sk->ei_rows * LEOPARD_RING_DEG) != 0 ||
     read_i64_array(fp,sk->k_coeffs,sk->k_rows * LEOPARD_RING_DEG) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_sk_data_clear(sk);
  return 1;
}

int leopard_write_reply_file(const char *path, const leopard_reply_data *rep) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_REP;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,rep->component_mask) != 0 ||
     write_u64(fp,rep->ux_rows) != 0 ||
     write_u64(fp,rep->proof_len) != 0 ||
     write_i64_array(fp,rep->ux_coeffs,rep->ux_rows * LEOPARD_RING_DEG) != 0 ||
     write_bytes(fp,rep->proof,rep->proof_len) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_reply_file(const char *path, leopard_reply_data *rep) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;

  memset(rep,0,sizeof(*rep));
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_REP ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  rep->component_mask = rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  rep->ux_rows = (size_t)rows64;
  if(read_u64(fp,&rows64) != 0) goto fail;
  rep->proof_len = (size_t)rows64;
  rep->ux_coeffs = alloc_i64(rep->ux_rows * LEOPARD_RING_DEG);
  rep->proof = alloc_u8(rep->proof_len);
  if(rep->ux_coeffs == NULL || rep->proof == NULL)
    goto fail;
  if(read_i64_array(fp,rep->ux_coeffs,rep->ux_rows * LEOPARD_RING_DEG) != 0 ||
     read_bytes(fp,rep->proof,rep->proof_len) != 0)
    goto fail;
  return fclose(fp);

fail:
  fclose(fp);
  leopard_reply_data_clear(rep);
  return 1;
}

int leopard_write_master_key_file(const char *path,
                                  const leopard_coeff *k_coeffs,
                                  size_t k_rows) {
  FILE *fp;
  uint32_t magic = LEOPARD_MAGIC_MSK;

  if(open_write(&fp,path) != 0)
    return 1;
  if(write_bytes(fp,&magic,sizeof(magic)) != 0 ||
     write_u64(fp,k_rows) != 0 ||
     write_i64_array(fp,k_coeffs,k_rows * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    return 1;
  }
  return fclose(fp);
}

int leopard_read_master_key_file(const char *path,
                                 leopard_coeff **k_coeffs,
                                 size_t *k_rows) {
  FILE *fp;
  uint32_t magic;
  uint64_t rows64;
  leopard_coeff *coeffs;

  *k_coeffs = NULL;
  *k_rows = 0;
  if(open_read(&fp,path) != 0)
    return 1;
  if(read_bytes(fp,&magic,sizeof(magic)) != 0 || magic != LEOPARD_MAGIC_MSK ||
     read_u64(fp,&rows64) != 0) {
    fclose(fp);
    return 1;
  }
  coeffs = alloc_i64((size_t)rows64 * LEOPARD_RING_DEG);
  if(coeffs == NULL ||
     read_i64_array(fp,coeffs,(size_t)rows64 * LEOPARD_RING_DEG) != 0) {
    fclose(fp);
    free(coeffs);
    return 1;
  }
  *k_coeffs = coeffs;
  *k_rows = (size_t)rows64;
  return fclose(fp);
}

static int hex_value(char c) {
  if(c >= '0' && c <= '9')
    return c - '0';
  if(c >= 'a' && c <= 'f')
    return c - 'a' + 10;
  if(c >= 'A' && c <= 'F')
    return c - 'A' + 10;
  return -1;
}

int leopard_parse_hex(uint8_t *out, size_t out_len, const char *hex) {
  size_t i;

  if(strlen(hex) != 2 * out_len)
    return 1;
  for(i=0;i<out_len;i++) {
    int hi = hex_value(hex[2*i]);
    int lo = hex_value(hex[2*i + 1]);

    if(hi < 0 || lo < 0)
      return 1;
    out[i] = (uint8_t)((hi << 4) | lo);
  }
  return 0;
}

void leopard_hex_encode(char *out, const uint8_t *in, size_t in_len) {
  static const char hex[16] = "0123456789abcdef";
  size_t i;

  for(i=0;i<in_len;i++) {
    out[2*i] = hex[in[i] >> 4];
    out[2*i + 1] = hex[in[i] & 0x0f];
  }
  out[2*in_len] = '\0';
}
