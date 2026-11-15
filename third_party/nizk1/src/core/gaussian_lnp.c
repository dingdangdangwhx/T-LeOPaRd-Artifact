#include "gaussian_lnp.h"

#include "aesctr.h"
#include "malloc.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double INV_LN2 = 1.4426950408889634;
static const double LN2 = 0.6931471805599453;
static const double DSS_155 = 1.0 / (2.0 * 1.55 * 1.55);

typedef struct {
  uint64_t hi;
  uint64_t lo;
} z128;

static const z128 CDF155[] = {
  {10894764499197476522ULL, 10804844707381617341ULL},
  {4761708367981796450ULL, 6209732027000382074ULL},
  {1476784279527800432ULL, 14108379346150813303ULL},
  {316388870594767345ULL, 17827298407763885637ULL},
  {46043503515468600ULL, 18385899657892021654ULL},
  {4503729779335039ULL, 3860889375818664979ULL},
  {294122444862326ULL, 13947176349836216550ULL},
  {12769598070895ULL, 14321894682751135119ULL},
  {367552986472ULL, 10286761328368440884ULL},
  {7001273393ULL, 2287787188970898528ULL},
  {88153536ULL, 17843977990435837663ULL},
  {733119ULL, 12174894787802692461ULL},
  {4024ULL, 18067426722645776197ULL},
  {14ULL, 10764017821655913055ULL},
  {0ULL, 643125733022530080ULL},
  {0ULL, 1014291014134832ULL},
  {0ULL, 1055215183460ULL},
  {0ULL, 724109373ULL},
  {0ULL, 327744ULL},
  {0ULL, 98ULL},
  {0ULL, 0ULL}
};

typedef struct {
  uint8_t buf[AES128CTR_BLOCKBYTES];
  size_t pos;
  uint64_t bits;
  size_t bpos;
} gaussian_randbuf;

static inline void randbuf_init(gaussian_randbuf *rb) {
  rb->pos = AES128CTR_BLOCKBYTES;
  rb->bits = 0;
  rb->bpos = 64;
}

static inline void randbuf_fill(gaussian_randbuf *rb, aes128ctr_ctx *state) {
  aes128ctr_squeezeblocks(rb->buf, 1, state);
  rb->pos = 0;
}

static inline const uint8_t *randbuf_bytes(gaussian_randbuf *rb,
                                           aes128ctr_ctx *state, size_t n) {
  if(rb->pos + n > AES128CTR_BLOCKBYTES)
    randbuf_fill(rb,state);
  const uint8_t *ret = &rb->buf[rb->pos];
  rb->pos += n;
  return ret;
}

static inline int randbuf_bit(gaussian_randbuf *rb, aes128ctr_ctx *state) {
  if(rb->bpos >= 64) {
    memcpy(&rb->bits,randbuf_bytes(rb,state,sizeof(rb->bits)),sizeof(rb->bits));
    rb->bpos = 0;
  }
  const int bit = (int)(rb->bits & 1);
  rb->bits >>= 1;
  rb->bpos++;
  return bit;
}

static inline double exp_small(double x) {
#define C1 (1.66666666666666019037e-01)
#define C2 (-2.77777777770155933842e-03)
#define C3 (6.61375632143793436117e-05)
#define C4 (-1.65339022054652515390e-06)
#define C5 (4.13813679705723846039e-08)
  double t = x*x;
  t = x - t*(C1 + t*(C2 + t*(C3 + t*(C4 + t*C5))));
  t = 1.0 - ((x*t)/(t - 2.0) - x);
  return t;
#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
}

static inline unsigned int cdf_sample(aes128ctr_ctx *state,
                                      gaussian_randbuf *rb) {
  uint64_t w[2];
  unsigned int z = 0;
  memcpy(w,randbuf_bytes(rb,state,sizeof(w)),sizeof(w));
  while(w[0] <= CDF155[z].hi &&
        (w[0] < CDF155[z].hi || w[1] < CDF155[z].lo))
    z++;
  return z;
}

static inline unsigned int ber_exp(aes128ctr_ctx *state,
                                   gaussian_randbuf *rb, double x) {
  uint64_t w[2];
  uint64_t t;
  uint64_t u;
  unsigned int bit;

  memcpy(w,randbuf_bytes(rb,state,sizeof(w)),sizeof(w));
  t = w[0];
  u = (uint64_t)(x*INV_LN2);
  x -= LN2*(double)u;
  u ^= (u ^ 63) & (uint64_t)((int64_t)(63 - u) >> 63);
  t ^= (t >> u) << u;
  bit = (unsigned int)(1 - ((t | -t) >> 63));

  t = w[1] & ((1ULL << 53) - 1);
  u = (uint64_t)(exp_small(-x)*(double)(1ULL << 53));
  bit &= (unsigned int)((t - u) >> 63);
  return bit;
}

static inline int gaussian155(aes128ctr_ctx *state, gaussian_randbuf *rb,
                              double center) {
  int k;
  int bit;
  double x;
  do {
    bit = randbuf_bit(rb,state);
    k = (int)cdf_sample(state,rb);
    k = (-bit & (2*k)) - k + bit;
    x = ((k - center)*(k - center) - (k - bit)*(k - bit))*DSS_155;
  } while(!ber_exp(state,rb,x));
  return k;
}

void gaussian_lnp_i64(int64_t *out, size_t count, const uint8_t seed[16],
                      uint64_t nonce, unsigned int log2sd) {
  const size_t nbits = (size_t)log2sd*count;
  const size_t nblocks = (nbits + AES128CTR_BLOCKBYTES*8 - 1)/
                         (AES128CTR_BLOCKBYTES*8) + 1;
  const uint32_t mask = ((uint32_t)1 << log2sd) - 1;
  const double inv2pow = 1.0/(double)((uint64_t)1 << log2sd);
  aes128ctr_ctx state;
  gaussian_randbuf rb;
  uint8_t *uniform;
  size_t i;

  if(log2sd == 0 || log2sd >= 32) {
    fprintf(stderr,"ERROR: unsupported Gaussian log2sd=%u\n",log2sd);
    exit(1);
  }

  uniform = _malloc(nblocks*AES128CTR_BLOCKBYTES);
  aes128ctr_init(&state,seed,nonce);
  aes128ctr_squeezeblocks(uniform,nblocks,&state);
  randbuf_init(&rb);

  for(i=0;i<count;i++) {
    const size_t bitoff = i*log2sd;
    uint64_t word;
    uint32_t u;
    int k;
    memcpy(&word,&uniform[bitoff >> 3],sizeof(word));
    u = (uint32_t)((word >> (bitoff & 7)) & mask);
    k = gaussian155(&state,&rb,(double)u*inv2pow);
    out[i] = (int64_t)k*((int64_t)1 << log2sd) - (int64_t)u;
  }

  free(uniform);
}
