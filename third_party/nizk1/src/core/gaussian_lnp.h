#ifndef GAUSSIAN_LNP_H
#define GAUSSIAN_LNP_H

#include <stddef.h>
#include <stdint.h>

/* Sample D_{1.55 * 2^log2sd} deterministically from (seed, nonce). */
void gaussian_lnp_i64(int64_t *out, size_t count, const uint8_t seed[16],
                      uint64_t nonce, unsigned int log2sd);

#endif
