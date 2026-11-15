#ifndef GREYHOUND_H
#define GREYHOUND_H

#include <stdint.h>
#include <stddef.h>
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"

#if LOGQ > 64
typedef __int128 polcom_scalar;
#else
typedef int64_t polcom_scalar;
#endif

typedef struct {
  size_t len;
  size_t m;
  size_t n;
  comparams cpp[1];
  const polz *s;
  polx *sx;
  poly *t;
  polz *u1;
  uint8_t h[16];
  uint64_t normsq;  // expected norm
} polcomctx;

typedef struct {
  size_t len;
  size_t m;
  size_t n;
  polcom_scalar x;
  polcom_scalar y;
  comparams cpp[1];
  polz *u1;
  polz *u2;
  uint64_t normsq;
} polcomprf;

void free_polcomctx(polcomctx *ctx);
void free_polcomprf(polcomprf *pi);

void print_polcomctx_pp(const polcomctx *ctx);
double print_polcomprf_pp(const polcomprf *pi);

int polcom_commit(polcomctx *ctx, const polz *s, size_t len);
polcom_scalar polzvec_eval(const polz *a, size_t len, polcom_scalar x);
void polcom_eval(witness *wt, polcomprf *pi, const polcomctx *ctx, polcom_scalar x, polcom_scalar y);
int polcom_reduce(prncplstmnt *st, const polcomprf *pi);

#endif
