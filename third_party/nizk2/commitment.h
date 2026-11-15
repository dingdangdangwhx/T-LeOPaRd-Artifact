#ifndef NIZK2_COMMITMENT_H
#define NIZK2_COMMITMENT_H

/* Shared commitment definitions; the artifact uses the full-message mode. */
#define NIZK2_COM_Q 268435399LL
#define NIZK2_COM_RANK 24
#define NIZK2_COM_ALPHABET 8192
#define NIZK2_COM_SCALE (NIZK2_COM_Q / NIZK2_COM_ALPHABET)
#define NIZK2_COM_BITS 28
#define NIZK2_COM_VERSION 1
#ifndef NIZK2_SMALL_COMMITMENT
#define NIZK2_SMALL_COMMITMENT 0
#endif
#ifndef NIZK2_DOUBLE_COMMITMENT
#define NIZK2_DOUBLE_COMMITMENT 0
#endif
#if NIZK2_SMALL_COMMITMENT && NIZK2_DOUBLE_COMMITMENT
#error "Small-message and full-message commitment modes are mutually exclusive"
#endif

#endif
