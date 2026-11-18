#ifndef NIZK2_BUDGETED_THREADING_H
#define NIZK2_BUDGETED_THREADING_H
#ifdef _OPENMP
#include <omp.h>
static inline int
nizk2_parallel_allowed(void)
{
#if defined(LAZER_AFFINE_NESTED) && LAZER_AFFINE_NESTED
  /* The reply scheduler divides the total budget before enabling level two. */
  return omp_get_active_level() < omp_get_max_active_levels()
         && omp_get_active_level() < 2 && omp_get_max_threads() > 1;
#else
  return !omp_in_parallel();
#endif
}
#endif
#endif
