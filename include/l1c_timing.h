#ifndef _L1C_TIMING_
#define _L1C_TIMING_
#include "config.h"

#include <sys/time.h>
#include <unistd.h>

static inline struct timeval l1c_get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);

  return tv;
}

static inline double l1c_get_time_diff(struct timeval tv_start, struct timeval tv_end) {

  return (double)(tv_end.tv_usec - tv_start.tv_usec) / 1000000.0 +
         (double)(tv_end.tv_sec - tv_start.tv_sec);
}
#endif
