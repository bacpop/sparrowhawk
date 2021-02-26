#pragma once

namespace sparrowhawk {

namespace utils {

const int warp_size = 32;
const int progress_bitshift = 10;

// Use atomic add to update a counter, so progress works regardless of
// dispatch order
__device__ inline void update_progress(long long count, long long total,
    volatile int *blocks_complete) {
  // Progress indicator
  // The >> progressBitshift is a divide by 1024 - update roughly every 0.1%
  if (count % (total >> progressBitshift) == 0) {
    atomicAdd((int *)blocks_complete, 1);
    __threadfence_system();
  }
}

// Writes a progress meter using the device int which keeps
// track of completed jobs
void report_progress(volatile int *blocks_complete, size_t total) {
  long long progress_blocks = 1 << progress_bitshift;
  int now_completed = 0;
  float kern_progress = 0;
  if (total > progress_blocks) {
    while (now_completed < progress_blocks - 1) {
      if (*blocks_complete > now_completed) {
        now_completed = *blocks_complete;
        kern_progress = now_completed / (float)progress_blocks;
        fprintf(stderr, "%cProgress: %.1lf%%", 13, kern_progress * 100);
      } else {
        usleep(1000);
      }
    }
  }
}

// Gives strides aligned to the warp size (32)
inline size_t warpPad(const size_t stride)
{
  return (stride + (stride % warp_size ? warp_size - stride % warp_size : 0));
}

}
}