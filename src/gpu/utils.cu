/*
 *
 * utils.cu
 * CUDA utility functions
 *
 */

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

// Get the blockSize and blockCount for CUDA call
std::tuple<size_t, size_t> getBlockSize(const size_t ref_samples,
                                        const size_t query_samples,
                                        const size_t dist_rows,
                                        const bool self) {
  // Each block processes a single query. As max size is 512 threads
  // per block, may need multiple blocks (non-exact multiples lead
  // to some wasted computation in threads)
  // We take the next multiple of 32 that is larger than the number of
  // reference sketches, up to a maximum of 512
  size_t blockSize =
      std::min(512, 32 * static_cast<int>((ref_samples + 32 - 1) / 32));
  size_t blockCount = 0;
  if (self) {
    for (int i = 0; i < ref_samples; i++) {
      blockCount += (ref_samples + blockSize - 2 - i) / blockSize;
    }
  } else {
    size_t blocksPerQuery = (ref_samples + blockSize - 1) / blockSize;
    blockCount = blocksPerQuery * query_samples;
  }
  return (std::make_tuple(blockSize, blockCount));
}

// Writes a progress meter using the device int which keeps
// track of completed jobs
void reportDistProgress(volatile int *blocks_complete, long long dist_rows) {
  long long progress_blocks = 1 << progressBitshift;
  int now_completed = 0;
  float kern_progress = 0;
  if (dist_rows > progress_blocks) {
    while (now_completed < progress_blocks - 1) {
      if (*blocks_complete > now_completed) {
        now_completed = *blocks_complete;
        kern_progress = now_completed / (float)progress_blocks;
        fprintf(stderr, "%cProgress (GPU): %.1lf%%", 13, kern_progress * 100);
      } else {
        usleep(1000);
      }
    }
  }
}

// Initialise device and return info on its memory
std::tuple<size_t, size_t, size_t> initialise_device(const int device_id) {
  CUDA_CALL(cudaSetDevice(device_id));
  CUDA_CALL(cudaDeviceReset());

  size_t mem_free = 0;
  size_t mem_total = 0;
  CUDA_CALL(cudaMemGetInfo(&mem_free, &mem_total));
  int shared_size = 0;
  CUDA_CALL(cudaDeviceGetAttribute(
      &shared_size, cudaDevAttrMaxSharedMemoryPerBlock, device_id));
  return (
      std::make_tuple(mem_free, mem_total, static_cast<size_t>(shared_size)));
}


