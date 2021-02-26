
/*
 *
 * sketch.cu
 * CUDA version of bindash sketch method
 *
 */

#include <stdint.h>

#include "cuda.cuh"
#include "gpu.hpp"

// nthash
#include "hash.hpp"


// hash iterator object
__global__ void process_reads(char *read_seq, const size_t n_reads,
                              const size_t read_length,
                              const size_t read_stride, const int k,
                              uint64_t *signs, const uint64_t binsize,
                              unsigned int *countmin_table, const bool use_rc,
                              const uint16_t min_count,
                              volatile int *blocks_complete) {
  int read_index = blockIdx.x * blockDim.x + threadIdx.x;
  if (read_index < n_reads) {
    uint64_t fhVal, rhVal, hVal;

    // Load reads in block into shared memory
    extern __shared__ char read_shared[];
    for (int base_idx = 0; base_idx < read_length; base_idx++) {
      read_shared[threadIdx.x + base_idx * blockDim.x] =
          read_seq[read_index + base_idx * read_stride];
    }
    __syncthreads();

    // Get first valid k-mer
    if (use_rc) {
      NTC64(read_shared, k, fhVal, rhVal, hVal, blockDim.x);
      binhash(signs, countmin_table, hVal, binsize, k, min_count);
    } else {
      NT64(read_shared, k, fhVal, blockDim.x);
      binhash(signs, countmin_table, hVal, binsize, k, min_count);
    }

    // Roll through remaining k-mers in the read
    for (int pos = 0; pos < read_length - k; pos++) {
      fhVal = NTF64(fhVal, k, read_shared[pos * blockDim.x],
                    read_shared[(pos + k) * blockDim.x]);
      if (use_rc) {
        rhVal = NTR64(rhVal, k, read_shared[pos * blockDim.x],
                      read_shared[(pos + k) * blockDim.x]);
        hVal = (rhVal < fhVal) ? rhVal : fhVal;
        binhash(signs, countmin_table, hVal, binsize, k, min_count);
      } else {
        binhash(signs, countmin_table, fhVal, binsize, k, min_count);
      }
    }

    // update progress meter
    update_progress(read_index, n_reads, blocks_complete);
    __syncwarp();
  }
}



// main function called here returns signs vector - rest can be done by
// sketch.cpp
std::vector<uint64_t>
get_signs(DeviceReads &reads, // use seqbuf.as_square_array() to get this
          GPUCountMin &countmin, const int k, const bool use_rc,
          const uint16_t min_count, const uint64_t binsize,
          const uint64_t nbins) {
  // Progress meter
  volatile int *blocks_complete;
  CUDA_CALL(cudaMallocManaged(&blocks_complete, sizeof(int)));
  *blocks_complete = 0;

  // Set countmin to zero (already on device)
  countmin.reset();

  // Signs
  std::vector<uint64_t> signs(nbins, UINT64_MAX);
  uint64_t *d_signs;
  CUDA_CALL(cudaMalloc((void **)&d_signs, nbins * sizeof(uint64_t)));
  CUDA_CALL(cudaMemcpy(d_signs, signs.data(), nbins * sizeof(uint64_t),
                       cudaMemcpyDefault));

  // Run process_read kernel
  //      This runs nthash on read sequence at all k-mer lengths
  //      Check vs signs and countmin on whether to add each
  //      (get this working for a single k-mer length first)
  const size_t blockSize = 64;
  const size_t blockCount = (reads.count() + blockSize - 1) / blockSize;
  process_reads<<<blockCount, blockSize,
                  reads.length() * blockSize * sizeof(char)>>>(
      reads.read_ptr(), reads.count(), reads.length(), reads.stride(), k,
      d_signs, binsize, countmin.get_table(), use_rc, min_count,
      blocks_complete);

  CUDA_CALL(cudaGetLastError());
  reportSketchProgress(blocks_complete, k, reads.count());

  // Copy signs back from device
  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaMemcpy(signs.data(), d_signs, nbins * sizeof(uint64_t),
                       cudaMemcpyDefault));
  CUDA_CALL(cudaFree(d_signs));

  fprintf(stderr, "%ck = %d (100%%)", 13, k);

  return (signs);
}

