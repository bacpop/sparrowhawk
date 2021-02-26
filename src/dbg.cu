
/*
 *
 * sketch.cu
 * CUDA version of bindash sketch method
 *
 */

#include <stdint.h>

#include "count_min.cuh"
#include "device_reads.cuh"
#include "seqio.hpp"

// Initialise device and return info on its memory
std::tuple<size_t, size_t, size_t> initialise_device(const int device_id) {
  CUDA_CALL(cudaSetDevice(device_id));

  size_t mem_free = 0;
  size_t mem_total = 0;
  CUDA_CALL(cudaMemGetInfo(&mem_free, &mem_total));
  int shared_size = 0;
  CUDA_CALL(cudaDeviceGetAttribute(
      &shared_size, cudaDevAttrMaxSharedMemoryPerBlock, device_id));
  return (std::make_tuple(mem_free, mem_total, static_cast<size_t>(shared_size)));
}

__global__ void count_kmers(char *read_seq,
                            const size_t n_reads,
                            const size_t read_length,
                            const size_t read_stride,
                            const int k,
                            unsigned int *countmin_table,
                            sparrowhawk::countmin::count_min_pars cm_pars,
                            volatile int *blocks_complete) {
  for (int read_index = blockIdx.x * blockDim.x + threadIdx.x; i < n_reads; i += blockDim.x * gridDim.x) {
    uint64_t fhVal, rhVal, hVal;

    // Load reads in block into shared memory
    extern __shared__ char read_shared[];
    auto block = cooperative_groups::this_thread_block();

    cuda::pipeline pipe;
    for (int base_idx = 0; base_idx < read_length; base_idx++) {
      memcpy_async(read_shared[threadIdx.x + base_idx * blockDim.x],
                   read_seq[read_index + base_idx * read_stride], pipe);
    }
    pipe.commit();
    pipe.wait_prior<0>();

    // Get first valid k-mer
    NTC64(read_shared, k, fhVal, rhVal, hVal, blockDim.x);
    sparrowhawk::countmin::probe(countmin_table, hVal, cm_pars, k, true)

    // Roll through remaining k-mers in the read
    for (int pos = 0; pos < read_length - k; pos++) {
      fhVal = NTF64(fhVal, k, read_shared[pos * blockDim.x],
                    read_shared[(pos + k) * blockDim.x]);
      rhVal = NTR64(rhVal, k, read_shared[pos * blockDim.x],
                    read_shared[(pos + k) * blockDim.x]);
      hVal = (rhVal < fhVal) ? rhVal : fhVal;
      sparrowhawk::countmin::probe(countmin_table, hVal, cm_pars, k, true)
    }

    // update progress meter
    update_progress(read_index, n_reads, blocks_complete);
    __syncwarp();
  }
}

// Call from python .make_dbg
void dbg_wrapper(const std::string &fwd_reads,
                 const std::string &rev_reads,
                 const int k,
                 const uint16_t min_count,
                 const int n_threads) {
  sequence = SeqBuf({fwd_reads, rev_reads}, const size_t kmer_len);
  sequence_device = sparrowhawk::device_reads(sequence, n_threads);
  filter = sparrowhawk::countmin::count_min_filter(27, 2, 4); // hardcoding pars for now
  build_dbg(sequence_device, filter, k, min_count);
}

void build_dbg(sparrowhawk::device_reads &reads,
               sparrowhawk::countmin::count_min_filter &countmin,
               const int k,
               const uint16_t min_count) {
  // Progress meter
  volatile int *blocks_complete;
  CUDA_CALL(cudaMallocManaged(&blocks_complete, sizeof(int)));
  *blocks_complete = 0;

  // Set countmin to zero (already on device)
  countmin.reset();

  // Run process_read kernel
  //      This runs nthash on read sequence at all k-mer lengths
  //      Check vs signs and countmin on whether to add each
  //      (get this working for a single k-mer length first)
  const size_t blockSize = 64;
  const size_t blockCount = (reads.count() + blockSize - 1) / blockSize;
  count_kmers<<<blockCount, blockSize,
                reads.length() * blockSize * sizeof(char)>>>(
      reads.read_ptr(), reads.count(), reads.length(), reads.stride(), k,
      countmin.get_table(), min_count,
      blocks_complete);

  CUDA_CALL(cudaGetLastError());
  sparrowhawk::utils::report_progress(blocks_complete, reads.count());
  CUDA_CALL(cudaDeviceSynchronize());
  fprintf(stderr, "%c100%%", 13);

  return (signs);
}

