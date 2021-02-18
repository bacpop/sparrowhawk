/*
 *
 * gpu.hpp
 * functions using CUDA
 *
 */
#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>

static const int warp_size = 32;

// Align structs
// https://stackoverflow.com/a/12779757
#if defined(__CUDACC__) // NVCC
#define ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
#define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
#define ALIGN(n) __declspec(align(n))
#else
#error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif

// Structure of flattened vectors
struct ALIGN(16) SketchStrides
{
  size_t bin_stride;
  size_t kmer_stride;
  size_t sample_stride;
  size_t sketchsize64;
  size_t bbits;
};

struct ALIGN(8) SketchSlice
{
  size_t ref_offset;
  size_t ref_size;
  size_t query_offset;
  size_t query_size;
};

// defined in utils.cu
std::tuple<size_t, size_t, size_t> initialise_device(const int device_id);

// defined in sketch.cu
class GPUCountMin
{
public:
  GPUCountMin();
  ~GPUCountMin();

  unsigned int *get_table() { return _d_countmin_table; }

  void reset();

private:
  // delete move and copy to avoid accidentally using them
  GPUCountMin(const GPUCountMin &) = delete;
  GPUCountMin(GPUCountMin &&) = delete;

  unsigned int *_d_countmin_table;

  const unsigned int _table_width_bits;
  const uint64_t _mask;
  const uint32_t _table_width;
  const int _hash_per_hash;
  const int _table_rows;
  const size_t _table_cells;
};

class DeviceReads
{
public:
  DeviceReads(const SeqBuf &seq_in, const size_t n_threads);
  ~DeviceReads();

  char *read_ptr() { return d_reads; }
  size_t count() const { return n_reads; }
  size_t length() const { return read_length; }
  size_t stride() const { return read_stride; }

private:
  // delete move and copy to avoid accidentally using them
  DeviceReads(const DeviceReads &) = delete;
  DeviceReads(DeviceReads &&) = delete;

  char *d_reads;
  size_t n_reads;
  size_t read_length;
  size_t read_stride;
};

void copyNtHashTablesToDevice();

std::vector<uint64_t> get_signs(DeviceReads &reads,
                                GPUCountMin &countmin,
                                const int k,
                                const bool use_rc,
                                const uint16_t min_count,
                                const uint64_t binsize,
                                const uint64_t nbins);
