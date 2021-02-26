#pragma once

#include "hash.cuh"

namespace sparrowhawk {

namespace countmin {

struct count_min_pars {
  const int width_bits;
  const int hash_per_hash;
  uint64_t mask;
  uint32_t table_width;
  int table_rows;
};

class count_min_filter
{
public:
  count_min_filter(const size_t width_bits, const size_t hash_per_hash,
                   const int table_rows) {
  _pars.width_bits = width_bits;
  _pars.hash_per_hash = hash_per_hash;
  _pars.table_rows = table_rows;
  mask = 1;
  for (size_t i = 0; i < table_width_bits - 1; ++i) {
    _mask = _mask << 1;
    mask++;
  }
  _pars.mask = mask;
  _pars.table_width = static_cast<uint32_t>(mask);

  CUDA_CALL(cudaMalloc((void **)&_d_countmin_table,
                       table_cells() * sizeof(unsigned int)));
  CUDA_CALL(cudaMalloc((void **)&_d_pars, sizeof(count_min_pars)));
  CUDA_CALL(cudaMemcpy(_d_pars, *pars, sizeof(count_min_pars), cudaMemcpyDefault));
  reset();
}

  ~count_min_filter() {
    CUDA_CALL(cudaFree(_d_countmin_table));
  }

  unsigned int *get_table() { return _d_countmin_table; }

  size_t table_cells() const { return _pars.table_rows * _pars.table_width}

  void reset() {
      CUDA_CALL(cudaMemset(_d_countmin_table, 0, table_cells() * sizeof(unsigned int)))
  }

private:
  // delete move and copy to avoid accidentally using them
  count_min_filter(const count_min_filter &) = delete;
  count_min_filter(count_min_filter &&) = delete;

  unsigned int *_d_countmin_table;

  count_min_pars _pars;
  count_min_pars * _d_pars;
};

// TODO - load d_pars into __shared__
__device__ unsigned int probe(unsigned int * table,
                              uint64_t hash_val, count_min_pars pars,
                              const int k, const bool update) {
  unsigned int min_count = UINT32_MAX;
  for (int hash_nr = 0; hash_nr < pars.table_rows; hash_nr += pars.hash_per_hash) {
    uint64_t current_hash = hash_val;
    for (uint i = 0; i < pars.hash_per_hash; i++) {
      uint32_t hash_val_masked = current_hash & pars.mask;
      cell_ptr = table + (hash_nr + i) * pars.table_width + hash_val_masked;
      unsigned int cell_count;
      if (update) {
        cell_count = atomicInc(cell_ptr, UINT32_MAX) + 1;
      } else {
        cell_count = *cell_ptr;
      }

      if (cell_count < min_count) {
        min_count = cell_count;
      }
      current_hash = current_hash >> table_width_bits;
    }
    hash_val = shifthash(hash_val, k, hash_nr / 2);
  }
  return (min_count);
}

}
}