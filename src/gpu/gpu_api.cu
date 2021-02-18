/*
 *
 * gpu_api.cpp
 * PopPUNK dists using CUDA
 * gcc compiled part (uses Eigen)
 *
 */

// std
#include <cstdint>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <future>

// internal headers
#include "sketch/bitfuncs.hpp"
#include "gpu.hpp"
#include "api.hpp"
#include "database/database.hpp"
#include "sketch/sketch.hpp"

static const float mem_epsilon = 0.05;

template <class T>
inline T samples_to_rows(const T samples)
{
  return ((samples * (samples - 1)) >> 1);
}

/*
*
*  Sketch functions
*
*/

// Read in a batch of sequence data (in parallel)
std::vector<SeqBuf> read_seq_batch(
    std::vector<std::vector<std::string>>::const_iterator &file_it,
    const size_t batch_size,
    const size_t max_kmer,
    const size_t cpu_threads)
{
  std::vector<SeqBuf> seq_in_batch(batch_size);
#pragma omp parallel for schedule(static) num_threads(cpu_threads)
  for (size_t j = 0; j < batch_size; j++)
  {
    seq_in_batch[j] = SeqBuf(*(file_it + j), max_kmer);
  }
  file_it += batch_size;
  return seq_in_batch;
}

inline size_t cap_batch_size(const size_t idx, const size_t total_size,
                             size_t batch_size)
{
  if (idx + batch_size >= total_size)
  {
    batch_size = total_size - idx;
  }
  return (batch_size);
}

std::vector<Reference> create_sketches_cuda(const std::string &db_name,
                                            const std::vector<std::string> &names,
                                            const std::vector<std::vector<std::string>> &files,
                                            const std::vector<size_t> &kmer_lengths,
                                            const size_t sketchsize64,
                                            const bool use_rc,
                                            size_t min_count,
                                            const size_t cpu_threads,
                                            const int device_id)
{
  // Try loading sketches from file
  std::vector<Reference> sketches;
  bool resketch = true;
  if (file_exists(db_name + ".h5"))
  {
    sketches = load_sketches(db_name, names, kmer_lengths);
    if (sketches.size() == names.size())
    {
      resketch = false;
    }
  }

  if (resketch)
  {
    Database sketch_db(db_name + ".h5", false);
    sketches.resize(names.size());

    initialise_device(device_id);
    std::cerr << "Sketching " << files.size() << " read sets on GPU device " << device_id << std::endl;
    std::cerr << "also using " << cpu_threads << " CPU cores" << std::endl;

    // memory for filter and nthash only need to be allocated once
    copyNtHashTablesToDevice();
    GPUCountMin countmin_filter;
    if (min_count > std::numeric_limits<unsigned int>::max())
    {
      min_count = std::numeric_limits<unsigned int>::max();
    }

    size_t worker_threads = MAX(1, cpu_threads - 1);
    size_t n_batches = files.size() / worker_threads +
                       (files.size() % worker_threads ? 1 : 0);
    size_t batch_size = cap_batch_size(0, files.size(), worker_threads);

    // CPU threads read in sequence (asynchronously)
    std::launch policy = std::launch::async;
    if (cpu_threads == 1)
    {
      // Gives serial behaviour if only one thread available
      policy = std::launch::deferred;
    }
    auto file_it = files.cbegin();
    std::future<std::vector<SeqBuf>> seq_reader =
        std::async(std::launch::deferred, &read_seq_batch, std::ref(file_it),
                   batch_size, kmer_lengths.back(), cpu_threads);

    for (size_t i = 0; i < files.size(); i += worker_threads)
    {
      std::cerr << "Sketching batch: "
                << i / worker_threads + 1
                << " of "
                << n_batches
                << std::endl;
      batch_size = cap_batch_size(i, files.size(), worker_threads);

      // Get the next batch asynchronously
      std::vector<SeqBuf> seq_in_batch = seq_reader.get();
      if (file_it != files.cend())
      {
        seq_reader = std::async(policy, &read_seq_batch, std::ref(file_it),
                                cap_batch_size(i + worker_threads, files.size(), worker_threads),
                                kmer_lengths.back(), worker_threads);
      }

      // Run the sketch on the GPU (serially over the batch)
      for (size_t j = 0; j < batch_size; j++)
      {
        robin_hood::unordered_map<int, std::vector<uint64_t>> usigs;
        size_t seq_length;
        bool densified;
        std::tie(usigs, seq_length, densified) =
            sketch_gpu(
                seq_in_batch[j],
                countmin_filter,
                sketchsize64,
                kmer_lengths,
                def_bbits,
                use_rc,
                min_count,
                cpu_threads);

        fprintf(stderr, "%ck = %d   (100%%)\n", 13,
                static_cast<int>(kmer_lengths.back()));

        // Make Reference object, and save in HDF5 DB
        sketches[i + j] = Reference(names[i + j], usigs, def_bbits, sketchsize64,
                                    seq_length, seq_in_batch[j].get_composition(),
                                    seq_in_batch[j].missing_bases(), use_rc, densified);
        sketch_db.add_sketch(sketches[i + j]);
        if (densified)
        {
          std::cerr << "NOTE: "
                    << names[i + j]
                    << " required densification"
                    << std::endl;
        }
      }
    }
  }
  return (sketches);
}

// Gives strides aligned to the warp size (32)
inline size_t warpPad(const size_t stride)
{
  return (stride + (stride % warp_size ? warp_size - stride % warp_size : 0));
}

std::tuple<robin_hood::unordered_map<int, std::vector<uint64_t>>, size_t, bool>
sketch_gpu(
    SeqBuf &seq,
    GPUCountMin &countmin,
    const uint64_t sketchsize,
    const std::vector<size_t> &kmer_lengths,
    const size_t bbits,
    const bool use_canonical,
    const uint8_t min_count,
    const size_t cpu_threads)
{
  const uint64_t nbins = sketchsize * NBITS(uint64_t);
  const uint64_t binsize = (SIGN_MOD + nbins - 1ULL) / nbins;
  robin_hood::unordered_map<int, std::vector<uint64_t>> sketch;

  DeviceReads reads(seq, cpu_threads);

  double minhash_sum = 0;
  bool densified = false;
  for (auto k : kmer_lengths)
  {
    fprintf(stderr, "%ck = %d (%.1lf%%)", 13, static_cast<int>(k), 0.);
    std::vector<uint64_t> usigs(sketchsize * bbits, 0);
    std::vector<uint64_t> signs = get_signs(reads, countmin, k,
                                            use_canonical, min_count,
                                            binsize, nbins);

    minhash_sum += inverse_minhash(signs);

    // Apply densifying function
    densified |= densifybin(signs);
    fillusigs(usigs, signs, bbits);
    sketch[k] = usigs;
  }
  size_t seq_size = static_cast<size_t>((double)kmer_lengths.size() / minhash_sum);
  return (std::make_tuple(sketch, seq_size, densified));
}