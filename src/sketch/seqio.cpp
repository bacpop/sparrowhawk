/*
 *
 * seqio.cpp
 * Sequence reader and buffer class
 *
 */

#include "seqio.hpp"
#include <stdint.h>
KSEQ_INIT(gzFile, gzread)

// C++ headers
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <utility>

#include "bitfuncs.hpp"

SeqBuf::SeqBuf() : _N_count(0), _max_length(0), _reads(false)
{
}

SeqBuf::SeqBuf(const std::vector<std::string> &filenames, const size_t kmer_len)
    : _N_count(0), _max_length(0), _reads(false)
{
  /*
    *   Reads entire sequence to memory
    */
  size_t seq_idx = 0;
  for (auto name_it = filenames.begin(); name_it != filenames.end(); name_it++)
  {
    // from kseq.h
    gzFile fp = gzopen(name_it->c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0)
    {
      size_t seq_len = strlen(seq->seq.s);
      if (seq_len > _max_length)
      {
        _max_length = seq_len;
      }
      if (seq_len >= kmer_len)
      {
        sequence.push_back(seq->seq.s);
        bool has_N = false;
        for (char &c : sequence.back())
        {
          if (c == 'N' || c == 'n')
          {
            _N_count++;
            has_N = true;
          }
        }
        if (!has_N)
        {
          _full_index.push_back(seq_idx);
        }
        seq_idx++;
      }

      // Presence of any quality scores - assume reads as input
      if (!_reads && seq->qual.l)
      {
        _reads = true;
      }
    }

    // If put back into object, move this to destructor below
    kseq_destroy(seq);
    gzclose(fp);
  }

  this->reset();
}

SeqBuf::SeqBuf(const std::vector<std::string> &sequence_in)
    : sequence(sequence_in), _reads(false)
{
  this->reset();
}

void SeqBuf::reset()
{
  /*
    *   Returns to start of sequences
    */
  if (sequence.size() > 0)
  {
    current_seq = sequence.begin();
    next_base = current_seq->begin();
    out_base = current_seq->end();
  }
  end = false;
}

std::vector<char> SeqBuf::as_square_array(const size_t n_threads) const
{
  if (!_reads)
  {
    throw std::runtime_error(
        "Square arrays (for GPU sketches) only supported with reads as input");
  }
  else if (n_full_seqs() == 0)
  {
    throw std::runtime_error("Input contains no sequence!");
  }

  // Pad so that the stride is a multiple of 4
  // (needed to align memory addresses)
  std::vector<char> read_array(max_length() * n_full_seqs_padded());
#pragma omp parallel for simd schedule(static) num_threads(n_threads)
  for (size_t read_idx = 0; read_idx < n_full_seqs(); read_idx++)
  {
    std::string seq = sequence[_full_index[read_idx]];
    for (size_t base_idx = 0; base_idx < seq.size(); base_idx++)
    {
      read_array[read_idx + base_idx * n_full_seqs_padded()] = seq[base_idx];
    }
    for (size_t base_idx = seq.size(); base_idx < _max_length; base_idx++)
    {
      read_array[read_idx + base_idx * n_full_seqs_padded()] = 'N';
    }
  }
  return read_array;
}

bool SeqBuf::move_next(size_t word_length)
{
  /*
    *   Moves along to next character in sequence and reverse complement
    *   Loops around to next sequence if end reached
    *   Keeps track of base before k-mer length
    */
  bool next_seq = false;
  if (!end)
  {
    next_base++;

    if (next_base == current_seq->end())
    {
      current_seq++;
      next_seq = true;
      if (current_seq == sequence.end())
      {
        end = true;
      }
      else
      {
        next_base = current_seq->begin();
        out_base = current_seq->end();
      }
    }
    else
    {
      if (out_base != current_seq->end())
      {
        out_base++;
      }
      else if ((next_base - word_length) >= current_seq->begin())
      {
        out_base = current_seq->begin();
      }
    }
  }
  return next_seq;
}
