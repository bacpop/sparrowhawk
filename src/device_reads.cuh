class device_reads
{
public:
  device_reads(const SeqBuf &seq_in, const size_t n_threads)
    : _n_reads(seq_in.n_full_seqs()),
      _read_length(seq_in.max_length()),
      _read_stride(seq_in.n_full_seqs_padded()) {
  std::vector<char> flattened_reads = seq_in.as_square_array(n_threads);
  CUDA_CALL(cudaMalloc((void **)&_d_reads, flattened_reads.size() * sizeof(char)));
  CUDA_CALL(cudaMemcpy(_d_reads, flattened_reads.data(),
                       flattened_reads.size() * sizeof(char),
                       cudaMemcpyDefault));
}
  ~device_reads() {
    CUDA_CALL(cudaFree(_d_reads));
  }

  char *read_ptr() { return _d_reads; }
  size_t count() const { return _n_reads; }
  size_t length() const { return _read_length; }
  size_t stride() const { return _read_stride; }

private:
  // delete move and copy to avoid accidentally using them
  device_reads(const device_reads &) = delete;
  device_reads(device_reads &&) = delete;

  char *_d_reads;
  size_t _n_reads;
  size_t _read_length;
  size_t _read_stride;
};



