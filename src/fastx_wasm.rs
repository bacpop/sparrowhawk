//! Common parsing functions for reading fasta or fastq files, taken from DATACIN
//! https://github.com/bacpop/DATACIN

use flate2::read::MultiGzDecoder;
use seq_io::fasta::Reader as FastaReader;
use seq_io::fastq::Reader as FastqReader;
use std::io::{self, Chain, Cursor, Read};

const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];

// TODO -- would be nice to improve this so the fasta/fastq is transparent
// like gzipped/not. See https://github.com/onecodex/needletail/blob/master/src/parser/mod.rs


/// Enum that allows for alternating between uncompressed and compressed fasta and fastq.
pub enum ReaderEnum<'a, F: Read + 'a> {
    /// Uncompressed fasta/fastq
    Plain(Chain<Cursor<[u8; 2]>, &'a mut F>),

    /// g-zipped compressed fasta/fastq
    Gzipped(MultiGzDecoder<Chain<Cursor<[u8; 2]>, &'a mut F>>),
}


impl<'a, F: Read + 'a> Read for ReaderEnum<'a, F> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            ReaderEnum::Plain(reader)   => reader.read(buf),
            ReaderEnum::Gzipped(reader) => reader.read(buf),
        }
    }
}


/// Returns a reader from a fasta file
pub fn open_fasta<'a, F>(file_in: &'a mut F) -> FastaReader<ReaderEnum<'a, F>>
where
    F: Read + 'a,
{
    let mut first_two_bytes = [0; 2];
    file_in
        .read_exact(&mut first_two_bytes)
        .expect("Empty input file");
    let first_two_cursor = Cursor::new(first_two_bytes);
    let new_reader = first_two_cursor.chain(file_in);
    match first_two_bytes {
        GZ_MAGIC => {
            let gz_reader = MultiGzDecoder::new(new_reader);
            FastaReader::new(ReaderEnum::Gzipped(gz_reader))
        }
        _ => FastaReader::new(ReaderEnum::Plain(new_reader)),
    }
}


/// Returns a reader from a fastq file
pub fn open_fastq<'a, F>(file_in: &'a mut F) -> FastqReader<ReaderEnum<'a, F>>
where
    F: Read + 'a,
{
    let mut first_two_bytes = [0; 2];
    file_in
        .read_exact(&mut first_two_bytes)
        .expect("Empty input file");
    let first_two_cursor = Cursor::new(first_two_bytes);
    let new_reader = first_two_cursor.chain(file_in);
    match first_two_bytes {
        GZ_MAGIC => {
            let gz_reader = MultiGzDecoder::new(new_reader);
            FastqReader::new(ReaderEnum::Gzipped(gz_reader))
        }
        _ => FastqReader::new(ReaderEnum::Plain(new_reader)),
    }
}
