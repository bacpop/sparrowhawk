use flate2::{read::DeflateDecoder, write::DeflateEncoder, Compression};
use std::io::{self, Read, Write};

// Max uncompressed bytes per BGZF block
const BGZF_BLOCK_SIZE: usize = 0xff00; // 65280

// BGZF block header template (18 bytes)
// Bytes 16–17 are BSIZE placeholder (total block size − 1), filled per block
const HEADER_TEMPLATE: [u8; 18] = [
    0x1f, 0x8b, 0x08, 0x04, // magic, method, FLG=FEXTRA
    0, 0, 0, 0, // MTIME
    0, 0xff, // XFL, OS=255 (unknown)
    0x06, 0x00, // XLEN=6
    b'B', b'C', 0x02, 0x00, // BC subfield id, length=2
    0, 0, // BSIZE placeholder (bytes 16–17)
];

// Standard BGZF EOF marker block (28 bytes)
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00,
    0x03, 0x00, 0, 0, 0, 0, 0, 0, 0, 0,
];

// ─────────────────────────────────────────────────────────────────────────────
// BgzfWriter
// ─────────────────────────────────────────────────────────────────────────────

pub struct BgzfWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    /// Compressed bytes written to inner so far.
    block_address: u64,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W) -> Self {
        BgzfWriter {
            inner,
            buf: Vec::with_capacity(BGZF_BLOCK_SIZE),
            block_address: 0,
        }
    }

    /// Virtual offset of the start of the next (unwritten) block.
    /// Between flushes the intra-block offset is always 0.
    pub fn virtual_offset(&self) -> u64 {
        self.block_address << 16
    }

    /// Compress and emit the current buffer as one BGZF block, then clear buf.
    fn flush_block(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }

        let crc = crc32fast::hash(&self.buf);
        let isize = self.buf.len() as u32;

        // Try deflate compression
        let compressed = {
            let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
            enc.write_all(&self.buf)?;
            enc.finish()?
        };

        // Total block size = 18 (header) + compressed_data + 8 (footer)
        // If it doesn't fit, fall back to a stored (non-compressed) block.
        let compressed_data: Vec<u8> = if compressed.len() + 26 > 65536 {
            // RFC 1951 stored block: [0x01][len_le][~len_le][data]
            let len = self.buf.len() as u16;
            let mut stored = Vec::with_capacity(5 + self.buf.len());
            stored.push(0x01); // BFINAL=1, BTYPE=00 (stored)
            stored.extend_from_slice(&len.to_le_bytes());
            stored.extend_from_slice(&(!len).to_le_bytes());
            stored.extend_from_slice(&self.buf);
            stored
        } else {
            compressed
        };

        // total = 18 header + data + 4 crc + 4 isize = data.len() + 26
        let total = compressed_data.len() + 26;
        debug_assert!(total <= 65536, "BGZF block exceeds 65536 bytes");

        let mut block = Vec::with_capacity(total);
        block.extend_from_slice(&HEADER_TEMPLATE);
        // BSIZE = total − 1 (little-endian u16 at bytes 16–17)
        let bsize = (total - 1) as u16;
        block[16] = bsize as u8;
        block[17] = (bsize >> 8) as u8;

        block.extend_from_slice(&compressed_data);
        block.extend_from_slice(&crc.to_le_bytes());
        block.extend_from_slice(&isize.to_le_bytes());

        self.inner.write_all(&block)?;
        self.block_address += block.len() as u64;
        self.buf.clear();
        Ok(())
    }

    /// Flush remaining data and append the BGZF EOF marker, returning the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.flush_block()?;
        self.inner.write_all(&EOF_BLOCK)?;
        Ok(self.inner)
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let mut written = 0;
        let mut remaining = data;
        while !remaining.is_empty() {
            let space = BGZF_BLOCK_SIZE - self.buf.len();
            let take = remaining.len().min(space);
            self.buf.extend_from_slice(&remaining[..take]);
            remaining = &remaining[take..];
            written += take;
            if self.buf.len() >= BGZF_BLOCK_SIZE {
                self.flush_block()?;
            }
        }
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.inner.flush()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BgzfReader
// ─────────────────────────────────────────────────────────────────────────────

pub struct BgzfReader<R: Read> {
    inner: R,
    /// Compressed bytes consumed so far (= start of the *next* unread block).
    pub block_address: u64,
    /// Compressed start offset of the block currently loaded in `block`.
    cur_block_start: u64,
    /// Decompressed contents of the current block.
    block: Vec<u8>,
    /// Read position within block.
    pos: usize,
    /// (compressed_offset, cumulative_uncompressed_offset) pairs – one per block.
    pub gzi: Vec<(u64, u64)>,
    /// Cumulative uncompressed bytes before the current block.
    pub uncompressed_addr: u64,
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R) -> Self {
        BgzfReader {
            inner,
            block_address: 0,
            cur_block_start: 0,
            block: Vec::new(),
            pos: 0,
            gzi: Vec::new(),
            uncompressed_addr: 0,
        }
    }

    /// Current virtual offset: (start_of_current_block << 16) | pos
    pub fn virtual_offset(&self) -> u64 {
        (self.cur_block_start << 16) | (self.pos as u64)
    }

    /// Current position as a plain uncompressed byte offset.
    pub fn uncompressed_offset(&self) -> u64 {
        self.uncompressed_addr - self.block.len() as u64 + self.pos as u64
    }

    /// GZI block entries collected during reading.
    pub fn gzi_entries(&self) -> &[(u64, u64)] {
        &self.gzi
    }

    /// Read and decompress the next BGZF block.
    /// Returns Ok(false) on clean EOF (empty read of header), Ok(true) on success.
    fn read_block(&mut self) -> io::Result<bool> {
        let caddr_before = self.block_address;
        let uaddr_before = self.uncompressed_addr;

        let mut header = [0u8; 18];
        match self.inner.read(&mut header[..1]) {
            Ok(0) => return Ok(false), // clean EOF
            Ok(_) => {}
            Err(e) => return Err(e),
        }
        read_exact_inner(&mut self.inner, &mut header[1..])?;

        // Validate magic and flags
        if header[0] != 0x1f || header[1] != 0x8b {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a gzip stream",
            ));
        }
        if header[2] != 0x08 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported gzip method",
            ));
        }

        let bsize = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        let deflate_len = bsize
            .checked_sub(26)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BGZF block too small"))?;

        let mut deflate_data = vec![0u8; deflate_len];
        read_exact_inner(&mut self.inner, &mut deflate_data)?;

        let mut footer = [0u8; 8];
        read_exact_inner(&mut self.inner, &mut footer)?;

        let expected_crc = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
        let expected_isize =
            u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;

        // Decompress
        self.block.clear();
        self.block.reserve(expected_isize);
        let mut dec = DeflateDecoder::new(&deflate_data[..]);
        dec.read_to_end(&mut self.block)?;

        if self.block.len() != expected_isize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "BGZF isize mismatch: got {} expected {}",
                    self.block.len(),
                    expected_isize
                ),
            ));
        }

        let actual_crc = crc32fast::hash(&self.block);
        if actual_crc != expected_crc {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "BGZF CRC32 mismatch",
            ));
        }

        self.cur_block_start = caddr_before;
        self.block_address += bsize as u64;
        self.pos = 0;

        // Record GZI entry: skip the implicit (0,0) first-block entry and skip
        // the empty EOF block (isize==0).
        if !self.block.is_empty() && (caddr_before > 0 || uaddr_before > 0) {
            self.gzi.push((caddr_before, uaddr_before));
        }
        self.uncompressed_addr += self.block.len() as u64;

        Ok(true)
    }

    /// Read bytes until `\n` (inclusive), appending to `buf`.
    /// Returns `(bytes_read, voff_at_line_start)`.
    /// Returns `(0, voff)` on EOF.
    pub fn read_line(&mut self, buf: &mut Vec<u8>) -> io::Result<(usize, u64)> {
        let voff_start = self.virtual_offset();
        if self.pos >= self.block.len() {
            let got = self.read_block()?;
            if !got || self.block.is_empty() {
                return Ok((0, voff_start));
            }
        }
        let mut total = 0usize;
        loop {
            // Scan for newline in current block
            let slice = &self.block[self.pos..];
            match slice.iter().position(|&b| b == b'\n') {
                Some(nl) => {
                    let end = nl + 1;
                    buf.extend_from_slice(&slice[..end]);
                    self.pos += end;
                    total += end;
                    if self.pos == self.block.len() {
                        let _ = self.read_block();
                    }
                    return Ok((total, voff_start));
                }
                None => {
                    buf.extend_from_slice(slice);
                    total += slice.len();
                    self.pos = self.block.len();
                    // Load next block and continue
                    let got = self.read_block()?;
                    if !got || self.block.is_empty() {
                        return Ok((total, voff_start));
                    }
                }
            }
        }
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        // Refill if needed
        while self.pos >= self.block.len() {
            let got = self.read_block()?;
            if !got || self.block.is_empty() {
                return Ok(0);
            }
        }
        let avail = &self.block[self.pos..];
        let take = avail.len().min(out.len());
        out[..take].copy_from_slice(&avail[..take]);
        self.pos += take;
        Ok(take)
    }
}

/// Like `read_exact` but for our inner reader (avoids naming conflicts).
fn read_exact_inner<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<()> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "unexpected EOF",
                ))
            }
            Ok(n) => filled += n,
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience function
// ─────────────────────────────────────────────────────────────────────────────

/// Compress all bytes from `input` into BGZF format, writing to `output`.
pub fn bgzf_compress<R: Read, W: Write>(mut input: R, output: W) -> io::Result<()> {
    let mut writer = BgzfWriter::new(output);
    let mut buf = vec![0u8; 65536];
    loop {
        let n = input.read(&mut buf)?;
        if n == 0 {
            break;
        }
        writer.write_all(&buf[..n])?;
    }
    writer.finish()?;
    Ok(())
}
