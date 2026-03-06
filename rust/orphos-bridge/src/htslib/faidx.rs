use std::io::{self, Read, Write};
use super::bgzf::BgzfReader;

// ---------------------------------------------------------------
// Public API
// ---------------------------------------------------------------

/// Build faidx (`.fai`) and GZI (`.gzi`) indexes for a BGZF-compressed FASTA.
///
/// - `bgzf_input`: BGZF-compressed FASTA byte stream.
/// - `fai_output`: receives the text `.fai` index.
/// - `gzi_output`: receives the binary `.gzi` block index.
pub fn faidx_index_fasta<R: Read, F: Write, G: Write>(
    bgzf_input: R,
    mut fai_output: F,
    mut gzi_output: G,
) -> io::Result<()> {
    let mut reader = BgzfReader::new(bgzf_input);

    // State for current sequence
    let mut cur_name: Option<String> = None;
    let mut cur_seq_offset: u64 = 0; // uncompressed byte position of first base
    let mut cur_seq_len: u64 = 0;
    let mut cur_line_blen: usize = 0; // raw bytes per line (including newline)
    let mut cur_line_len: usize = 0;  // bases per line (excluding newline)
    let mut first_data_line: bool = false;

    let write_record = |fai: &mut F,
                        name: &str,
                        seq_len: u64,
                        seq_offset: u64,
                        line_blen: usize,
                        line_len: usize| -> io::Result<()> {
        let line = format!("{}\t{}\t{}\t{}\t{}\n", name, seq_len, seq_offset, line_len, line_blen);
        fai.write_all(line.as_bytes())
    };

    let mut line_buf = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let (n, _voff_start) = reader.read_line(&mut line_buf)?;
        if n == 0 {
            // EOF – flush last sequence
            if let Some(ref name) = cur_name {
                write_record(
                    &mut fai_output,
                    name,
                    cur_seq_len,
                    cur_seq_offset,
                    cur_line_blen,
                    cur_line_len,
                )?;
            }
            break;
        }

        if line_buf.is_empty() || line_buf[0] == b'\n' || line_buf[0] == b'\r' {
            continue;
        }

        if line_buf[0] == b'>' {
            // Flush previous sequence
            if let Some(ref name) = cur_name {
                write_record(
                    &mut fai_output,
                    name,
                    cur_seq_len,
                    cur_seq_offset,
                    cur_line_blen,
                    cur_line_len,
                )?;
            }

            // Parse new header – name ends at first whitespace
            let header = strip_newline(&line_buf[1..]); // skip '>'
            let name_end = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(header.len());
            let name = std::str::from_utf8(&header[..name_end])
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF8 sequence name"))?
                .to_owned();

            cur_name = Some(name);
            cur_seq_len = 0;
            cur_line_blen = 0;
            cur_line_len = 0;
            // seq_offset is the uncompressed byte position of the first base
            cur_seq_offset = reader.uncompressed_offset();
            first_data_line = true;
        } else {
            // Data line
            let raw_len = line_buf.len(); // includes newline chars
            let base_count = line_buf.iter().filter(|&&b| b.is_ascii_graphic()).count();

            if first_data_line {
                cur_line_blen = raw_len;
                cur_line_len = base_count;
                first_data_line = false;
            }

            cur_seq_len += base_count as u64;
        }
    }

    // ---------------------------------------------------------------
    // Write GZI
    // ---------------------------------------------------------------
    // Format:
    //    n_blocks: u64
    //    For each block: caddr: u64, uaddr: u64
    // The implicit (0,0) block is NOT written.
    let entries = reader.gzi_entries();
    gzi_output.write_all(&(entries.len() as u64).to_le_bytes())?;
    for &(caddr, uaddr) in entries {
        gzi_output.write_all(&caddr.to_le_bytes())?;
        gzi_output.write_all(&uaddr.to_le_bytes())?;
    }

    Ok(())
}

// ---------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------

fn strip_newline(buf: &[u8]) -> &[u8] {
    let mut end = buf.len();
    while end > 0 && (buf[end - 1] == b'\n' || buf[end - 1] == b'\r') {
        end -= 1;
    }
    &buf[..end]
}
