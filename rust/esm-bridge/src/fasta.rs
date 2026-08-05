//! Minimal protein FASTA parser. Input must already be decompressed.

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// The header up to the first whitespace.
    pub id: String,
    /// The remainder of the header, if any.
    pub desc: String,
    pub seq: Vec<u8>,
}

// TODO: see if uniform with needletail, or perhaps we should reduce the imported crates for size and keep alternatives like these ones?
pub fn parse(bytes: &[u8]) -> Result<Vec<Record>, String> {
    let mut out: Vec<Record> = Vec::new();
    let mut cur: Option<Record> = None;

    for line in bytes.split(|&b| b == b'\n') {
        let line = match line.last() {
            Some(&b'\r') => &line[..line.len() - 1],
            _ => line,
        };
        if line.is_empty() {
            continue;
        }
        match line[0] {
            b'>' => {
                if let Some(r) = cur.take() {
                    out.push(r);
                }
                let header = String::from_utf8_lossy(&line[1..]).into_owned();
                let mut it = header.splitn(2, char::is_whitespace);
                let id = it.next().unwrap_or("").to_string();
                let desc = it.next().unwrap_or("").trim().to_string();
                if id.is_empty() {
                    return Err(format!(
                        "FASTA record {} has an empty identifier",
                        out.len() + 1
                    ));
                }
                cur = Some(Record {
                    id,
                    desc,
                    seq: Vec::new(),
                });
            }
            // Legacy FASTA comment line.
            b';' => continue,
            _ => match cur.as_mut() {
                Some(r) => r.seq.extend_from_slice(line),
                None => return Err("FASTA content before the first '>' header".into()),
            },
        }
    }
    if let Some(r) = cur.take() {
        out.push(r);
    }
    if out.is_empty() {
        return Err("No FASTA records found. Is this a protein FASTA file?".into());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiline_records() {
        let f = b">p1 first protein\nMKAL\nVGT\n>p2\nWWW\n";
        let r = parse(f).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].id, "p1");
        assert_eq!(r[0].desc, "first protein");
        assert_eq!(r[0].seq, b"MKALVGT");
        assert_eq!(r[1].id, "p2");
        assert_eq!(r[1].desc, "");
        assert_eq!(r[1].seq, b"WWW");
    }

    #[test]
    fn handles_crlf_and_blank_lines_and_no_trailing_newline() {
        let f = b">p1\r\nMK\r\n\r\nAL\r\n>p2\r\nWW";
        let r = parse(f).unwrap();
        assert_eq!(r[0].seq, b"MKAL");
        assert_eq!(r[1].seq, b"WW");
    }

    #[test]
    fn skips_legacy_comment_lines() {
        let r = parse(b">p1\n; a comment\nMK\n").unwrap();
        assert_eq!(r[0].seq, b"MK");
    }

    #[test]
    fn tabs_separate_id_from_description() {
        let r = parse(b">p1\tsome desc\nMK\n").unwrap();
        assert_eq!(r[0].id, "p1");
        assert_eq!(r[0].desc, "some desc");
    }

    #[test]
    fn empty_sequence_records_are_kept() {
        // A header with no residues is preserved so downstream row counts stay aligned
        // with the input file.
        let r = parse(b">p1\n>p2\nMK\n").unwrap();
        assert_eq!(r.len(), 2);
        assert!(r[0].seq.is_empty());
    }

    #[test]
    fn rejects_content_before_first_header() {
        assert!(parse(b"MKAL\n>p1\nMK\n").is_err());
    }

    #[test]
    fn rejects_empty_identifier() {
        assert!(parse(b"> \nMK\n").is_err());
        assert!(parse(b">\nMK\n").is_err());
    }

    #[test]
    fn rejects_input_with_no_records() {
        assert!(parse(b"").is_err());
        assert!(parse(b"\n\n").is_err());
    }
}
