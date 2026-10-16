#!/usr/bin/env python3
"""Rebuild tests/wasm/data/ from the site's example inputs: `python3 tests/wasm/make_examples.py`.

Downloads the site's K. pneumoniae genome, kept whole, and its smaller proteome, cut to its first 50 proteins,
and cuts mandrake's HIV alignment to its first 300 sequences. Only needed to change the data: the WASM checks
themselves never use the network. Standard library only.
"""

import gzip
import io
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "wasm" / "data"
MANDRAKE_HIV = ROOT / "rust" / "mandrake" / "tests" / "fixtures" / "sub5k_hiv_refs_prrt_trim.fas.gz"
GENOME = "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true"
PROTEOME = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000807/UP000000807_243273.fasta.gz"
PROTEINS = 50
HIV_SEQUENCES = 300


def download(url: str, path: Path) -> None:
    print(f"downloading {path.name}", flush=True)
    partial = path.with_name(f"{path.name}.part")
    urllib.request.urlretrieve(url, partial)
    partial.rename(path)


def first_records(source: Path, count: int, output: Path) -> None:
    """The first `count` records of the gzipped FASTA `source`, gzipped without a timestamp, so that rebuilding
    unchanged data changes no bytes."""
    gz = gzip.GzipFile(output, "wb", compresslevel=9, mtime=0)
    with gzip.open(source, "rt", encoding="utf-8") as handle, io.TextIOWrapper(gz, encoding="utf-8") as out:
        records = 0
        for line in handle:
            if line.startswith(">"):
                records += 1
                if records > count:
                    break
            out.write(line)
    print(f"{output.name}: the first {count} records of {source.name}")


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    download(GENOME, DATA / "GCA_004138665.1.fasta.gz")
    proteome = DATA / "UP000000807_243273.fasta.gz"
    download(PROTEOME, proteome)
    first_records(proteome, PROTEINS, DATA / f"UP000000807_{PROTEINS}.fasta.gz")
    proteome.unlink()
    first_records(MANDRAKE_HIV, HIV_SEQUENCES, DATA / f"hiv_{HIV_SEQUENCES}.fas.gz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
