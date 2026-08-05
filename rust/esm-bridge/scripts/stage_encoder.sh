#!/usr/bin/env bash
# Copy the trained UMAP encoder into www/public/, gzipped.
#
# Produce the input first with `train_umap embed`, scripts/umap_coords.py, then
# `train_umap fit`. Roughly 460 KB raw.

set -euo pipefail
cd "$(dirname "$0")/.."

SRC="${1:-model/generated/esm2_umap_encoder.bpk}"
DEST="../../www/public/esm2_umap_encoder.bpkz"

if [ ! -f "$SRC" ]; then
    echo "$SRC not found; produce it with the three stages first:" >&2
    echo "  cargo run --release --bin train_umap -- embed --fasta-list <corpus.tsv> --work-dir <dir>" >&2
    echo "  scripts/umap_coords.py --work-dir <dir>" >&2
    echo "  cargo run --release --bin train_umap -- fit --work-dir <dir> --out $SRC" >&2
    exit 1
fi

mkdir -p "$(dirname "$DEST")"
gzip -9 -c "$SRC" > "$DEST.tmp"
mv "$DEST.tmp" "$DEST"

printf 'staged %s\n  raw     %s\n  gzipped %s\n' \
    "$DEST" "$(du -h "$SRC" | cut -f1)" "$(du -h "$DEST" | cut -f1)"
