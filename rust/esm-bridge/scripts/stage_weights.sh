#!/usr/bin/env bash
# Copy the burnpack weights produced by build.rs into www/public/, as f16, gzipped.
#
# Deliberately NOT embedded in the .wasm: build.rs uses LoadStrategy::Bytes and the worker
# fetches this asset at runtime, as the AMR tab does its .amridx.gz index.
#
# f16 is lossless for this model, not a trade — `halve_weights` re-checks that on every run
# and refuses to write if embeddings differ. Loading detects the dtype, so f32 still works.
#
# Sizes: 29.7 MB f32 raw -> 14.9 MB f16 raw -> 13.7 MB gzipped.

set -euo pipefail
cd "$(dirname "$0")/.."

SRC="model/generated/esm2_t6_8M_UR50D.bpk"
HALF="model/generated/esm2_t6_8M_UR50D.f16.bpk"
DEST="../../www/public/esm2_t6_8M_UR50D.bpkz"

if [ ! -f "$SRC" ]; then
    echo "$SRC not found; building to generate it..."
    cargo build --release
fi

# Regenerates and self-verifies; exits non-zero if f16 ever stops being lossless.
cargo run --release --bin halve_weights -- "$SRC" "$HALF"

mkdir -p "$(dirname "$DEST")"
gzip -9 -c "$HALF" > "$DEST.tmp"
mv "$DEST.tmp" "$DEST"

printf 'staged %s\n  f32 raw %s\n  f16 raw %s\n  gzipped %s\n' \
    "$DEST" \
    "$(du -h "$SRC" | cut -f1)" \
    "$(du -h "$HALF" | cut -f1)" \
    "$(du -h "$DEST" | cut -f1)"
