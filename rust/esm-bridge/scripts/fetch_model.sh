#!/usr/bin/env bash
# Obtain the ESM-2 (8M) ONNX graph that build.rs feeds to burn-onnx.
#
# Xue-Jun/esm2-8M `onnx/model.onnx`: opset 18, so LayerNormalization is a single node, and a
# feature-extraction export, so the output is last_hidden_state[batch,seq,320] not the LM head.
#
# Provenance is checked by test, not asserted here: `matches_official_meta_weights` compares
# our embeddings against facebook/esm2_t6_8M_UR50D run in torch (scripts/reference_official.py).
#
# Re-exporting from the original weights, if ever needed (heavier: pulls in torch):
#   pip install "optimum[exporters]" transformers torch onnx
#   optimum-cli export onnx --model facebook/esm2_t6_8M_UR50D \
#       --task feature-extraction --opset 18 --monolith model/onnx_export/
#   mv model/onnx_export/model.onnx model/esm2_t6_8M_UR50D.onnx

set -euo pipefail
cd "$(dirname "$0")/.."

URL="https://huggingface.co/Xue-Jun/esm2-8M/resolve/main/onnx/model.onnx"
DEST="model/esm2_t6_8M_UR50D.onnx"
SHA256="848b8f09bec0bb610b5a4bb4a739e96418934b475faf276f91b35cefbb4d0b36"

mkdir -p model
if [ -f "$DEST" ] && echo "$SHA256  $DEST" | sha256sum -c --status; then
    echo "$DEST already present and matches the expected checksum."
    exit 0
fi

echo "Downloading $URL ..."
curl -sSfL -o "$DEST.tmp" "$URL"
echo "$SHA256  $DEST.tmp" | sed "s|$DEST.tmp|$DEST.tmp|" | sha256sum -c -
mv "$DEST.tmp" "$DEST"
echo "Wrote $DEST"
