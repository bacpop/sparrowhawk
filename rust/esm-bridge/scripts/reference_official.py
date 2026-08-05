#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "torch", "transformers", "safetensors"]
#
# [tool.uv.sources]
# torch = { index = "pytorch-cpu" }
#
# [[tool.uv.index]]
# name = "pytorch-cpu"
# url = "https://download.pytorch.org/whl/cpu"
# explicit = true
# ///
"""Reference ESM-2 embeddings from Meta's own weights.

Unlike reference.py, which runs our own ONNX graph and so cannot tell a good export from a
bad one, a match here validates the weights themselves.

The tokenizer and pooling are re-derived here rather than shared, so the two stay independent
cross-checks. They must agree on:

  * a single trailing `*` dropped, whitespace removed, residues uppercased;
  * tokens are `<cls>` + residues[:1022] + `<eos>`, padded to the batch maximum;
  * the mean is over *residue* positions only — `<cls>`, `<eos>` and `<pad>` excluded.

Output is a committed fixture, so neither the build nor CI needs this script.

Usage:
    scripts/reference_official.py --fasta tests/data/proteins20.faa \
        --out tests/data/reference_official.npy
"""

import argparse
import pathlib

import numpy as np
import torch
from transformers import AutoModel

# Pinned by commit, not by tag: a tag can be moved, and the whole point here is provenance.
MODEL = "facebook/esm2_t6_8M_UR50D"
REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"

# Fixed 33-token ESM-1b/ESM-2 alphabet, in the published order.
VOCAB = (
    "<cls> <pad> <eos> <unk> L A G V S E R T I D P K Q N F Y M H W C X B U Z O . - "
    "<null_1> <mask>"
).split()
TOK = {t: i for i, t in enumerate(VOCAB)}
CLS, PAD, EOS, UNK = TOK["<cls>"], TOK["<pad>"], TOK["<eos>"], TOK["<unk>"]
MAX_RESIDUES = 1022


def parse_fasta(path):
    records, ident, seq = [], None, []
    for raw in pathlib.Path(path).read_bytes().decode("utf-8", "replace").splitlines():
        line = raw.rstrip("\r")
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if ident is not None:
                records.append((ident, "".join(seq)))
            ident = (line[1:].split(None, 1) or [""])[0]
            seq = []
        else:
            seq.append(line)
    if ident is not None:
        records.append((ident, "".join(seq)))
    return records


def residues(seq, max_residues=MAX_RESIDUES):
    s = "".join(seq.split())
    if s.endswith("*"):
        s = s[:-1]
    return s.upper()[:max_residues]


def encode(seqs, max_residues=MAX_RESIDUES):
    kept = [residues(s, max_residues) for s in seqs]
    length = max((len(k) for k in kept), default=0) + 2
    ids = np.full((len(kept), length), PAD, dtype=np.int64)
    att = np.zeros((len(kept), length), dtype=np.int64)
    pool = np.zeros((len(kept), length), dtype=np.float32)
    for i, k in enumerate(kept):
        ids[i, 0], att[i, 0] = CLS, 1
        for j, c in enumerate(k):
            ids[i, j + 1] = TOK.get(c, UNK)
            att[i, j + 1] = 1
            pool[i, j + 1] = 1.0
        ids[i, len(k) + 1], att[i, len(k) + 1] = EOS, 1
    return ids, att, pool


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=MAX_RESIDUES)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--perturb", type=float, default=0.0,
                    help="add N(0, x) noise to every float parameter, to test if the test would work in case the weights were actually different")
    args = ap.parse_args()

    records = parse_fasta(args.fasta)
    print(f"sequences : {len(records)}")

    # add_pooling_layer=False: the pooler is absent from the checkpoint and would be randomly
    # initialised. We read last_hidden_state, so it is unused either way
    model = AutoModel.from_pretrained(
        MODEL, revision=REVISION, add_pooling_layer=False
    ).eval()
    print(f"model     : {MODEL}@{REVISION[:12]}")

    if args.perturb:
        torch.manual_seed(0)
        with torch.no_grad():
            for p in model.parameters():
                p.add_(torch.randn_like(p) * args.perturb)
        print(f"PERTURBED : N(0, {args.perturb}) added to every parameter")

    dim = model.config.hidden_size
    embeddings = np.zeros((len(records), dim), dtype=np.float32)

    with torch.no_grad():
        for start in range(0, len(records), args.batch_size):
            chunk = records[start:start + args.batch_size]
            ids, att, pool = encode([r[1] for r in chunk], args.max_len)
            hidden = model(
                input_ids=torch.from_numpy(ids),
                attention_mask=torch.from_numpy(att),
            ).last_hidden_state.numpy()
            denom = np.maximum(pool.sum(axis=1, keepdims=True), 1.0)
            embeddings[start:start + len(chunk)] = (
                (hidden * pool[..., None]).sum(axis=1) / denom
            ).astype(np.float32)

    np.save(args.out, embeddings)
    print(f"wrote     : {args.out} {embeddings.shape}")
    for (ident, seq), vec in list(zip(records, embeddings))[:5]:
        print(f"  {ident:<40} len={len(residues(seq, args.max_len)):<5} "
              f"|v|={np.linalg.norm(vec):.4f} v[0:3]={np.round(vec[:3], 4).tolist()}")


if __name__ == "__main__":
    main()
