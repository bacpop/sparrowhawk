#!/usr/bin/env python3
"""Generate reference ESM-2 embeddings with onnxruntime, for comparison with the burn port.

This runs the *same* ONNX graph that build.rs feeds to burn-onnx, so a match proves that
burn-onnx imported the graph faithfully. It deliberately does not need PyTorch.

That is a narrower claim than it may look: it cannot tell a good export from a bad one,
because the reference and the model under test share the same weights. For the check against
Meta's published weights, see scripts/reference_official.py.

The tokenizer here is re-derived from the published ESM alphabet rather than shared with
the Rust implementation, so the two are genuine cross-checks of each other.

Usage:
    ./onnxenv/bin/python scripts/reference.py \
        --onnx model/esm2_t6_8M_UR50D.onnx \
        --fasta tests/data/seqs.faa \
        --out tests/data/reference.npy
"""

import argparse
import numpy as np
import onnxruntime as ort

# Fixed 33-token ESM-1b/ESM-2 alphabet.
VOCAB = (
    "<cls> <pad> <eos> <unk> L A G V S E R T I D P K Q N F Y M H W C X B U Z O . - "
    "<null_1> <mask>"
).split()
TOK = {t: i for i, t in enumerate(VOCAB)}
CLS, PAD, EOS, UNK = TOK["<cls>"], TOK["<pad>"], TOK["<eos>"], TOK["<unk>"]
MAX_RESIDUES = 1022


def parse_fasta(path):
    records, ident, desc, seq = [], None, "", []
    for raw in open(path, "rb").read().decode("utf-8", "replace").splitlines():
        line = raw.rstrip("\r")
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if ident is not None:
                records.append((ident, desc, "".join(seq)))
            head = line[1:]
            parts = head.split(None, 1)
            ident = parts[0] if parts else ""
            desc = parts[1].strip() if len(parts) > 1 else ""
            seq = []
        else:
            seq.append(line)
    if ident is not None:
        records.append((ident, desc, "".join(seq)))
    return records


def residues(seq, max_residues=MAX_RESIDUES):
    s = "".join(seq.split())          # drop all whitespace
    if s.endswith("*"):               # drop one trailing stop codon
        s = s[:-1]
    return s.upper()[:max_residues]   # uppercase, then head-truncate


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--onnx", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=MAX_RESIDUES)
    ap.add_argument("--batch-size", type=int, default=8)
    args = ap.parse_args()

    records = parse_fasta(args.fasta)
    print(f"sequences : {len(records)}")

    sess = ort.InferenceSession(args.onnx, providers=["CPUExecutionProvider"])
    out_name = sess.get_outputs()[0].name
    print(f"onnx out  : {out_name} {sess.get_outputs()[0].shape}")

    dim = sess.get_outputs()[0].shape[-1]
    embeddings = np.zeros((len(records), dim), dtype=np.float32)

    for start in range(0, len(records), args.batch_size):
        chunk = records[start:start + args.batch_size]
        ids, att, pool = encode([r[2] for r in chunk], args.max_len)
        hidden = sess.run([out_name], {"input_ids": ids, "attention_mask": att})[0]
        denom = np.maximum(pool.sum(axis=1, keepdims=True), 1.0)
        embeddings[start:start + len(chunk)] = (
            (hidden * pool[..., None]).sum(axis=1) / denom
        ).astype(np.float32)

    np.save(args.out, embeddings)
    print(f"wrote     : {args.out} {embeddings.shape}")
    for (ident, _, seq), vec in zip(records, embeddings):
        print(f"  {ident:<40} len={len(residues(seq, args.max_len)):<5} "
              f"|v|={np.linalg.norm(vec):.4f} v[0:3]={np.round(vec[:3], 4).tolist()}")


if __name__ == "__main__":
    main()
