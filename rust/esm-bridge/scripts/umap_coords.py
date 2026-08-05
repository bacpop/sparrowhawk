#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "umap-learn>=0.5.7", "numba>=0.60"]
# ///
"""2-D UMAP coordinates for the corpus embedded by `train_umap embed`.

Stage 2 of 3. umap-learn is the reference implementation, so the layout the encoder is
fitted to is the one everything else is validated against.

Usage:
    scripts/umap_coords.py --work-dir model/generated/umap_corpus
    cargo run --release --bin train_umap -- fit --work-dir model/generated/umap_corpus
"""

import argparse
import pathlib
import time

import numpy as np
import umap

HIDDEN = 320


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work-dir", type=pathlib.Path, required=True)
    ap.add_argument("--n-neighbours", type=int, default=15)
    ap.add_argument("--min-dist", type=float, default=0.1)
    # Cosine, because that is how mean-pooled language-model embeddings are compared, and
    # what the crate's accuracy figures are quoted in.
    ap.add_argument("--metric", default="cosine")
    # Unset by default: umap-learn forces n_jobs=1 whenever random_state is given
    # (umap_.py:1950), which on a large corpus means one core instead of all of them.
    ap.add_argument("--seed", type=int, default=None,
                    help="fix the layout for reproducibility, at the cost of running "
                         "single-threaded")
    ap.add_argument("--n-jobs", type=int, default=-1,
                    help="-1 for every core. Ignored when --seed is given")
    ap.add_argument("--low-memory", action="store_true",
                    help="lower peak RAM in the kNN step, at some speed cost")
    args = ap.parse_args()

    src = args.work_dir / "embeddings.f32" # this is cool I didn't know it!
    if not src.exists():
        raise SystemExit(f"{src} not found; run `train_umap embed` first")

    x = np.fromfile(src, dtype=np.float32)
    if x.size % HIDDEN:
        raise SystemExit(f"{src} holds {x.size} floats, not a multiple of {HIDDEN}")
    x = x.reshape(-1, HIDDEN)
    print(f"{x.shape[0]} proteins x {x.shape[1]} dimensions")

    started = time.monotonic()
    coords = umap.UMAP(
        n_components=2,
        n_neighbors=min(args.n_neighbours, x.shape[0] - 1),
        min_dist=args.min_dist,
        metric=args.metric,
        random_state=args.seed,
        n_jobs=args.n_jobs,
        low_memory=args.low_memory,
    ).fit_transform(x)
    mode = "seeded, single-threaded" if args.seed is not None else f"unseeded, n_jobs={args.n_jobs}"
    print(f"UMAP in {time.monotonic() - started:.1f} s ({mode})")

    dest = args.work_dir / "coords.f32"
    np.ascontiguousarray(coords, dtype=np.float32).tofile(dest)
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
