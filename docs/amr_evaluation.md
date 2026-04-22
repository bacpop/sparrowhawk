# AMR Evaluation Workflow

The AMR bridge includes a native `x86_64` CLI so the same Rust detection core can be
evaluated against native ResFinder before shipping the reduced WebAssembly profile.

## Build the Rust index

```bash
cargo run --manifest-path rust/amr-bridge/Cargo.toml -- build-index \
  --db-root /home/vrbouza/Projects/assembler_development/amr_tests/resfinder_db \
  --out /tmp/resfinder.amridx \
  --k 31
```

## Build KMA

```bash
make -C /home/vrbouza/Projects/assembler_development/amr_tests/kma
```

## Run evaluation against native ResFinder

```bash
cargo run --manifest-path rust/amr-bridge/Cargo.toml -- eval \
  --index /tmp/resfinder.amridx \
  --fasta /home/vrbouza/Projects/assembler_development/amr_tests/GCF_020579515.1_ASM2057951v1_genomic.fna \
  --resfinder-root /home/vrbouza/Projects/assembler_development/amr_tests/resfinder \
  --db-root /home/vrbouza/Projects/assembler_development/amr_tests/resfinder_db \
  --out-dir /tmp/resfinder_eval \
  --blastn-path /usr/bin/blastn \
  --kma-path /home/vrbouza/Projects/assembler_development/amr_tests/kma/kma
```

If you have already produced a ResFinder JSON file, pass `--baseline-json` instead of re-running ResFinder.
