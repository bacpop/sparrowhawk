# WASM checks

One Node 24 script per tool, each loading the site's WASM package (a web-target build) and replaying its
worker's calls at the page's defaults. `browser.mjs` stands in for `File`, `FileReaderSync` and `postMessage`, and the inputs are in the repository, mostly in
`data/` (see `data/README.txt`, and `make_examples.py` to rebuild it from the site's examples).

From `www/`:

    npm run test:wasm                                  # build every package, run every check
    npm run test:wasm -- --only sparrowhawk --no-build # one check, reusing the last packages
    npm run test:wasm -- --update-golden               # record new golden values after an intended change

Deterministic outputs are compared with `golden/<tool>-<target>.json`. A missing golden file is written on
the first run, to be reviewed and committed. They check correctness, not performance.

Not covered: host depletion (deacon-bridge); esm-bridge's WebGPU path, as Node has none; and mandrake's
embedding coordinates, which are seeded from the clock, so only its counts and shapes are checked.

The assembler's check also takes `--native-single` and `--native-multik`, contigs from
`sparrowhawk-asm rust/sparrowhawk-asm/tests/data/ecoli_window_{1,2}.fastq.gz --threads 2
--no-graphs --no-histo` (with `-k 31` for the first), and reports whether the browser builds match them.

The assembler check decompresses its four gzip downloads. It compares the FASTA's contigs with the
goldens; for single-k, it checks only that the DOT and GFA files have a header and are not empty. It also
checks that the output getters hand each buffer over only once, and covers fixed quality floors and empty
FASTA output.
