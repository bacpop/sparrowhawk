# Gene calling

Gene calling predicts protein-coding genes (open reading frames) in bacterial genome sequences. Sparrowhawk offers this through [Orphos](https://github.com/FullHuman/orphos), a Rust port of [Prodigal](https://github.com/hyattpd/Prodigal) ([paper](https://doi.org/10.1186/1471-2105-11-119)), a widely-used prokaryotic gene predictor.

## General considerations

- Accepts FASTA files (plain or gzipped). Multiple files can be processed in parallel across workers.
- Two operating modes:
  - **Single-genome mode** (default, `metag=false`): all contigs in a file are first concatenated (separated by stop-codon bridges) to train the gene model on the whole genome, then each contig is analysed independently with that trained model. Best for complete or near-complete assemblies.
  - **Metagenomic mode** (`metag=true`), also called anonymous: each contig is analysed independently with a pre-trained general model. Recommended for metagenome-assembled genomes (MAGs) or any file with many short, unrelated sequences.
- Output files (downloadable after analysis):
  - GFF file with called genes.

## Parameters
The values in brackets are the default ones:
- **Workers** \[4\]: integer 1–8; number of parallel web workers. Multiple input files are distributed across workers.
- **Translation table** \[Default/Auto\]: NCBI translation table (1–25, or 0 for auto-detection). See the [official NCBI list](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for more information.
- **Use metagenomic mode** \[false\]: see General considerations above.
- **Ignore truncated genes** \[false\]: suppress reporting of genes that run off the edge of a contig (open ends).
- **Break calling on N subsequences** \[false\]: do not bridge over runs of unknown (`N`) bases when predicting genes.
- **Ignore Shine-Dalgarno sequences** \[false\]: force the algorithm to not use Shine-Dalgarno ribosome-binding site signals in gene calling.

## Example

The following data can be used to try out gene calling:

- Species: *Klebsiella pneumoniae*
- Assembly: [GCA_004138665.1](https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true) (FASTA, gzipped) from ENA
