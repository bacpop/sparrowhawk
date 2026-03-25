# Taxonomic identification
This tab will try to assess what bacterial species correspond to any reads or genomic sequences that you "upload" (load into memory really) in it. Sparrowhawk offers this service through [sketchlib.rust](https://github.com/bacpop/sketchlib.rust): a further developed Rust port of [pp-sketchlib](https://github.com/bacpop/pp-sketchlib), that allows for fast comparisons of sequences. The data for the reference database comes from the [AllTheBacteria dataset](https://allthebacteria.org) ([preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059)), that comprises close to 2.5 million bacterial genomes uniformly assembled. These have been clustered using gemsparcl, to obtain high-quality groups of sequences whose species we know.

## General considerations
- The underlying code uses an inverted index (already implemented inside sketchlib.rust), built with sketches (akin to "fingerprints")[^notesketch] of representatives of high-quality clusters extracted from the AllTheBacteria dataset. When one sample is uploaded to this web page (either as a FASTA or FASTQ file), it is sketched as well, and its sketch is compared with the representative ones using sketchlib.rust's inverted index functionality. The sketches with highest [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) are returned.
- Multiple files are processed simultaneously, one per worker, using parallel Web Workers in the browser.

## Parameters
The values in brackets are the default ones:
- **Min Illumina read quality** \[20\]: integer 0–33; filters nucleotides below this quality score when processing FASTQ reads.
- **Min counts for k-mer filtering** \[5\]: integer 1–30; only k-mers appearing more than this threshold are used. Analogous to the assembly parameter of the same name.
- **Proportion of reads** \[1.0\]: a real number between 0 and 1 that allows you to choose what fraction of reads to sample.
- **Workers** \[4\]: integer 1–8; number of parallel web workers processing files simultaneously. Higher values speed up batch identification but use more memory.

## Example

The following data can be used to try out taxonomic identification:

- Species: *Klebsiella pneumoniae*
- Assembly: [GCA_004138665.1](https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true) (FASTA, gzipped) from ENA


[^notesketch]: See the [sketchlib.rust code repository](https://github.com/bacpop/sketchlib.rust) for more documentation on how sketches are obtained from a genome.
