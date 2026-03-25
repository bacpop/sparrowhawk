# Host depletion

Host depletion filters sequencing reads to remove (or isolate) reads that match a host genome index, retaining the reads of interest for downstream analysis. Sparrowhawk uses [Deacon](https://github.com/bede/deacon) ([preprint](https://doi.org/10.1101/2025.06.09.658732)) for this purpose.

## General considerations

- Requires uploading a pre-built Deacon index (`.dci` file) for the host genome of interest. The Deacon repository contains pre-built indexes for common hosts (human, etc.).
- Accepts FASTQ reads (plain or gzipped); output is written in the same format (optionally gzip-compressed).
- The two thresholds (absolute and relative) are applied together: a read must satisfy both to be classified as host-derived.

## Parameters
The values in brackets are the default ones:
- **Deplete mode** \[enabled\]: when enabled, reads matching the host index are *removed* from the output. When disabled, matching reads are *kept* instead (useful for isolating host reads).
- **Absolute threshold** \[1\]: integer 1–50; minimum number of k-mers that must match the host index for a read to be classified as host-derived.
- **Relative threshold** \[0.05\]: real 0–1; minimum proportion of a read's k-mers that must match the host index.

## Example

The following data can be used to try out host depletion:

- Human host index (pre-built, panhuman-1, k=31 w=61): [panhuman-1.k31w61.idx](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w61.idx)
- Reads to deplete (as a test FASTA): *Klebsiella pneumoniae* assembly [GCA_004138665.1](https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true)
