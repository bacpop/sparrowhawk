# Mapping

Mapping is often used as an alternative way of reconstructing a genome from sequencing raw data, instead of trying to re-assemble it from scratch. Sparrowhawk offers this service through [ska.rust](https://github.com/bacpop/ska.rust) ([paper](https://genome.cshlp.org/content/34/10/1661.abstract)): a toolkit for prokaryotic (and other haploid and small genome beings) sequence analysis with split kmers.

## General considerations

After processing, the mapped sequences are displayed in an interactive [MSEABOARD](https://github.com/lucanest/mseaboard) sequence viewer directly in the browser.

## Parameters
The values in brackets are the default ones:
- **k** \[31\]: controls the size of subsequences used in the algorithm. It is an odd integer that should be between ~17 up until 63.
- **Proportion of reads** \[1.0\]: a real number between 0 and 1 that allows you to choose, when you upload reads, what fraction of those you want to sample (if you want to do so).
- **Min Illumina read quality** \[20\]: integer 0–33; filters nucleotides below this quality score.
- **Quality filter type** \[All bases\]: dropdown; "No filter" / "Middle base" (only the central split-kmer base) / "All bases" (any base in the k-mer).
- **Use canonical k-mers** \[false\]: accounts for both strand orientations; recommended when using raw reads. Can only be set before uploading any file.
- **Mask ambiguous bases** \[false\]: replaces ambiguous bases in the output with `N`. Can only be set before uploading any file.
- **Mask repeats** \[false\]: masks repeated regions in the output with `N`. Can only be set before uploading any file.

## Example

The following data can be used to try out the mapper:

- Species: *Mycobacterium tuberculosis*
- Reference genome: [GCA_000195955.2](https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000195955.2?download=true&gzip=true) (FASTA, gzipped)
- Reads to map (SRR30941327): [forward](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/027/SRR30941327/SRR30941327_1.fastq.gz) and [reverse](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/027/SRR30941327/SRR30941327_2.fastq.gz)
