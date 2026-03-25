# Alignment
Alignment is a way of showing the differences between two or more DNA sequences, often used in evolution/phylogenetic studies. Sparrowhawk offers this service through [ska.rust](https://github.com/bacpop/ska.rust) ([paper](https://genome.cshlp.org/content/34/10/1661.abstract)): a toolkit for prokaryotic (and other haploid and small genome beings) sequence analysis with split kmers.

## General considerations
- Take into account that ska.rust (the underlying methodology) is designed to work with closely related genomes (less than species level). If you use very different sequences, the results might not be accurate.
- You will need three or more sequences to perform the alignment.
- After processing, the alignment is displayed in an interactive multiple sequence alignment (MSA) viewer in the browser, with the phylogenetic tree shown alongside it.
- The phylogenetic tree is computed with the [speedytree crate](https://crates.io/crates/speedytree) and rendered with [Taxonium](https://github.com/theosanderson/taxonium) ([paper](https://doi.org/10.7554/eLife.82392)).

## Parameters
The values in brackets are the default ones:
- **k** \[31\]: controls the size of subsequences used in the algorithm. It is an odd integer that should be between ~17 up until 63.
- **Proportion of reads** \[1.0\]: a real number between 0 and 1 that allows you to choose, when you upload reads, what fraction of those you want to sample (if you want to do so).
- **Min Illumina read quality** \[20\]: integer 0–33; filters nucleotides below this quality score.
- **Quality filter type** \[All bases\]: dropdown; "No filter" / "Middle base" (only the central split-kmer base) / "All bases" (any base in the k-mer).
- **Use canonical k-mers** \[false\]: accounts for both strand orientations; recommended when using raw reads. Can only be set before uploading any file.

## Example

The following data can be used to try out the aligner:

- Species: *Escherichia coli* (three or more assemblies required)
- Dataset: ENA study [PRJEB23541](https://www.ebi.ac.uk/ena/browser/view/PRJEB23541) — Tanzanian hospital isolates. Click any "Accession" value in the study table, then click "SET FASTA" in the right panel to load each assembly.
