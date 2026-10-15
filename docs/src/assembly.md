# Assembly

Genome assembly consists on recovering a genome (the DNA sequence) from organisms that have been sequenced. Sparrowhawk offers this service through a novel assembler, written in Rust, that uses as main structure a node-based [de Bruijn graph](https://en.wikipedia.org/wiki/De_Bruijn_graph).

## General considerations
- The assembler has been designed with bacterial genomes (\\( O(10^0-10^1)\\,\text{Mbp} \\)) on mind, so large read files might not be work properly with it.
- The aim for its development was for it to be efficient and able to be loaded into a website (using WebAssembly), so it is not the most thorough and detailed genomic assembler. Its performance is overall similar/superior to Velvet, and similar to that of Minia (with single k-value).
- The assembler runs as 64-bit WebAssembly[^memory64] in browsers that support it, which lets it use up to 16 GB of memory. Other browsers fall back to the 32-bit build, limited to 4 GB, and the page warns when this happens; there, high-coverage reads or those coming from large-genome species might reach that limit.
- If you want to run this assembler locally in the command line, you can find more information in [its GitHub repository](https://github.com/bacpop/sparrowhawk-asm).

## Parameters
The values in brackets are the default ones:
- **k** \[31\]: controls the size of subsequences used in the algorithm. It is an odd integer that should be between ~21 up until 255.
- **Multi-k** \[false\]: assembles at several k in turn, each seeded with the contigs of the previous one, instead of at the single k above. Leave the k values empty for the default ladder (21, 31, 41, 55, 71, and then up to the read length), or give odd values in ascending order, separated by commas. Slower, and usually more contiguous.
- **Automatic read quality filtering** \[true\]: at each k, the assembler chooses the minimum base quality a k-mer must have from the qualities in the reads, at most 20, and can lower it when the k-mer spectrum asks for it.
- **Minimum read quality** \[20\]: shown when the previous option is disabled. The minimum base quality, an integer from zero to 33, used as is at every k.
- **Automatic minimum count setting** \[true\]: the assembler fits the k-mer spectrum distribution and extracts the ideal minimum count, separately for each k. If using, it is recommended to check in the k-mer spectrum plot the chosen value.
- **Min counts for k-mer filtering** \[5\]: shown when the previous option is disabled. This value controls the cutoff for filtering low-count k-mers, that are usually associated to sequencing errors. Low values will imply larger memory usage. Its value is correlated with k. It is recommended to check in the k-mer spectrum histogram how much of it is filtered. It should not take the main bell part of the distribution, only the peak that should be near one or two.
- **Use Bloom filter for preprocessing** \[true\]: counts k-mers through a Bloom filter first, as the command-line tool does by default. This significantly reduces the memory used, while allowing for some degree of overcounting some k-mers.
- **Minimum contig length** \[500\]: contigs shorter than this many nucleotides are left out of the FASTA output; contigs of exactly this length are kept.


## Example

The following data can be used to try out the assembler:

- Species: *Mycobacterium tuberculosis*
- Reads (ENA accession ERR8158023): [forward](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_1.fastq.gz) and [reverse](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_2.fastq.gz)


[^memory64]: The memory64 proposal gives WebAssembly 64-bit memory addresses, removing the 4 GB limit of 32-bit ones. The assembler is compiled for both `wasm64-unknown-unknown` and `wasm32-unknown-unknown`, and the page loads the 64-bit build first. The 16 GB ceiling is the maximum memory the 64-bit build declares.