# Assembly

Genome assembly consists on recovering a genome (the DNA sequence) from organisms that have been sequenced. Sparrowhawk offers this service through a novel assembler, written in Rust, that uses as main structure a node-based [de Bruijn graph](https://en.wikipedia.org/wiki/De_Bruijn_graph).

## General considerations
- The assembler has been designed with bacterial genomes (\\( O(10^0-10^1)\\,\text{Mbp} \\)) on mind, so large read files might not be work properly with it.
- The aim for its development was for it to be efficient and able to be loaded into a website (using WebAssembly), so it is not the most thorough and detailed genomic assembler. Its performance is overall similar/superior to Velvet, and similar to that of Minia (with single k-value).
- Given a current technical limitation[^notememory64] in the amount of memory that can be used (4 GB), high coverage reads or those coming from large genome species might reach that limit.
- If you want to run this assembler locally in the command line, you can find more information in [its GitHub repository](https://github.com/bacpop/sparrowhawk-asm).

## Parameters
The values in brackets are the default ones:
- **k** \[31\]: controls the size of subsequences used in the algorithm. It is an odd integer that should be between ~21 up until 255.
- **Minimum Illumina read quality** \[20\]: the quality of the input fastq files, it is an integer from zero to 33.
- **Minimum counts for k-mer filtering** \[5\]: this value controls the cutoff for filtering low-count k-mers, that are usually associated to sequencing errors. Low values will imply larger memory usage. Its value is correlated with k. It is recommended to check in the k-mer spectrum histogram how much of it is filtered. It should not take the main bell part of the distribution, only the peak that should be near one or two.
- **Automatically set the minimum counts for k-mer filtering** \[false\]: the assembler has an automatic setting for fixing the previous parameter, that fits the k-mer spectrum distribution and extracts the ideal value. If using, it is recommended to check in the k-mer spectrum plot the chosen value.
- **Chunk processing size** \[150000\]: reads are by default processed in chunks, reducing the memory usage. Their size can be modified here, or disabled completely (by setting it to zero).
- **Use Bloom filter for preprocessing** \[false\]: if facing issues with memory usage, try using a Bloom filter for the preprocessing step. This significantly reduces the memory used, while allowing for some degree of overcounting some k-mers.
- **Do not remove dead-ends** \[false\]: this removes the dead-end removal step in the error correction algorithm of the assembler. 
- **Do not collapse bubbles** \[false\]: this removes the simple bubble removal step in the error correction algorithm of the assembler.


## Example

The following data can be used to try out the assembler:

- Species: *Mycobacterium tuberculosis*
- Reads (ENA accession ERR8158023): [forward](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_1.fastq.gz) and [reverse](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_2.fastq.gz)


[^notememory64]: The intermediate libraries for compiling Rust code into WebAssembly have not implemented the memory64 proposal, that allows to have 64 bit memory addresses (thus removing the 4GB memory limit). Once this is done, compilation to `wasm64-unknown-unknown` should be easily possible with them and thus we would remove this annoying limitation.