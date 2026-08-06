# Introduction

**Sparrowhawk** is a web platform for doing basic bioinformatic analyses on bacterial genomes locally. It offers genome assembly, sequence mapping and alignment, taxonomic identification, gene calling, and host depletion.

You can see the different sections of this guide in the left, where each of the tabs/methodologies usage is detailed.


## Why use Sparrowhawk?
All analyses are run in your own system, with nothing being uploaded anywhere, nor a job running in a random foreign computing cluster: this makes Sparrowhawk ideal for resource-constrained setups (e.g. bad Internet connection), or for sensitive data studies. The methods are provided as web pages, being thus intuitive to be used, as well as available when you are using a tablet (or even a smartphone!).


<!-- ## Citation -->


## Code repositories
All code used is open source and most of it has been developed inside the [bacterial population genomics group (bacpop)](https://www.bacpop.org) at [EMBL's European Bioinformatics Institute (EMBL-EBI)](https://www.ebi.ac.uk). The main code repositories are the following:

- **Assembly**: [github.com/bacpop/sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm).
- **Mapping and alignment**: the software used is ska.rust, whose main repository is [github.com/bacpop/ska.rust](https://github.com/bacpop/ska.rust).
- **Taxonomic identification**: the software used is sketchlib.rust, whose main repository is [github.com/bacpop/sketchlib.rust](https://github.com/bacpop/sketchlib.rust).
- **Gene calling**: the software used is Orphos, a Rust port of [Prodigal](https://github.com/hyattpd/Prodigal), whose main repository is [github.com/FullHuman/orphos](https://github.com/FullHuman/orphos).
- **Host depletion**: the software used is Deacon, whose main repository is [github.com/bede/deacon](https://github.com/bede/deacon) ([preprint](https://doi.org/10.1101/2025.06.09.658732)).
- **Website**: [github.com/bacpop/sparrowhawk](https://github.com/bacpop/sparrowhawk).

Developers who want to adapt their own tools to run in the browser should consult the [Developer guide](./wasm_guide.md).
