# sparrowhawk <img src='sparrowhawk_logo.png' align="right" height="250" />
Access Sparrowhawk at https://sparrowhawk.bacpop.org

Sparrowhawk is a WebAssembly toolkit for analysis of bacterial genomes including:

- Genome assembly with [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm).
- Mapping with [ska.rust](https://github.com/bacpop/ska.rust) and
  visualisation with [MSEABOARD](https://github.com/lucanest/mseaboard).
- Reference-free alignment with [ska.rust](https://github.com/bacpop/ska.rust) and phylogeny with [taxonium](https://github.com/theosanderson/taxonium).
- Single-linkage transmission networks with [ska.rust](https://github.com/bacpop/ska.rust) and [transmission_estimator](https://github.com/klebgenomics/transmission_estimator).
- Taxonomic identification (species-identification) with [sketchlib.rust](https://github.com/bacpop/sketchlib.rust) and [gemsparcl](https://github.com/johannahelene/gemsparcl).
- Annotation of genes with [prodigal](https://github.com/hyattpd/Prodigal) via [orphos](https://github.com/FullHuman/orphos).
- Host read depletion with [deacon](https://github.com/bede/deacon).
- AMR gene detection with a k-mer-ized [AMRFinderPlus](https://github.com/ncbi/amr) database.
- Protein embedding with [ESM-2 (8M)](https://biolm.ai/models/esm2-8m/).

## Documentation

Docs are here: https://sparrowhawk-docs.bacpop.org

## WebAssembly?

- Easy: WebAssembly (WASM) runs directly in your browser, so no installation is
  needed.
- Private: Everything is processed locally, on your own machine, so no data is
  transferred to external servers.
- Fast: All code is optimised for rapid analysis, and some modules even
  support using your local GPU for increased speed.
- Up-to-date: Our packages build both native and WASM code, so the
  latest features and updates continue to be included.

For larger analyses, some browsers (Chrome, Safari, Edge) are faster than others.

## Sparrowhawk?
Sparrowhawk was at one time the Archmage of [Earthsea](https://en.wikipedia.org/wiki/Earthsea).
Also, the [sparrowhawk](https://en.wikipedia.org/wiki/Eurasian_sparrowhawk) (*Accipiter nisus*) is a bird of prey native to Europe (and the island of Gont).

For fans and haters of backronyms: **S**peedy **P**rivate **AR**chitecture **R**unning **O**n-device **WH**ich **A**nalyses **W**ith **K**-mers

## Contributors

- [Víctor Rodríguez Bouza](https://github.com/vrbouza)
- [John Lees](https://github.com/johnlees)
- [Bede Constantinides](https://github.com/bede/)
- [Luca Nesterenko](https://github.com/lucanest/)
- [Erikson Odih](https://github.com/Erkison)
- [Andrea Epifani](https://github.com/andyepx)
- [Antoine Andreoletti](https://github.com/apollis44)

<br>
<br>
<br>
<br>

---

## Disclaimer :warning: :construction:
This is a **work in progress** project. This in particular implies:

- Not all the main features we want are yet implemented.
- Code might be messy, and not even documented.
- General documentation on how to install and use the tool might be short or even missing.
- Finding unexpected errors/behaviour or bugs should not be a surprise.
- Some features might be partially hardcoded.

These (and potentially other) items will be progressively fixed before version 1.0.

---

# Developers (developers, developers, developers)

**Note:** this repository is for the web interface of the assembler. If you are looking for the assembler itself (that can be run locally) see [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm).

This web interface aims to offer a way of having a simple website that offers the WebAssembly compiled [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm) assembler. It has been developed taking advantage/inspiration from other WebAssembly projects from our group (such as [DATACIN](https://github.com/bacpop/DATACIN)).

Current **main features** (see the [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm) repository for details on the assembler itself):
- Simple, working, web-interface.
- Allows to drag-and-drop (or select from a file browser) the Illumina paired-end reads from your computer.
- Download of the assembled contigs in FASTA format, as well as the de Bruijn graph before collapse as DOT, GFAv1.1, and GFAv2 formats.
- Customised parameter setting.
- Due to the current 32bit memory addresses restriction, there is a **total 4GB RAM limit**. This implies that some reads won't be able to be assembled using the standard preprocessing, due to their size. Some additional options have been implemented to try to circumvent this restriction, though keep in mind that they might still not be enough. These are:
    - Chunking of the preprocessing.
    - Alternative preprocessing using a Bloom filter.


## Installation
Development has been done only on x86_64 GNU/Linux-based systems, and most surely will probably stay that way (i.e. no other systems have been tested). To use it you will need to have [NPM](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) installed, as well as a working Rust installation (check the [official Rust website](https://rust-lang.org/tools/install/) for information on how to get it) and then download and install the required packages as follows (for example, to get the v0.1.2 version)

```
git clone --branch v0.1.2 --recurse-submodules https://github.com/bacpop/sparrowhawk.git
cd sparrowhawk/www
npm install
```
### Note for macOS enjoyers
It has been found that there is an issue when compiling some crates to wasm32-unknown-unknown in macOS due to the weird clang packaging there. Following [this issue](https://github.com/gyscos/zstd-rs/issues/93), it seems that installing llvm fixes it.

## Usage
Once you have it installed, you can run the following, that will automatically compile to WebAssembly [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm) and run the development server locally (from the www folder!):

```
npm run serve
```

## Customize configuration
See [Configuration Reference](https://cli.vuejs.org/config/).
