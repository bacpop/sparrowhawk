# sparrowhawk <img src='sparrowhawk_logo.png' align="right" height="250" />
Access Sparrowhawk at https://sparrowhawk.bacpop.org

Sparrowhawk is a WebAssembly toolkit for analysis of bacterial genomes including:

- Genome assembly with [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm).
- Mapping with [ska.rust](https://github.com/bacpop/ska.rust) and
  visualisation with [MSEABOARD](https://github.com/lucanest/mseaboard).
- Reference-free alignment with [ska.rust](https://github.com/bacpop/ska.rust) and phylogeny with [taxonium](https://github.com/theosanderson/taxonium).
- Single-linkage transmission networks with [ska.rust](https://github.com/bacpop/ska.rust) and [transmission_estimator](https://github.com/klebgenomics/transmission_estimator).
- Clustering and embedding with [mandrake.rust](https://github.com/bacpop/mandrake.rust).
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
## I want to create something similar to this with my method, what do I do?
We have created a [guide](https://sparrowhawk-docs.bacpop.org/wasm_guide) that details an example on how would such a method be created. Additionally, we have prepared a `SKILLS.md` and a `AGENTS.md` in the case you want to draft or outline such a project with the help of LLMs. Point them to the `SKILLS.md` file, that is designed to guide them developing both the interface of your method to the web part, and the web itself. We have improved this file by doing two exercises manually. We have adapted the nice bioinformatics read mapping crate [rammap-rs](https://github.com/jwanglab/rammap) to an Sparrowhawk analogous web in [rammap-web](https://github.com/vrbouza/rammap-web). Additionally, and to show the applicability on other scientific fields, we have done a small port of the [Combine](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit) toolkit, in C++ and Python, to Rust, and taking advantage of the Rust port of Minuit2 of [minuit2-rs](https://github.com/ricardofrantz/minuit2-rs), and a slightly extended [oxyroot-rs](https://github.com/m-dupont/oxyroot), we created [combine-web](https://github.com/vrbouza/combine-web) using Sparrowhawk as scaffold.



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
