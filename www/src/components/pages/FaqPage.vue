<template>
  <div class="w-1/3 max-w-[350px]">
    <h1 class="text-2xl font-medium mb-4">
      FAQ
    </h1>
  </div>
  <div class="max-w-3/4">
    <h3 class="mt-2 mb-2 font-medium">
      Where is the full documentation?
    </h3>
    <p class="text-sm">
        <a href="https://sparrowhawk-docs.bacpop.org/">https://sparrowhawk-docs.bacpop.org/</a>
    </p>
    <h3 class="mt-2 mb-2 font-medium">
      Who is this website for?
    </h3>
    <p class="text-sm">
      Anyone with whole genome sequence data of bacteria who wishes to do genome assembly, sequence mapping and
      alignment, and taxonomic identification from their Internet browser, without installing anything in their computer
      nor using command line software. It is also useful for sensitive data analysis, as everything is run locally and
      nothing leaves your computer. This also makes it applicable for low resource contexts, such as having a bad
      Internet connection, or even not having a laptop at your disposal, as this website can be run in your
      smartphone/tablet browser!
    </p>

    <h3 class="mt-2 mb-2 font-medium">
      What can I do in the "Assembly" tab?
    </h3>
    <p class="text-sm">
      Here you can perform genome assembly of short read sequencing data. You will need your sequencing read files (in
      FASTQ format, compressed or not), that can be single files or paired-end ones (keep in mind that these tools were
      designed for bacterial genomes, so large reads files might not work!). You can then "upload" them (though actually
      you are just loading them into your computer's memory), starting the assembly process. At the end, you can
      "download" (from memory) the final contigs in FASTA format, as well as the intermediate de Bruijn graph in
      different formats. You can check the documentation for more information on what the different parameters do.
    </p>

    <h3 class="mt-2 mb-2 font-medium">
      What can I do in the "Mapping" tab?
    </h3>
    <p class="text-sm">
      Here you can map a sequence (reads or assemblies) or more to a reference one, normally a reference assembly. This algorithm only works if the genomes are closely related (within-species distances). You
      will need first your reference sequence (in FASTA format, can be compressed), that you can "upload" (load in
      memory) and afterwards you will need to do the same with the sequence(s) you want to map. At the end you can
      download the mapping alignment as a FASTA file and also see it with a viewer below. You can check the
      documentation for more information on what the different parameters do.
    </p>

    <h3 class="mt-2 mb-2 font-medium">What can I do in the "Alignment" tab?</h3>
    <p class="text-sm">
      Here you can align three or more sequences (reads or assemblies) between themselves. This algorithm only works if the genomes are closely related (within-species distances). After "uploading" (loading
      into memory) them, and after they are processed, you will see a quick phylogenetic tree built from their
      differences. Take into account that the underlying software, ska.rust, is designed to work with small taxonomic
      distances (inside a species), so inputting sequences genomically very different might make the web show absurd
      phylogenetic trees. You can check the documentation for more information on what the different parameters do and
      how this works.
    </p>

    <h3 class="mt-2 mb-2 font-medium">What can I do in the "Taxonomic ID" tab?</h3>
    <p class="text-sm">
      Here you can identify the taxonomic classification at the level of species of bacterial genomic sequences (either
      reads in FASTQ format, or any FASTA sequence). For that, you just need to "upload" (load into memory) your
      sequence files, and you will obtain the three most closest results, with its probability, that the methodology
      obtains for your sequences. You can download the results in TSV format. For analysing several files, you can increase
      the number of threads as multi-file processing is parallelised. You can check the documentation for more information on how this works.
    </p>

    <h3 class="mt-2 mb-2 font-medium">What can I do in the "Gene calling" tab?</h3>
    <p class="text-sm">
      Here you can call genes from any bacterial genomic sequences in FASTA format. For that, you just need to "upload" (load into memory) your sequence files, and you will obtain a GFF file to download. We have also included a genome viewer so that you 
      will be able to see your recently called genes. You can analyse simultaneously multiple files, and download in bulk the results.
      This algorithm is parallelised when you analyse more than one file, you can set the amount of parallel threads in the arguments. You can check the documentation for more information on how this works.
    </p>

    <h3 class="mt-2 mb-2 font-medium">What can I do in the "Host depletion" tab?</h3>
    <p class="text-sm">
      Here you can call either filter or select reads or sequences in general. This is particularly useful when dealing with human bacterial samples, that might be contaminated with host DNA, and you usually want to filter as much as you can that part, for data protection and privacy reasons.
      To use this, first upload an index. There is one example provided in the "How to use Host depletion" block that is a compact human read index. You can create your own using Deacon. Afterwards, you can upload your reads (in pairs or individual files) by drag and drop or selecting them manually. You can process multiple sequentially (this is not parallelised). At the end, all the results will be available for download individually or in bulk as .zip or .tar.gz files. You can check the documentation for more information on how this works.
    </p>

    <h3 class="mt-2 mb-2 font-medium">What can I do in the "Embed proteins" tab?</h3>
    <p class="text-sm">
      Here you can turn protein sequences into embeddings: fixed-length numerical vectors produced by a protein language model, which capture something of a protein's biochemistry and likely structure. Proteins that are functionally or structurally related tend to have similar embeddings even when their sequences have diverged too far for alignment to help, which makes them a useful starting point for clustering, similarity search, or visualisation.
      Just "upload" (load into memory) a protein FASTA file and each sequence is passed through ESM-2 (8M), the smallest ESM-2 model, running entirely in your browser. You get one 320-dimensional vector per protein, averaged over its residues, downloadable as a TSV. The first use downloads the model (about 14 MB), which your browser then caches. If your browser and hardware support WebGPU the model runs on your GPU, which is much faster than the CPU fallback. You can check the documentation for more information on how this works.
    </p>

    <h3 class="mt-2 mb-2 font-medium">Where can I find information on what parameters to use?</h3>
    <p class="text-sm">
      You can hover over the parameters and get a small description of each one.
    </p>

    <h3 class="mt-2 mb-2 font-medium">
      How do I add tools or develop my own WASM applications with Sparrowhawk?
    </h3>
    <p class="text-sm">
        More information is at <a href="https://sparrowhawk-docs.bacpop.org/wasm_guide/">https://sparrowhawk-docs.bacpop.org/wasm_guide</a>
    </p>

    <h3 class="mt-2 mb-2 font-medium">Where can I find the code and references for each method?</h3>
    <p class="text-sm">
      All the code is fully open-source and you can access it on the following public repositories:
    </p>
    <ul class="">
      <li class="text-sm list-disc ml-5">
        Underlying assembly code: sparrowhawk-asm. <a target="_blank" href="https://github.com/bacpop/sparrowhawk-asm" class="text-blue-500 underline">Link to repository</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Underlying mapping and alignment code: ska.rust. <a target="_blank" href="https://github.com/bacpop/ska.rust" class="text-blue-500 underline">Link to repository</a>. 
        <a target="_blank" href="https://genome.cshlp.org/content/34/10/1661.abstract" class="text-blue-500 underline">Link to paper</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        MSA viewer: MSEABOARD. <a target="_blank" href="https://github.com/lucanest/mseaboard"
                                 class="text-blue-500 underline">Link to repository</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Phylogenetic tree viewer: Taxonium. <a target="_blank" href="https://github.com/theosanderson/taxonium"
                                 class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://elifesciences.org/articles/82392" class="text-blue-500 underline">Link to paper</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Underlying taxonomic identification code: sketchlib.rust. <a target="_blank" href="https://github.com/bacpop/sketchlib.rust"
                                                     class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://doi.org/10.1098/rstb.2021.0237" class="text-blue-500 underline">Link to related paper where the methodology is detailed</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Cluster information for taxonomic identification extracted with Gemsparcl. <a target="_blank" href="https://www.biorxiv.org/content/10.64898/2025.12.30.695181" class="text-blue-500 underline">Link to preprint</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Data use to create reference index for taxonomic identification: All The Bacteria (AtB) dataset. <a target="_blank" href="https://allthebacteria.org/" class="text-blue-500 underline">Link to website</a>. <a target="_blank" href="https://github.com/AllTheBacteria/AllTheBacteria" class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://www.biorxiv.org/content/10.1101/2024.03.08.584059" class="text-blue-500 underline">Link to preprint</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Underlying gene calling code: Orphos, a Prodigal Rust port. <a target="_blank" href="https://github.com/FullHuman/orphos" class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://github.com/hyattpd/Prodigal" class="text-blue-500 underline">Link to Prodigal's repository</a>. <a target="_blank" href="https://doi.org/10.1186/1471-2105-11-119" class="text-blue-500 underline">Link to Prodigal's paper</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Protein embedding model: ESM-2. <a target="_blank" href="https://huggingface.co/facebook/esm2_t6_8M_UR50D" class="text-blue-500 underline">Link to model</a>. <a target="_blank" href="https://github.com/facebookresearch/esm" class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://www.science.org/doi/10.1126/science.ade2574" class="text-blue-500 underline">Link to paper</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Deep-learning framework running the embedding model in the browser: Burn. <a target="_blank" href="https://github.com/tracel-ai/burn" class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://github.com/tracel-ai/burn-onnx" class="text-blue-500 underline">Link to its ONNX importer</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Underlying host depletion code: Deacon. <a target="_blank" href="https://github.com/bede/deacon" class="text-blue-500 underline">Link to repository</a>. <a target="_blank" href="https://doi.org/10.1101/2025.06.09.658732" class="text-blue-500 underline">Link to preprint</a>.
      </li>
      <li class="text-sm list-disc ml-5">
        Website: <a target="_blank" href="https://github.com/bacpop/sparrowhawk" class="text-blue-500 underline">link to repository</a>.
      </li>
    </ul>
  </div>
</template>
<script setup lang="ts">
</script>