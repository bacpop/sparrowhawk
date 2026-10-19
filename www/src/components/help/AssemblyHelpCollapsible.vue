<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Assembly</span>
      </div>
      <ChevronDown
        class="w-4 h-4 text-gray-500 transition-transform duration-200"
        :class="{ 'rotate-180': isOpen }"
      />
    </CollapsibleTrigger>

    <CollapsibleContent class="px-3 pb-3">
      <Tabs default-value="overview" class="mt-2">
        <TabsList class="grid w-full grid-cols-2 h-8">
          <TabsTrigger value="overview" class="text-xs">Overview</TabsTrigger>
          <TabsTrigger value="parameters" class="text-xs">Parameters</TabsTrigger>
        </TabsList>

        <div class="mt-3 text-sm text-gray-600 max-h-64 overflow-y-auto">
          <TabsContent value="overview" class="space-y-4">
            <p>
              Here you can perform genome assembly of short read sequencing data. You will need your sequencing read files (in
              FASTQ format, compressed or not), that can be single files or paired-end ones (keep in mind that these tools were
              designed for bacterial genomes, so large reads files might not work!). This tab uses <a href="https://github.com/bacpop/sparrowhawk-asm" target="_blank" class="text-blue-600 hover:underline">sparrowhawk-asm</a> under the hood.
            </p>

            <p>
              You can then "upload" them (though actually
              you are just loading them into your computer's memory), starting the assembly process. At the end, you can
              "download" (from memory) the final contigs as gzipped FASTA, as well as the intermediate de Bruijn graph in
              different gzipped formats. You can check the documentation for more information on what the different parameters do.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <p class="text-sm mb-1">
              <a href="https://www.ebi.ac.uk/ena/browser/view/ERR8158023" target="_blank" class="text-blue-600 hover:underline"><i>Mycobacterium tuberculosis</i> reads from ENA</a>
              </p>
              <ul class="list-disc list-inside space-y-1">
                <li>
                  <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_1.fastq.gz" class="text-blue-600 hover:underline">ERR8158023_1.fastq.gz</a>
                </li>
                <li>
                  <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR815/003/ERR8158023/ERR8158023_2.fastq.gz" class="text-blue-600 hover:underline">ERR8158023_2.fastq.gz</a>
                </li>
              </ul>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">k</dt>
                <dd class="ml-4">Most important factor, controls the size of subsequences used in the algorithm. It is an odd integer between 21 and 255. By default, 31.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Multi-k</dt>
                <dd class="ml-4">Assembles at several k in turn, each seeded with the contigs of the previous one, instead of at the single k above. Leave the k values empty for the default ladder (21, 31, 41, 55, 71, and up to the read length), or give odd values in ascending order, separated by commas. Slower, and usually more contiguous. Off by default.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Automatic read quality filtering / Minimum read quality</dt>
                <dd class="ml-4">On by default: at each k, the assembler chooses the minimum base quality a k-mer must have from the qualities in the reads, at most 20. Disable it to set the minimum read quality yourself, an integer from 0 to 33 (by default 20), used as is.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Automatic minimum count setting / Min counts for k-mer filtering</dt>
                <dd class="ml-4">On by default: the minimum count a k-mer needs is fitted from the k-mer spectrum, separately for each k. Disable it to set the minimum count yourself (by default 5). Memory usage could increase.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Use Bloom filter for preprocessing</dt>
                <dd class="ml-4">Counts k-mers with much less memory, at the cost of rare false-positive k-mers. On by default, as in the command-line tool.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Minimum contig length</dt>
                <dd class="ml-4">Contigs shorter than this many nucleotides are left out of the FASTA output; contigs of exactly this length are kept. By default, 500.</dd>
              </div>
            </dl>
          </TabsContent>
        </div>
      </Tabs>
    </CollapsibleContent>
  </Collapsible>
</template>

<script lang="ts">
import { defineComponent, ref } from "vue";
import { HelpCircle, ChevronDown } from "@lucide/vue";
import {
  Collapsible,
  CollapsibleContent,
  CollapsibleTrigger,
} from "@/components/ui/collapsible";
import {
  Tabs,
  TabsContent,
  TabsList,
  TabsTrigger,
} from "@/components/ui/tabs";

export default defineComponent({
  name: "AssemblyHelpCollapsible",
  components: {
    HelpCircle,
    ChevronDown,
    Collapsible,
    CollapsibleContent,
    CollapsibleTrigger,
    Tabs,
    TabsContent,
    TabsList,
    TabsTrigger,
  },
  setup() {
    const isOpen = ref(false);
    return { isOpen };
  },
});
</script>
