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
              designed for bacterial genomes, so large reads files might not work!). This tab uses <a href="https://github.com/bacpop/sparrowhawk" target="_blank" class="text-blue-600 hover:underline">sparrowhawk</a> under the hood.
            </p>

            <p>
              You can then "upload" them (though actually
              you are just loading them into your computer's memory), starting the assembly process. At the end, you can
              "download" (from memory) the final configs in FASTA format, as well as the intermediate de Bruijn graph in
              different formats. You can check the documentation for more information on what the different parameters do.
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
                <dd class="ml-4">Most important factor, controls the size of subsequences used in the algorithm. It is an odd integer that should be between ~21 up until 255 (currently until 89). By default, 31.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Minimum Illumina read quality</dt>
                <dd class="ml-4">The quality of the input fastq files, it is an integer from zero to 33. By default, 20.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Minimum counts for k-mer filtering</dt>
                <dd class="ml-4">Another integer for a filtering step, by default 5.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Automatically set the minimum counts for k-mer filtering</dt>
                <dd class="ml-4">This box does what it says, it automatically disables the slider for the minimum counts. Memory usage could increase.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Use Bloom filter for preprocessing</dt>
                <dd class="ml-4">Recommended for non-small read files. This box disables chunking automatically, removing the previous text box.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Keep dead-ends / Do not collapse bubbles</dt>
                <dd class="ml-4">The last two options are booleans that modify two parameters from the algorithm. They do not affect the rest of the parameters.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Chunk processing size (set zero for no chunking)</dt>
                <dd class="ml-4">A value to be put, by default 150000, anything inputted should be a positive integer (including zero).</dd>
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
