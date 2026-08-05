<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Alignment</span>
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
              Here you can align three or more sequences (reads or assemblies) between themselves. After "uploading" (loading
              into memory) them, and after they are processed, you will see a quick phylogenetic tree built from their
              differences using <a href="https://github.com/theosanderson/taxonium" target="_blank" class="text-blue-600 hover:underline">Taxonium</a>. This tab uses <a href="https://github.com/bacpop/ska.rust" target="_blank" class="text-blue-600 hover:underline">ska.rust</a> under the hood.
            </p>
            <p>
              Take into account that the underlying software, ska.rust, is designed to work with small taxonomic
              distances (inside a species), so inputting sequences genomically very different might make the web show absurd
              phylogenetic trees. You can check the documentation for more information on what the different parameters do and
              how this works.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <p class="text-sm mb-1">
                Ten (or more) E. coli assemblies from a study of a Tanzanian hospital:
              </p>
              <p class="text-sm">
                <a href="https://www.ebi.ac.uk/ena/browser/view/PRJEB23541" target="_blank" class="text-blue-600 hover:underline">Study: PRJEB23541</a>
              </p>
              <p class="text-xs text-gray-500 mt-1">
                To download: click on each "Accession" value in the table, then click "SET FASTA" on the right.
              </p>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">k</dt>
                <dd class="ml-4">Controls the size of the subsequences (k-mers) used to process the reads or assemblies. Similar to that of the assembly, but in this case only up until 63.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Minimum Illumina read quality</dt>
                <dd class="ml-4">The quality of the input fastq files, it is an integer from zero to 33. By default, 20.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Proportion of reads</dt>
                <dd class="ml-4">When using reads, this allows to subsample them (making thus the analysis faster). A real number between 0 and 1.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Quality filter type</dt>
                <dd class="ml-4">Allows to set the Illumina filter type for reads that the algorithm will do. By default, it will check all nucleotides in a k-mer. You can set it to check only the middle one, or noone.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Use canonical k-mers</dt>
                <dd class="ml-4">It allows to identify identically one k-mer and its reverse-complement. It is recommended when uploading reads.</dd>
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
  name: "AlignmentHelpCollapsible",
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
