<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Gene calling</span>
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
              Here you can call genes from an uploaded FASTA file using the <a href="https://github.com/hyattpd/Prodigal" target="_blank" class="text-blue-600 hover:underline">Prodigal</a> algorithm through <a href="https://github.com/fullhuman/orphos" target="_blank" class="text-blue-600 hover:underline">Orphos</a> (a port to Rust). Afterwards, you can download the output GFF file, as well as seeing your annotated genome in <a href="https://github.com/EBI-Metagenomics/mgnify-jbrowse" target="_blank" class="text-blue-600 hover:underline">an interactive browser</a>. Just set your preferred settings, and then drag-and-drop or click and select your genome(s) as FASTA, and the gene calling algorithm will automatically start. You can check the documentation for more information on how this works.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <p class="text-sm">
                <a href="https://www.ebi.ac.uk/ena/browser/view/GCA_004138665.1" target="_blank" class="text-blue-600 hover:underline"><i>Klebsiella pneumoniae</i> assembly: GCA_004138665.1</a>
              </p>
              <ul class="list-disc list-inside mt-1">
                <li>
                  <a href="https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true" class="text-blue-600 hover:underline">GCA_004138665.1.fasta.gz</a>
                </li>
              </ul>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">Workers</dt>
                <dd class="ml-4">When uploading multiple files, they are processed in different threads. This allows you to select how many you want to create. More will speed up processing of multiple files, but also consume more memory.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Translation table</dt>
                <dd class="ml-4">Here you can choose a different codon translation table or genetic code from the default one. All <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes" class="text-blue-600 hover:underline">NCBI tables</a> supported natively by Prodigal are available.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Use metagenomic mode</dt>
                <dd class="ml-4">This allows to activate Prodigal metagenomic mode, ideal for nonrelated, and/or fragmented sequences.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Ignore truncated genes</dt>
                <dd class="ml-4">By default the algorithm will try to see if a gene that started before the end of a contig could be continued in another contig. Activating this checkbox will disable that behaviour.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Break calling on N subsequences</dt>
                <dd class="ml-4">By default the algorithm will try to see if it is possible to call a gene with one or multiple uncertain nucleotides (usually represented by <code>N</code>). Activating this checkbox will disable that behaviour.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Ignore Shine-Dalgarno sequences</dt>
                <dd class="ml-4">If activated, the algorithm will not try to search and use these sequences to call genes.</dd>
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
  name: "GeneCallingHelpCollapsible",
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
