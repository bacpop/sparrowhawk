<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Embed proteins</span>
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
              Here you can turn protein sequences into <b>embeddings</b>: fixed-length numerical vectors that capture
              something of a protein's biochemistry and likely structure, learnt by a protein language model. Proteins
              with similar embeddings tend to be functionally or structurally related, even when their sequences have
              diverged past the point where alignment is informative.
            </p>
            <p>
              Load a protein FASTA and each sequence is passed through
              <a href="https://huggingface.co/facebook/esm2_t6_8M_UR50D" target="_blank" class="text-blue-600 hover:underline">ESM-2 (8M)</a>,
              the smallest of the ESM-2 family, running entirely in your browser. You get one 320-dimensional vector per
              protein, averaged over its residues, which you can download as a TSV.
            </p>
            <p>
              The first use downloads the model (about 14 MB). Everything run on your CPU by default, which might be the faster choice in some occasions. You can enable optional GPU acceleration through WebAssembly
            </p>

            <p>
              Each protein is then placed on a two-dimensional plot, coloured by length, using a
              small network trained offline to reproduce a <b>UMAP</b> layout. Read it only as
              "these proteins are near those": the distance <i>between</i> two groups, and how
              tight or large a group looks, carry no meaning at all.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <p class="text-sm">
                <a href="https://www.uniprot.org/proteomes/UP000000625" target="_blank" class="text-blue-600 hover:underline"><i>Escherichia coli</i> K-12 reference proteome: UP000000625</a>
              </p>
              <ul class="list-disc list-inside mt-1">
                <li>
                  <a href="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz" class="text-blue-600 hover:underline">UP000000625_83333.fasta.gz</a>
                </li>
              </ul>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">GPU acceleration</dt>
                <dd class="ml-4">
                  Optional, and off by default. Whether WebGPU actually helps depends on the browser, operating system, and hardware.
                </dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Automatic device selection</dt>
                <dd class="ml-4">
                  On by default, and only relevant with GPU acceleration enabled. A page cannot enumerate the GPUs in
                  your machine, but it can ask for a high-performance or a low-power one, and the browser honours that
                  hint. Uncheck it to pick between them yourself. Those two are the only choices available, however many GPUs you have.
                </dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">GPU tasks per submission</dt>
                <dd class="ml-4">
                  How many pieces of GPU work are queued up into a single submission, and only relevant
                  with GPU acceleration enabled. The best value depends on your browser, which is why it
                  is adjustable.
                </dd>
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
  name: "ProteinEmbeddingsHelpCollapsible",
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
