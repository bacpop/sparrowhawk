<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Mapping</span>
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
              Here you can map a sequence (reads or assemblies) or more to a reference one, normally a reference assembly. This method only works if all genomes are closely related (at the range of within-species distances). You
              will need first your reference sequence (in FASTA format, can be compressed), that you can "upload" (load in
              memory) and afterwards you will need to do the same with the sequence(s) you want to map. This uses <a href="https://github.com/bacpop/ska.rust" target="_blank" class="text-blue-600 hover:underline">ska.rust</a> under the hood.
            </p>
            <p>
              At the end you can
              download the mapping alignment as a FASTA file and also see it with an adapted <a href="https://github.com/lucanest/mseaboard" target="_blank" class="text-blue-600 hover:underline">MSEABOARD viewer</a> below. You can check the
              documentation for more information on what the different parameters do.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <div class="space-y-2">
                <div>
                    <p class="text-sm font-medium">Reference (<i>Mycobacterium tuberculosis</i> reference):</p>
                  <p class="text-sm">
                    <a href="https://www.ebi.ac.uk/ena/browser/view/GCA_000195955.2" target="_blank" class="text-blue-600 hover:underline">Record: GCA_000195955.2</a>
                  </p>
                  <ul class="list-disc list-inside">
                    <li>
                        <a href="https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000195955.2?download=true&gzip=true" class="text-blue-600 hover:underline">GCA_000195955.2.fasta.gz</a>
                    </li>
                  </ul>
                </div>
                <div>
                    <p class="text-sm font-medium">Files to map (<i>Mycobacterium tuberculosis</i> reads):</p>
                  <p class="text-sm">
                    <a href="https://www.ebi.ac.uk/ena/browser/view/ERR8158023" target="_blank" class="text-blue-600 hover:underline">Record: ERR8158023</a>
                  </p>
                  <ul class="list-disc list-inside">
                    <li>
                      <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/027/SRR30941327/SRR30941327_1.fastq.gz" class="text-blue-600 hover:underline">SRR30941327_1.fastq.gz</a>
                    </li>
                    <li>
                      <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/027/SRR30941327/SRR30941327_2.fastq.gz" class="text-blue-600 hover:underline">SRR30941327_2.fastq.gz</a>
                    </li>
                  </ul>
                </div>
              </div>
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
              <div>
                <dt class="font-medium text-gray-900">Mask ambiguous bases</dt>
                <dd class="ml-4">Replaces ambiguous bases in the alignment by "N".</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Mask repeats</dt>
                <dd class="ml-4">Replaces repeats in the alignment by "N".</dd>
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
  name: "MappingHelpCollapsible",
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
