<template>
  <div class="flex flex-col gap-2">
    <div class="flex flex-col gap-2">
      Download Assembly
      <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadFASTA">
        .fasta
      </Button>
    </div>
    <div class="flex flex-col gap-2">
      Download De Bruijn graph
      <div class="flex flex-row gap-2">
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadDOT">
          .dot
        </Button>
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadGFA">
          .gfa
        </Button>
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadGFAv2">
          .gfa2
        </Button>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {defineComponent, Ref} from "vue";
import {useState} from "vuex-composition-helpers";
import type {AllResults} from "@/types";
import {Button} from "@/components/ui/button";
import {saveTextFile} from "@/platform/files";

export default defineComponent({
  name: 'DownloadButton',
  setup() {
    const {allResults} = useState(["allResults"]) as { allResults: Ref<AllResults> };

    return {
      allResults
    }
  },
  components: {
    Button
  },
  methods: {
    async downloadFASTA(): Promise<void> {
      await saveTextFile(this.allResults.fastaOutput, "assembly.fasta", "text/x-fasta");
    },
    async downloadDOT(): Promise<void> {
      await saveTextFile(this.allResults.dotOutput, "graph.dot", "text/vnd.graphviz");
    },
    async downloadGFA(): Promise<void> {
      await saveTextFile(this.allResults.gfaOutput, "graph.gfa");
    },
    async downloadGFAv2(): Promise<void> {
      await saveTextFile(this.allResults.gfav2Output, "graph.gfa2");
    },
  }
});
</script>

<style scoped>
#Download {
  float: left;
}
</style>
