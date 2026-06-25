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
import {downloadTextFile} from "@/downloadUtils";

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
    downloadFASTA(): void {
      downloadTextFile(this.allResults.fastaOutput, "assembly.fasta", "text/x-fasta");
    },
    downloadDOT(): void {
      downloadTextFile(this.allResults.dotOutput, "graph.dot", "text/vnd.graphviz");
    },
    downloadGFA(): void {
      downloadTextFile(this.allResults.gfaOutput, "graph.gfa");
    },
    downloadGFAv2(): void {
      downloadTextFile(this.allResults.gfav2Output, "graph.gfa2");
    },
  }
});
</script>

<style scoped>
#Download {
  float: left;
}
</style>
