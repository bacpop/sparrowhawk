<template>
  <div class="flex flex-col gap-2">
    <div class="flex flex-col gap-2">
      Download assembly
      <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadFASTA">
        .fasta.gz
      </Button>
    </div>
    <div class="flex flex-col gap-2">
      Download de Bruijn graph
      <div class="flex flex-row gap-2">
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadDOT">
          .dot.gz
        </Button>
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadGFA">
          .gfa.gz
        </Button>
        <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadGFAv2">
          .gfa2.gz
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
import {saveBinaryFile} from "@/platform/files";

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
      if (this.allResults.fastaGz) await saveBinaryFile(this.allResults.fastaGz, "assembly.fasta.gz", "application/gzip");
    },
    async downloadDOT(): Promise<void> {
      if (this.allResults.dotGz) await saveBinaryFile(this.allResults.dotGz, "graph.dot.gz", "application/gzip");
    },
    async downloadGFA(): Promise<void> {
      if (this.allResults.gfaGz) await saveBinaryFile(this.allResults.gfaGz, "graph.gfa.gz", "application/gzip");
    },
    async downloadGFAv2(): Promise<void> {
      if (this.allResults.gfa2Gz) await saveBinaryFile(this.allResults.gfa2Gz, "graph.gfa2.gz", "application/gzip");
    },
  }
});
</script>

<style scoped>
#Download {
  float: left;
}
</style>
