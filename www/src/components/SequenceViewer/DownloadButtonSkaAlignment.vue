
<template>
  <div class="flex flex-wrap gap-6 items-start">
    <div class="flex flex-col gap-2">
      <span>Download alignment</span>
      <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadALN">
        .aln
      </Button>
    </div>

    <div class="flex flex-col gap-2">
      <span>Download tree</span>
      <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadTREE">
        .tree
      </Button>
    </div>

    <div class="flex flex-col gap-2">
      <span>Download distances</span>
      <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadCSV">
        .csv
      </Button>
    </div>
  </div>
</template>


<script lang="ts">
import {defineComponent} from "vue";
import {useState} from "vuex-composition-helpers";
import {Button} from "@/components/ui/button";
import {downloadTextFile} from "@/downloadUtils";

export default defineComponent({
  name: 'DownloadButtonSkaAlignment',
  components: {
    Button
  },
  setup() {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const {allResults_ska} = useState(["allResults_ska"]) as any;

    return {
      allResults_ska
    }
  },
  methods: {
    downloadALN(): void {
      downloadTextFile(this.allResults_ska.alignResults[0].alignment, "alignment.aln", "text/x-fasta");
    },
    downloadTREE(): void {
      downloadTextFile(this.allResults_ska.alignResults[0].newick, "alignmenttree.tree", "text/x-newick");
    },
    downloadCSV(): void {
      const csv = this.allResults_ska.alignResults[0]?.distances_csv;
      if (!csv) return;
      downloadTextFile(csv, "distances.csv", "text/csv");
    }
  }
});
</script>

<style scoped>
</style>
