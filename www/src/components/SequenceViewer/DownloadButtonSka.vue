<template>
  <div class="flex flex-col gap-2">
    Download alignment
    <Button class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadALN">
      .aln
    </Button>
  </div>
</template>

<script lang="ts">
import {defineComponent} from "vue";
import {useState} from "vuex-composition-helpers";
import {Button} from "@/components/ui/button";
import {downloadTextFile} from "@/downloadUtils";

export default defineComponent({
  name: 'DownloadButtonSka',
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
      let text = "";
      let mapping_alignment = "";

      for (const file of Object.keys(this.allResults_ska.mapResults)) {
        text += `>${file}\n`;
        mapping_alignment = "";
        const sequences = this.allResults_ska.mapResults[file].mapped_sequences;
        if (sequences) {
          for (const sequence of sequences) {
            mapping_alignment += sequence;
          }
        }
        text += mapping_alignment + "\n";
      }

      downloadTextFile(text, "alignment.aln", "text/x-fasta");
    }
  }
});
</script>

<style scoped>
</style>
