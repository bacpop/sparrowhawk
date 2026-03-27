
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
      const filename = 'alignment.aln';
      const element = document.createElement('a');
      element.setAttribute('href', 'data:application/json;charset=utf-8,' + encodeURIComponent(this.allResults_ska.alignResults[0].alignment));
      element.setAttribute('download', filename);

      element.style.display = 'none';
      document.body.appendChild(element);

      element.click();
      document.body.removeChild(element);
    },
    downloadTREE(): void {
      const filename = 'alignmenttree.tree';
      const element = document.createElement('a');
      element.setAttribute('href', 'data:application/json;charset=utf-8,' + encodeURIComponent(this.allResults_ska.alignResults[0].newick));
      element.setAttribute('download', filename);

      element.style.display = 'none';
      document.body.appendChild(element);

      element.click();
      document.body.removeChild(element);
    },
    downloadCSV(): void {
      const csv = this.allResults_ska.alignResults[0]?.distances_csv;
      if (!csv) return;
      const element = document.createElement('a');
      element.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv));
      element.setAttribute('download', 'distances.csv');
      element.style.display = 'none';
      document.body.appendChild(element);
      element.click();
      document.body.removeChild(element);
    }
  }
});
</script>

<style scoped>
</style>
