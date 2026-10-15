<template>
  <div :class="showPlot ? 'my-6 border-t border-t-gray-200 pt-6' : 'hidden'">
    <div ref="plotElement"
         id="plotElement"
         class="w-full">
    </div>
    <p v-if="readsPreprocessing.multik" class="mx-6 mt-2 text-sm text-gray-500">
      Multi-k: this is the spectrum of the ladder's first k ({{ readsPreprocessing.k }}). The later k values
      are counted during the assembly, and are not plotted.
    </p>
  </div>
</template>

<script lang="ts">
import {defineComponent, ref, Ref, useTemplateRef} from "vue";
import {useStore} from "vuex";
import {useState} from "vuex-composition-helpers";
// @ts-ignore - plotly.js-dist doesn't have types
import Plotly from 'plotly.js-dist';

export default defineComponent({
  name: 'KmerHistogram',

  setup() {
    const store = useStore();
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const {readsPreprocessing} = useState(["readsPreprocessing"]) as any;
    const skip: Ref<boolean> = ref(true);
    const plotElement = useTemplateRef<HTMLDivElement>('plotElement')
    return {
      readsPreprocessing,
      skip,
      store,
      plotElement
    }
  },

  data() {
    return {};
  },

  methods: {
    renderChart(): void {
      console.log("Drawing!");

      const elWidth = Number(getComputedStyle(this.plotElement as HTMLDivElement).width.replace('px', ''));

      const inputdata = [{
        x: Array.from({length: 500}, (_, i) => i + 1),
        y: this.readsPreprocessing.histo,
        type: "bar",
        marker: {
          color: "#18974C",
        },
      }];
      const chartLayout = {
        title: {text: this.readsPreprocessing.k ? `k-mer spectrum (k = ${this.readsPreprocessing.k})` : "k-mer spectrum"},
        font: {
          family: "IBM Plex Sans",
        },
        autosize: true,
        xaxis: {
          title: {text: "k-mer frequency"},
          rangemode: 'nonnegative',
        },
        yaxis: {
          title: {text: "Counts"},
          rangemode: 'nonnegative',
        },
        width: elWidth,
        shapes: [{
          type: 'rect',
          x0: 0,
          y0: 0,
          x1: (this.readsPreprocessing.used_min_count || 0) + 0.5,
          y1: 1,
          xref: "x",
          yref: "paper",
          line: {
            color: '#373A36',
            width: 1.5,
            dash: 'dot'
          },
          fillcolor: "rgba(55,58,54,0.5)",
        }],
      };

      Plotly.newPlot("plotElement", inputdata, chartLayout as Partial<Plotly.Layout>);
    },
  },

  computed: {
    queryPreprocessed(): boolean {
      return this.store.getters.readsPreprocessed;
    },
    isAssembling(): boolean {
      return this.store.getters.isAssembling;
    },
    showPlot(): boolean {
      return this.queryPreprocessed && !this.isAssembling;
    },
  },

  mounted(): void {
  },

  beforeUnmount(): void {
  },

  watch: {
    // After the DOM update, so the plot is drawn into a visible element of known width. Clearing the
    // results removes the plot rather than drawing an empty one.
    showPlot: {
      handler(show: boolean): void {
        if (show) {
          this.renderChart();
        } else if (!this.queryPreprocessed) {
          Plotly.purge(this.plotElement as HTMLDivElement);
        }
      },
      flush: 'post',
    },
  }
});
</script>

<style scoped>
</style>
