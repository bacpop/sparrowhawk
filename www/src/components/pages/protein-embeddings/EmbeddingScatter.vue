<template>
  <div ref="plotElement" class="w-full" />
</template>

<script lang="ts">
import { defineComponent, onBeforeUnmount, onMounted, ref, watch } from "vue";
import { useResizeObserver } from "@vueuse/core";
// @ts-ignore - plotly.js-dist ships no types (same as KmerHistogram.vue)
import Plotly from "plotly.js-dist";

// Light-only on purpose: the app has no dark mode anywhere.
const SURFACE = "#fcfcfb";
const INK_PRIMARY = "#0b0b0b";
const INK_MUTED = "#898781";
const GRIDLINE = "#e1e0d9";
const AXIS = "#c3c2b7";

/**
 * Length is a magnitude, so the ramp is sequential: one hue, light to dark, never a rainbow.
 * Starts at step 250 rather than 100 because in a scatter the dot *is* the datum — a paler
 * step would be invisible against the surface.
 */
const LENGTH_RAMP: [number, string][] = [
  [0.0, "#86b6ef"],
  [0.333, "#3987e5"],
  [0.667, "#1c5cab"],
  [1.0, "#0d366b"],
];

export default defineComponent({
  name: "EmbeddingScatter",

  props: {
    /** Row-major [n, 2]. */
    coords: {
      type: Object as () => Float32Array,
      required: true,
    },
    ids: {
      type: Array as () => string[],
      required: true,
    },
    /** Residues per protein — the colour channel. */
    lengths: {
      type: Array as () => number[],
      required: true,
    },
  },

  setup(props) {
    const plotElement = ref<HTMLDivElement | null>(null);

    function render(): void {
      const el = plotElement.value;
      if (!el || props.coords.length === 0) return;

      const n = Math.floor(props.coords.length / 2);
      const x = new Array<number>(n);
      const y = new Array<number>(n);
      for (let i = 0; i < n; i++) {
        x[i] = props.coords[i * 2];
        y[i] = props.coords[i * 2 + 1];
      }

      // Copies, not the props: Plotly writes into the arrays it is given, and a readonly
      // reactive proxy refuses every write with one console warning per index.
      const ids = Array.from(props.ids);
      const lengths = Array.from(props.lengths);

      const trace = {
        // WebGL, not SVG: a proteome is several thousand points.
        type: "scattergl",
        mode: "markers",
        x,
        y,
        text: ids,
        customdata: lengths,
        hovertemplate: "<b>%{text}</b><br>%{customdata} aa<extra></extra>",
        marker: {
          size: 8,
          color: lengths,
          // Copied too: LENGTH_RAMP is module-level, so a mutation would outlive the plot.
          colorscale: LENGTH_RAMP.map((stop) => [...stop]),
          // The scale legend, and the only legend needed: there is one series.
          colorbar: {
            title: { text: "Length (aa)", font: { size: 11, color: INK_MUTED } },
            thickness: 10,
            outlinewidth: 0,
            tickfont: { size: 10, color: INK_MUTED },
          },
          // A surface-coloured ring, not a border: keeps overlapping dots legible.
          line: { color: SURFACE, width: 2 },
        },
      };

      // Tick values hidden on purpose: UMAP axes have no units, and showing them invites
      // reading distances off the plot.
      const axis = {
        zeroline: false,
        showticklabels: false,
        ticks: "" as const,
        showgrid: true,
        gridcolor: GRIDLINE,
        gridwidth: 1,
        linecolor: AXIS,
        linewidth: 1,
        title: { font: { size: 11, color: INK_MUTED } },
      };

      const layout = {
        // Includes the axis band, so the card never grows a nested scrollbar.
        height: 440,
        margin: { l: 48, r: 16, t: 8, b: 44 },
        paper_bgcolor: SURFACE,
        plot_bgcolor: SURFACE,
        font: { family: "IBM Plex Sans", color: INK_PRIMARY },
        showlegend: false,
        // Nearest-point, so an 8px dot need not be landed on dead centre.
        hovermode: "closest",
        xaxis: { ...axis, title: { ...axis.title, text: "UMAP 1" } },
        yaxis: { ...axis, title: { ...axis.title, text: "UMAP 2" }, scaleanchor: "x", scaleratio: 1 },
      };

      Plotly.react(el, [trace], layout, {
        displaylogo: false,
        responsive: false,
        modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
      });
    }

    onMounted(render);
    // `Plotly.react` diffs rather than rebuilding the WebGL context.
    watch(() => [props.coords, props.ids, props.lengths], render, { deep: false });

    useResizeObserver(plotElement, () => {
      if (plotElement.value) Plotly.Plots.resize(plotElement.value);
    });

    onBeforeUnmount(() => {
      // Plotly holds a WebGL context per plot; without this it outlives the component.
      if (plotElement.value) Plotly.purge(plotElement.value);
    });

    return { plotElement };
  },
});
</script>
