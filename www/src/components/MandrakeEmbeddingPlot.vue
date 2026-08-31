<template>
  <div ref="plotElement" class="w-full" role="img" aria-label="Mandrake two-dimensional embedding" />
</template>

<script lang="ts">
import { defineComponent, onBeforeUnmount, onMounted, ref, watch } from "vue";
import { useResizeObserver } from "@vueuse/core";
// @ts-ignore - plotly.js-dist ships no types (as in the existing scatter plots).
import Plotly from "plotly.js-dist";

const SURFACE = "#fcfcfb";
const INK_PRIMARY = "#0b0b0b";
const INK_MUTED = "#898781";
const GRIDLINE = "#e1e0d9";
const AXIS = "#c3c2b7";
const PALETTE = [
  "#2563eb", "#db2777", "#059669", "#d97706", "#7c3aed",
  "#0891b2", "#dc2626", "#4f46e5", "#65a30d", "#c026d3",
];

export default defineComponent({
  name: "MandrakeEmbeddingPlot",
  props: {
    embedding: { type: Object as () => Float64Array, required: true },
    names: { type: Array as () => string[], required: true },
    labels: { type: Array as () => string[] | undefined, default: undefined },
    noiseLabel: { type: String, default: undefined },
    runKey: { type: Number, default: 0 },
  },
  setup(props) {
    const plotElement = ref<HTMLDivElement | null>(null);

    function render(): void {
      const element = plotElement.value;
      if (!element || props.embedding.length < 2) return;

      const pointCount = Math.floor(props.embedding.length / 2);
      const hasLabels = Boolean(props.labels && props.labels.length === pointCount);
      const groups = new Map<string, { x: number[]; y: number[]; names: string[] }>();
      const noise = { x: [] as number[], y: [] as number[], names: [] as string[] };

      for (let index = 0; index < pointCount; index += 1) {
        const x = props.embedding[index * 2];
        const y = props.embedding[index * 2 + 1];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        const label = hasLabels ? (props.labels?.[index] ?? "") : "Samples";
        const name = props.names[index] ?? `Sample ${index + 1}`;
        if (props.noiseLabel !== undefined && label === props.noiseLabel) {
          noise.x.push(x);
          noise.y.push(y);
          noise.names.push(name);
          continue;
        }
        const group = groups.get(label) ?? { x: [], y: [], names: [] };
        group.x.push(x);
        group.y.push(y);
        group.names.push(name);
        groups.set(label, group);
      }
      if (!groups.size && !noise.x.length) return;

      const labels = Array.from(groups.keys()).sort((left, right) => left.localeCompare(right));
      const traces: Array<Record<string, unknown>> = labels.map((label, index) => {
        const group = groups.get(label)!;
        return {
          type: "scattergl",
          mode: "markers",
          name: label || "(empty label)",
          x: [...group.x],
          y: [...group.y],
          text: [...group.names],
          hovertemplate: "<b>%{text}</b><br>%{fullData.name}<extra></extra>",
          marker: {
            size: 8,
            color: PALETTE[index % PALETTE.length],
            line: { color: SURFACE, width: 2 },
          },
        };
      });
      if (noise.x.length) {
        traces.push({
          type: "scattergl",
          mode: "markers",
          name: props.noiseLabel ?? "Noise",
          x: [...noise.x],
          y: [...noise.y],
          text: [...noise.names],
          hovertemplate: "<b>%{text}</b><br>%{fullData.name}<extra></extra>",
          marker: {
            size: 8,
            color: "#18212b",
            opacity: 0.45,
            line: { color: SURFACE, width: 1.5 },
          },
          showlegend: false,
        });
      }

      const axis = {
        zeroline: false,
        showgrid: true,
        gridcolor: GRIDLINE,
        gridwidth: 1,
        linecolor: AXIS,
        linewidth: 1,
        tickfont: { size: 10, color: INK_MUTED },
        title: { font: { size: 11, color: INK_MUTED } },
      };
      const layout = {
        height: 440,
        margin: { l: 52, r: hasLabels ? 18 : 16, t: 8, b: 46 },
        paper_bgcolor: SURFACE,
        plot_bgcolor: SURFACE,
        font: { family: "DM Sans", color: INK_PRIMARY },
        showlegend: hasLabels,
        hovermode: "closest",
        legend: { font: { size: 10 }, orientation: "h", y: -0.14 },
        uirevision: props.runKey,
        xaxis: { ...axis, title: { ...axis.title, text: "SCE dimension 1" } },
        yaxis: {
          ...axis,
          title: { ...axis.title, text: "SCE dimension 2" },
          scaleanchor: "x",
          scaleratio: 1,
        },
      };
      Plotly.react(element, traces, layout, {
        displaylogo: false,
        responsive: false,
        modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
      });
    }

    onMounted(render);
    watch(() => [props.embedding, props.names, props.labels, props.noiseLabel, props.runKey], render);
    useResizeObserver(plotElement, () => {
      if (plotElement.value) Plotly.Plots.resize(plotElement.value);
    });
    onBeforeUnmount(() => {
      if (plotElement.value) Plotly.purge(plotElement.value);
    });

    return { plotElement };
  },
});
</script>
