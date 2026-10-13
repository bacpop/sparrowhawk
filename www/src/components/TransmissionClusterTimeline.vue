<template>
  <div ref="plotEl" class="w-full" style="min-height: 300px;"></div>
</template>

<script lang="ts">
import { defineComponent, PropType } from "vue";
import type { MetadataRow } from "@/types";
import { buildSampleMetadataIndex } from "@/utils/transmissionMetadata";
import { transmissionClusterColor } from "@/utils/transmissionClusterColors";
// @ts-ignore - plotly.js-dist doesn't have types
import Plotly from 'plotly.js-dist';

const BASE_HEIGHT = 300;
const CROWDED_CLUSTERS = 32;

// Linear growth from 1x at 32 clusters up to 2x at 64 or more.
function timelineHeight(clusterCount: number): number {
  return Math.round(BASE_HEIGHT * Math.min(2, Math.max(1, clusterCount / CROWDED_CLUSTERS)));
}

const MIN_MARKER_DIAMETER = 8;
const MAX_MARKER_DIAMETER = 36;

function markerDiameter(cases: number, maximumCases: number): number {
  if (maximumCases <= 1) return MIN_MARKER_DIAMETER;
  const relativeCases = (cases - 1) / (maximumCases - 1);
  return MIN_MARKER_DIAMETER + (MAX_MARKER_DIAMETER - MIN_MARKER_DIAMETER)
    * Math.sqrt(Math.max(0, Math.min(1, relativeCases)));
}

interface ClusterRow {
  name: string;
  cluster: number;
  singleton?: boolean;
}

interface PreparedPoint {
  cluster: number;
  date: string;
  cluster_size: number;
  cases: number;
  text: string;
}

export default defineComponent({
  name: 'TransmissionClusterTimeline',

  props: {
    clusterRows: {
      type: Array as PropType<ClusterRow[]>,
      required: true,
    },
    metadataRows: {
      type: Array as PropType<MetadataRow[]>,
      required: true,
    },
  },

  computed: {
    preparedData(): PreparedPoint[] {
      const metadataIndex = buildSampleMetadataIndex(this.clusterRows.map(row => row.name), this.metadataRows);
      const joined: { cluster: number; date: string }[] = [];
      // The index covers every sample so ID ambiguity is judged as in the table; singletons are only skipped here.
      for (const row of this.clusterRows) {
        if (row.singleton) continue;
        const meta = metadataIndex.rowsBySample.get(row.name);
        if (meta?.date) {
          joined.push({ cluster: row.cluster, date: meta.date });
        }
      }
      if (joined.length === 0) return [];

      // cluster_size: count per cluster
      const clusterSize = new Map<number, number>();
      for (const { cluster } of joined) {
        clusterSize.set(cluster, (clusterSize.get(cluster) ?? 0) + 1);
      }

      // cases: count per (cluster, date)
      const casesMap = new Map<string, number>();
      for (const { cluster, date } of joined) {
        const key = `${cluster}::${date}`;
        casesMap.set(key, (casesMap.get(key) ?? 0) + 1);
      }

      // Deduplicate to one row per (cluster, date)
      const seen = new Set<string>();
      const result: PreparedPoint[] = [];
      for (const { cluster, date } of joined) {
        const key = `${cluster}::${date}`;
        if (seen.has(key)) continue;
        seen.add(key);
        const cluster_size = clusterSize.get(cluster) ?? 0;
        const cases = casesMap.get(key) ?? 0;
        const text = `Cluster ${cluster}\nCluster size: ${cluster_size}\nDate: ${date}\nCases: ${cases}`;
        result.push({ cluster, date, cluster_size, cases, text });
      }

      return result.sort((a, b) => a.cluster - b.cluster || a.date.localeCompare(b.date));
    },
  },

  methods: {
    renderChart(): void {
      const data = this.preparedData;
      if (!this.$refs.plotEl) return;

      if (!data.length) {
        Plotly.purge(this.$refs.plotEl as HTMLElement);
        return;
      }

      const clusters = [...new Set(data.map(d => d.cluster))].sort((a, b) => a - b);
      const maximumCases = data.reduce((maximum, point) => Math.max(maximum, point.cases), 1);
      const traces: object[] = [];

      for (const cl of clusters) {
        const points = data.filter(d => d.cluster === cl).sort((a, b) => a.date.localeCompare(b.date));
        const xs = points.map(p => p.date);
        const ys = points.map(() => String(cl));

        // Line trace connecting points within the cluster
        traces.push({
          x: xs,
          y: ys,
          type: 'scatter',
          mode: 'lines',
          line: { color: transmissionClusterColor(cl), width: 1.8 },
          showlegend: false,
          hoverinfo: 'skip',
        });

        // Marker trace sized by cases
        traces.push({
          x: xs,
          y: ys,
          type: 'scatter',
          mode: 'markers',
          marker: {
            size: points.map(p => markerDiameter(p.cases, maximumCases)),
            sizemode: 'diameter',
            sizeref: 1,
            color: transmissionClusterColor(cl),
            showscale: false,
          },
          text: points.map(p => p.text),
          hoverinfo: 'text',
          name: `Cluster ${cl}`,
          showlegend: false,
        });
      }

      const layout = {
        xaxis: { title: 'Isolation date', type: 'date', tickformat: '%Y-%m', tickangle: 45 },
        yaxis: { title: 'Cluster', type: 'category' },
        showlegend: false,
        font: { family: 'IBM Plex Sans' },
        autosize: true,
        height: timelineHeight(clusters.length),
        margin: { t: 30, r: 20 },
      };

      Plotly.react(this.$refs.plotEl as HTMLElement, traces, layout as Partial<Plotly.Layout>);
    },
  },

  mounted(): void {
    this.$nextTick(() => this.renderChart());
  },

  watch: {
    preparedData: {
      handler(): void {
        this.$nextTick(() => this.renderChart());
      },
      deep: true,
    },
  },
});
</script>
