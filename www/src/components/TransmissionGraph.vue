<template>
  <div class="relative w-full h-full border border-gray-200 rounded-md bg-white overflow-hidden">
    <svg ref="svgRef" class="w-full h-full" />
    <!-- Legend -->
    <div v-if="clusterIds.length > 0"
         class="absolute top-3 right-3 bg-white bg-opacity-90 border border-gray-200 rounded-md p-2 text-xs space-y-1 max-h-48 overflow-y-auto">
      <div v-for="clusterId in clusterIds" :key="clusterId" class="flex items-center gap-1.5">
        <span class="inline-block w-3 h-3 rounded-full flex-shrink-0" :style="{ backgroundColor: clusterColor(clusterId) }" />
        <span>Cluster {{ clusterId }}</span>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, watch, onMounted, onUnmounted, computed, PropType } from "vue";
import * as d3 from "d3";
import type { TransmissionGraphNode, TransmissionGraphLink } from "@/types";

export default defineComponent({
  name: "TransmissionGraph",
  props: {
    nodes: {
      type: Array as PropType<TransmissionGraphNode[]>,
      required: true,
    },
    links: {
      type: Array as PropType<TransmissionGraphLink[]>,
      required: true,
    },
  },
  setup(props) {
    const svgRef = ref<SVGSVGElement | null>(null);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    let simulation: any = null;

    const palette = d3.schemeTableau10.concat(d3.schemePastel1 as readonly string[]);

    const clusterIds = computed(() => {
      const ids = new Set(props.nodes.map(n => n.cluster));
      return Array.from(ids).sort((a, b) => a - b);
    });

    function clusterColor(clusterId: number): string {
      return palette[clusterId % palette.length];
    }

    function buildGraph() {
      if (!svgRef.value) return;

      const svg = d3.select(svgRef.value);
      svg.selectAll("*").remove();

      if (props.nodes.length === 0) return;

      const width = svgRef.value.clientWidth || 600;
      const height = svgRef.value.clientHeight || 500;

      // Deep-copy so D3 can mutate positions without touching Vue props
      // fx/fy are D3 simulation fixed-position fields
      const nodesCopy = props.nodes.map(n => ({ ...n, x: width / 2, y: height / 2, fx: null as number | null, fy: null as number | null }));
      const nodeById = new Map(nodesCopy.map(n => [n.id, n]));
      const linksCopy = props.links
        .filter(l => nodeById.has(l.source) && nodeById.has(l.target))
        .map(l => ({ ...l, source: nodeById.get(l.source)!, target: nodeById.get(l.target)! }));

      const maxDist = Math.max(...props.links.map(l => l.snp_distance), 1);

      const g = svg.append("g");

      // Zoom/pan
      const zoomBehavior = d3.zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.2, 5])
        .on("zoom", (event) => g.attr("transform", event.transform));
      svg.call(zoomBehavior);

      // Links
      const link = g.append("g")
        .selectAll("line")
        .data(linksCopy)
        .join("line")
        .attr("stroke", "#94a3b8")
        .attr("stroke-width", 1.5)
        .attr("stroke-opacity", d => 0.3 + 0.7 * (1 - d.snp_distance / maxDist));

      // Link labels
      const linkLabel = g.append("g")
        .selectAll("text")
        .data(linksCopy)
        .join("text")
        .attr("font-size", 9)
        .attr("fill", "#64748b")
        .attr("text-anchor", "middle")
        .text(d => d.snp_distance);

      // Node circles
      const node = g.append("g")
        .selectAll("circle")
        .data(nodesCopy)
        .join("circle")
        .attr("r", 12)
        .attr("fill", d => clusterColor(d.cluster))
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5)
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        .call(d3.drag<SVGCircleElement, typeof nodesCopy[0]>()
            .on("start", (event, d) => {
              if (!event.active) simulation.alphaTarget(0.3).restart();
              d.fx = d.x; d.fy = d.y;
            })
            .on("drag", (event, d) => { d.fx = event.x; d.fy = event.y; })
            .on("end", (event, d) => {
              if (!event.active) simulation.alphaTarget(0);
              d.fx = null; d.fy = null;
            }) as any
        );

      // Node labels
      const label = g.append("g")
        .selectAll("text")
        .data(nodesCopy)
        .join("text")
        .attr("font-size", 10)
        .attr("fill", "#1e293b")
        .attr("dy", "0.35em")
        .attr("x", 11)
        .text(d => d.id.length > 20 ? d.id.slice(0, 18) + "…" : d.id);

      // Tooltip
      node.append("title").text(d => `${d.id} (Cluster ${d.cluster})`);

      simulation = d3.forceSimulation(nodesCopy)
        .force("link", d3.forceLink(linksCopy).id((d: any) => d.id).distance(40))
        .force("charge", d3.forceManyBody().strength(-50).distanceMax(100))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("collision", d3.forceCollide(24))
        .on("tick", () => {
          link
            .attr("x1", d => (d.source as any).x)
            .attr("y1", d => (d.source as any).y)
            .attr("x2", d => (d.target as any).x)
            .attr("y2", d => (d.target as any).y);

          linkLabel
            .attr("x", d => ((d.source as any).x + (d.target as any).x) / 2)
            .attr("y", d => ((d.source as any).y + (d.target as any).y) / 2);

          node
            .attr("cx", d => d.x!)
            .attr("cy", d => d.y!);

          label
            .attr("x", d => d.x! + 11)
            .attr("y", d => d.y!);
        })
        .on("end", () => {
          if (nodesCopy.length === 0) return;
          const xs = nodesCopy.map(d => d.x ?? 0);
          const ys = nodesCopy.map(d => d.y ?? 0);
          const minX = Math.min(...xs), maxX = Math.max(...xs);
          const minY = Math.min(...ys), maxY = Math.max(...ys);
          const padding = 50;
          const scale = Math.min(
            (width - 2 * padding) / (maxX - minX || 1),
            (height - 2 * padding) / (maxY - minY || 1),
            3
          );
          const tx = width / 2 - scale * (minX + maxX) / 2;
          const ty = height / 2 - scale * (minY + maxY) / 2;
          svg.call(zoomBehavior.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
        });
    }

    onMounted(() => { buildGraph(); });
    onUnmounted(() => { if (simulation) simulation.stop(); });
    watch([() => props.nodes, () => props.links], () => {
      if (simulation) simulation.stop();
      buildGraph();
    });

    return { svgRef, clusterIds, clusterColor };
  },
});
</script>
