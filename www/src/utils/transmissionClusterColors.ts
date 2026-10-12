import * as d3 from "d3";

export const TRANSMISSION_CLUSTER_PALETTE: readonly string[] = [
  ...d3.schemeTableau10,
  ...(d3.schemePastel1 as readonly string[]),
];

export function transmissionClusterColor(clusterId: number): string {
  const index = ((clusterId % TRANSMISSION_CLUSTER_PALETTE.length) + TRANSMISSION_CLUSTER_PALETTE.length)
    % TRANSMISSION_CLUSTER_PALETTE.length;
  return TRANSMISSION_CLUSTER_PALETTE[index];
}

type TaxoniumRgb = [number, number, number];

export function buildTaxoniumClusterColorMapping(clusterIds: number[]): Record<string, TaxoniumRgb> {
  const mapping: Record<string, TaxoniumRgb> = {};
  for (const clusterId of new Set(clusterIds)) {
    const hex = transmissionClusterColor(clusterId);
    mapping[`Cluster ${clusterId}`] = [
      Number.parseInt(hex.slice(1, 3), 16),
      Number.parseInt(hex.slice(3, 5), 16),
      Number.parseInt(hex.slice(5, 7), 16),
    ];
  }
  return mapping;
}
