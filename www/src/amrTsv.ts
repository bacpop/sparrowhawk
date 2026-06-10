import { AmrDetectionResult, GeneMetadataMap } from "@/types";

const BASE_HEADERS = [
    "Sample", "Database version", "Index alphabet", "Index k",
    "Query", "Unit ID", "Unit label", "Call",
    "Category", "Subtype",
    "Element symbol", "Gene symbol", "Allele symbol", "Gene group",
    "Hierarchy node", "Class", "Subclass", "Call fraction",
    "Diagnostic k-mers matched", "Diagnostic k-mers total",
];

function baseRow(result: AmrDetectionResult, hit: AmrDetectionResult["hits"][number]): string[] {
    return [
        result.sample_name,
        result.database_version,
        result.index_alphabet,
        String(result.index_k),
        hit.query_id,
        hit.unit_id,
        hit.unit_label,
        hit.call_type,
        hit.type_name ?? '',
        hit.subtype ?? '',
        hit.element_symbol ?? '',
        hit.gene_symbol ?? '',
        hit.allele_symbol ?? '',
        hit.gene_group,
        hit.hierarchy_node ?? '',
        hit.class_name ?? '',
        hit.subclass ?? '',
        hit.call_fraction.toFixed(4),
        String(hit.first_pass_distinct),
        String(hit.first_pass_diagnostic_total),
    ];
}

export function buildAmrTsv(result: AmrDetectionResult): string {
    const headers = [
        ...BASE_HEADERS,
        "Start", "End", "Member count",
    ];
    const rows = result.hits.map((hit) => [
        ...baseRow(result, hit),
        String(hit.start),
        String(hit.end),
        String(hit.member_count),
    ]);
    return [headers, ...rows].map(r => r.join("\t")).join("\n");
}

export function buildGeneCallingAmrTsv(result: AmrDetectionResult, geneMetadata: GeneMetadataMap): string {
    const headers = [
        ...BASE_HEADERS,
        "Contig", "Gene start", "Gene end", "Strand",
        "Match start in query", "Match end in query", "Member count",
    ];
    const rows = result.hits.map((hit) => {
        const metadata = geneMetadata[hit.query_id];
        return [
            ...baseRow(result, hit),
            metadata?.contig ?? '',
            metadata ? String(metadata.start) : '',
            metadata ? String(metadata.end) : '',
            metadata?.strand ?? '',
            String(hit.start),
            String(hit.end),
            String(hit.member_count),
        ];
    });
    return [headers, ...rows].map(r => r.join("\t")).join("\n");
}
