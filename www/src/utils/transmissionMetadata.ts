import type { Dict, MetadataRow } from "@/types";

export function parseTransmissionMetadataCsv(text: string): MetadataRow[] {
    const lines = text.split(/\r?\n/).filter((line) => line.trim());
    if (lines.length < 2) {
        throw new Error("Metadata file is empty or has no data rows.");
    }

    const headers = lines[0].split(',').map((header) => header.trim().toLowerCase());
    const idIdx = headers.indexOf('id');
    const dateIdx = headers.indexOf('date');
    const locIdx = headers.indexOf('location');
    if (idIdx === -1 || dateIdx === -1 || locIdx === -1) {
        throw new Error("Metadata CSV must have columns: ID, date, location.");
    }

    const datePattern = /^\d{4}-\d{2}-\d{2}$/;
    const rows: MetadataRow[] = [];
    const errors: string[] = [];
    const ids = new Set<string>();

    for (let i = 1; i < lines.length; i++) {
        const cols = lines[i].split(',');
        const id = cols[idIdx]?.trim();
        const date = cols[dateIdx]?.trim() ?? '';
        const location = cols[locIdx]?.trim() ?? '';

        if (!id) {
            continue;
        }
        if (ids.has(id)) errors.push(`Row ${i + 1}: duplicate metadata ID "${id}"`);
        ids.add(id);
        if (date && !datePattern.test(date)) {
            errors.push(`Row ${i + 1}: invalid date "${date}" (expected yyyy-mm-dd)`);
        }
        rows.push({ id, date, location });
    }

    if (errors.length > 0) {
        throw new Error(errors.join('; '));
    }

    return rows;
}

// Strip only recognised terminal sequence extensions; internal dots are significant.
export function sequenceFileStem(name: string): string {
    return name.replace(/\.(?:fa|fna|fasta|fas|fq|fnq|fastq|aln)(?:\.gz)?$/i, "");
}

export interface SampleMetadataIndex {
    hasMetadata: boolean;
    rowsBySample: Map<string, MetadataRow>;
    ambiguousSamples: string[];
}

export function buildSampleMetadataIndex(sampleNames: string[], metadataRows: MetadataRow[]): SampleMetadataIndex {
    const exactRows = new Map(metadataRows.map(row => [row.id, row]));
    const samplesByStem = new Map<string, string[]>();
    const rowsByStem = new Map<string, MetadataRow[]>();
    for (const name of sampleNames) {
        const stem = sequenceFileStem(name);
        const names = samplesByStem.get(stem) ?? [];
        names.push(name);
        samplesByStem.set(stem, names);
    }
    for (const row of metadataRows) {
        const stem = sequenceFileStem(row.id);
        const rows = rowsByStem.get(stem) ?? [];
        rows.push(row);
        rowsByStem.set(stem, rows);
    }
    const rowsBySample = new Map<string, MetadataRow>();
    const ambiguousSamples: string[] = [];
    for (const name of sampleNames) {
        const exact = exactRows.get(name);
        if (exact) { rowsBySample.set(name, exact); continue; }
        const stem = sequenceFileStem(name);
        const matches = rowsByStem.get(stem) ?? [];
        if (!matches.length) continue;
        if (matches.length === 1 && samplesByStem.get(stem)?.length === 1) rowsBySample.set(name, matches[0]);
        else ambiguousSamples.push(name);
    }
    return { rowsBySample, ambiguousSamples, hasMetadata: metadataRows.length > 0 };
}

export function adjustClusterAssignments(
    clusterResults: Dict<number> | null,
    metadataIndex: SampleMetadataIndex,
    enableLocationMatching: boolean,
): Dict<number> | null {
    if (!clusterResults || !enableLocationMatching || !metadataIndex.hasMetadata) return clusterResults;
    const keyToId = new Map<string, number>();
    const result: Dict<number> = Object.create(null);
    let nextId = 1;
    for (const [name, cluster] of Object.entries(clusterResults)) {
        const location = metadataIndex.rowsBySample.get(name)?.location ?? "";
        const key = `${cluster}:${location}`;
        if (!keyToId.has(key)) keyToId.set(key, nextId++);
        result[name] = keyToId.get(key)!;
    }
    return result;
}
