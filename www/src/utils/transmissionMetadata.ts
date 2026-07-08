import type { MetadataRow } from "@/types";

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

    for (let i = 1; i < lines.length; i++) {
        const cols = lines[i].split(',');
        const id = cols[idIdx]?.trim();
        const date = cols[dateIdx]?.trim() ?? '';
        const location = cols[locIdx]?.trim() ?? '';

        if (!id) {
            continue;
        }
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
