#!/usr/bin/env node
// orphos-bridge (Gene calling tab), with sparrowhawk-amr annotating the genes as Caller.ts does: the
// worker's calls at the page's defaults on S. pneumoniae R6 and on the tab's K. pneumoniae example.
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";
import { DATA, REPO, fastaRecords, fileFrom, golden, loadPackage, options, report, sha256 } from "./browser.mjs";

const values = options({ "amr-package": { type: "string" } });
assert.ok(values["amr-package"], "Pass --amr-package DIR: gene calling annotates with the AMR detector.");
const BGZF_EOF = Buffer.from("1f8b08040000000000ff0600424302001b0003000000000000000000", "hex");
const amr = await loadPackage(values["amr-package"]);
const detector = new amr.api.AmrDetector(gunzipSync(readFileSync(join(REPO, "www", "public",
    "amrfinderplus_2026-05-15.1_dna_k23_amr-stress-virulence.amridx.gz"))));

function callGenes(api, genome) {
    const orphos = api.OrphosData.new(false, "gff", false, false, false, 0);
    try {
        orphos.read_fasta(genome);
        orphos.index_fasta();
        orphos.call_genes();
        const metadata = JSON.parse(orphos.get_gene_metadata_json());
        const cds = orphos.get_cds_fasta();
        const amrJson = cds ? detector.detect_cds(genome.name, new TextEncoder().encode(cds), 0.3, 0.1) : null;
        const result = JSON.parse(amrJson ? orphos.get_annotated_results("gff", amrJson) : orphos.get_results("gff"));
        assert.ok(result.output_file.startsWith("##gff-version"), `${genome.name}: GFF header`);
        const cdsLines = result.output_file.split("\n").filter((line) => line.split("\t")[2] === "CDS").length;
        assert.equal(cdsLines, result.gene_count, `${genome.name}: CDS lines`);
        assert.equal(Object.keys(metadata).length, result.gene_count, `${genome.name}: gene metadata`);
        assert.equal(fastaRecords(cds).length, result.gene_count, `${genome.name}: CDS records`);
        for (const [label, bgz] of [["FASTA", orphos.take_fasta_bgz()], ["GFF", orphos.take_gff_bgz()]]) {
            assert.ok(Buffer.from(bgz).subarray(-BGZF_EOF.length).equals(BGZF_EOF), `${genome.name}: ${label} BGZF end-of-file block`);
        }
        for (const take of ["take_fasta_fai", "take_fasta_gzi", "take_gff_csi"]) {
            assert.ok(orphos[take]().length > 0, `${genome.name}: ${take}`);
        }
        const hits = amrJson ? JSON.parse(amrJson).hits.map((h) => h.gene_symbol).sort() : [];
        return { genome: genome.name, genes: result.gene_count, sequences: result.sequence_count,
            gff: sha256(result.output_file), amr: hits };
    } finally {
        orphos.free();
    }
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    golden(`genecalling-${target}`, [
        callGenes(api, fileFrom(join(REPO, "rust", "sketchlib.rust", "tests", "test_files_in", "R6.fa.gz"))),
        callGenes(api, fileFrom(join(DATA, "GCA_004138665.1.fasta.gz"))),
    ], values["update-golden"]);
    // An unknown translation table is a panic, which ends the instance: last.
    assert.throws(() => api.OrphosData.new(false, "gff", false, false, false, 7));
    report(`PASS genecalling ${target} (${root})`);
}
