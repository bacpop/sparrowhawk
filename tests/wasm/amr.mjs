#!/usr/bin/env node
// sparrowhawk-amr (AMR detection tab): the worker's calls at the page's defaults, with the site's index,
// on S. pneumoniae R6 and on the K. pneumoniae genome the other tabs offer.
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";
import { DATA, REPO, golden, loadPackage, options, report } from "./browser.mjs";

const values = options();

function detect(detector, path, sampleName) {
    // AmrDetector.ts gunzips .gz files with fflate before detect_direct.
    const result = JSON.parse(detector.detect_direct(sampleName, gunzipSync(readFileSync(path)), 0.3, 0.1));
    assert.equal(result.query_kind, "Direct");
    assert.equal(result.gene_count, result.hits.filter((h) => h.call_type === "gene").length, `${sampleName}: gene count`);
    assert.equal(result.gene_group_count, result.hits.filter((h) => h.call_type === "gene_group").length,
        `${sampleName}: gene group count`);
    return { sample: sampleName, hits: result.hits.map((h) => [h.gene_symbol, h.call_type, h.class_name]) };
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    const detector = new api.AmrDetector(gunzipSync(readFileSync(join(REPO, "www", "public",
        "amrfinderplus_2026-05-15.1_dna_k23_amr-stress-virulence.amridx.gz"))));
    const info = detector.info();
    assert.match(info, /alphabet=dna/);
    assert.match(info, /k=23/);
    golden(`amr-${target}`, [
        detect(detector, join(REPO, "rust", "sketchlib.rust", "tests", "test_files_in", "R6.fa.gz"), "R6"),
        detect(detector, join(DATA, "GCA_004138665.1.fasta.gz"), "GCA_004138665.1"),
    ], values["update-golden"]);
    assert.throws(() => detector.detect_direct("x", new TextEncoder().encode(">x\nACGT\n"), 1.5, 0.1), /between 0 and 1/);
    report(`PASS amr ${target} (${root})`);
}
