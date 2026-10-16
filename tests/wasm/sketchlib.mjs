#!/usr/bin/env node
// sketchlib.rust (Taxonomic ID tab): the worker's calls at the page's defaults, with the site's index, on two
// S. pneumoniae genomes from sketchlib's tests and on the tab's K. pneumoniae example.
import assert from "node:assert/strict";
import { join } from "node:path";
import { DATA, REPO, fileFrom, golden, loadPackage, options, report } from "./browser.mjs";

const values = options();
const tests = join(REPO, "rust", "sketchlib.rust", "tests", "test_files_in");

function identify(data, sample) {
    data.query(sample, null, 1, 5, 20);
    const { ani, ranks, names, metadata } = JSON.parse(data.get_ani(3));
    assert.ok(ani.length >= 1 && [ranks, names, metadata].every((a) => a.length === ani.length), `${sample.name}: hits`);
    assert.equal(ranks[0], 1, `${sample.name}: the first hit's rank`);
    assert.ok(ani.every((a, i) => a > 0 && a <= 1 && (i === 0 || a <= ani[i - 1])), `${sample.name}: ANI order`);
    return { sample: sample.name, top: metadata[0], ani: ani.map((a) => Number(a.toFixed(6))), ranks };
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    // Sketcher.ts builds one SketchlibData per worker and queries it for every sample.
    const data = api.SketchlibData.new(fileFrom(join(REPO, "www", "public", "inverted_k_17_ss_50.ski")));
    const observed = ["R6.fa.gz", "TIGR4.fa.gz"].map((name) => identify(data, fileFrom(join(tests, name))));
    const klebsiella = identify(data, fileFrom(join(DATA, "GCA_004138665.1.fasta.gz")));
    assert.match(klebsiella.top, /Klebsiella/, "the K. pneumoniae example's top hit");
    golden(`sketchlib-${target}`, [...observed, klebsiella], values["update-golden"]);
    // The page also accepts .fna, which the crate rejects with a panic: this ends the instance, so it runs last.
    assert.throws(() => data.query(fileFrom(join(tests, "R6.fa.gz"), "R6.fna"), null, 1, 5, 20));
    report(`PASS sketchlib ${target} (${root})`);
}
