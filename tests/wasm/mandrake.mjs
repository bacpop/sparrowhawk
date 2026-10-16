#!/usr/bin/env node
// mandrake (Clustering tab): the worker's calls at the page's defaults, except far fewer updates, on the
// examples its help panel points to (mandrake's fixtures; 300 of its 5000 HIV sequences). The embedding is
// seeded from the clock, so only counts and shapes are checked.
import assert from "node:assert/strict";
import { join } from "node:path";
import { DATA, REPO, fileFrom, golden, loadPackage, options, report } from "./browser.mjs";

const values = options();
const fixtures = join(REPO, "rust", "mandrake", "tests", "fixtures");
// mode, value, perplexity, max updates (the page uses 1,000,000), repulsion samples, learning rate,
// initial exaggeration.
const SETTINGS = ["knn", 15, 30, 20_000, 5, 1, false];

// Mandrake.ts: distances in rounds of a hundredth, then the embedding in rounds of at least 64 updates.
function run(api, operation) {
    try {
        let distances = operation.advanceDistances(0);
        while (!distances.complete) distances = operation.advanceDistances(Math.max(1, Math.ceil(distances.maximum / 100)));
        operation.beginEmbedding();
        let progress = operation.advance(0);
        while (!progress.complete) progress = operation.advance(Math.max(64, Math.ceil(progress.maximum / 100)));
        assert.equal(progress.completed, progress.maximum, "the embedding ran every update");
        const names = operation.names().split("\n");
        const embedding = operation.embedding();
        assert.equal(embedding.length, 2 * names.length, "two coordinates per sample");
        assert.ok(embedding.every(Number.isFinite), "finite coordinates");
        const clusters = api.clusterEmbedding(embedding);
        assert.equal(clusters.length, names.length, "one cluster label per sample");
        assert.ok(clusters.every((label) => Number.isInteger(label) && label >= -1), "cluster labels");
        return names.length;
    } finally {
        operation.free();
    }
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    const skm = fileFrom(join(fixtures, "sketches.skm"));
    const skd = fileFrom(join(fixtures, "sketches.skd"));
    assert.deepEqual(Array.from(api.sketchKmerLengths(skm)), [21]);
    // One k-mer length allows only Jaccard distances, as the page chooses.
    assert.throws(() => api.MandrakeOperation.fromSketchFiles(skm, skd, ...SETTINGS, "core", 21), /at least two k-mer lengths/);
    const observed = {
        sketches: run(api, api.MandrakeOperation.fromSketchFiles(skm, skd, ...SETTINGS, "jaccard", 21)),
        accessory: run(api, api.MandrakeOperation.fromAccessoryFile(fileFrom(join(fixtures, "gene_presence_absence.Rtab.gz")), ...SETTINGS)),
        alignment: run(api, api.MandrakeOperation.fromAlignmentFile(fileFrom(join(DATA, "hiv_300.fas.gz")), ...SETTINGS)),
    };
    assert.equal(observed.alignment, 300, "HIV sequences");
    golden(`mandrake-${target}`, observed, values["update-golden"]);
    report(`PASS mandrake ${target} (${root})`);
}
