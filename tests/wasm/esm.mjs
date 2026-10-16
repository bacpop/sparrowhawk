#!/usr/bin/env node
// esm-bridge (Protein embeddings tab): the worker's calls on the CPU, as the page defaults to, with the
// site's weights and UMAP encoder, on esm-bridge's test proteins against their reference embeddings and
// on the first 50 proteins of the tab's smaller proteome (UP000000807, frozen in tests/wasm/data).
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";
import { DATA, REPO, golden, loadPackage, messages, options, report, sha256 } from "./browser.mjs";

const values = options();
const tests = join(REPO, "rust", "esm-bridge", "tests", "data");
const asset = (name) => gunzipSync(readFileSync(join(REPO, "www", "public", name)));

/** A float32 .npy, row-major: the header's length is the u16 at bytes 8-9. */
function readNpy(path) {
    const bytes = readFileSync(path);
    const body = bytes.subarray(10 + bytes.readUInt16LE(8));
    return new Float32Array(body.buffer.slice(body.byteOffset, body.byteOffset + body.byteLength));
}

/** Residues per record, without a trailing stop. */
const residues = (fasta) => fasta.toString().split(">").slice(1)
    .map((record) => record.split("\n").slice(1).join("").replace(/\*$/, "").length);

async function embed(embedder, fasta) {
    messages.length = 0;
    const batch = await embedder.embed_fasta(fasta);
    try {
        const meta = JSON.parse(batch.meta);
        const vectors = batch.vectors;
        const coords = batch.coords;
        const headers = fasta.toString().split("\n").filter((l) => l.startsWith(">")).map((l) => l.slice(1).split(/\s/)[0]);
        assert.equal(meta.n, headers.length, "one embedding per protein");
        assert.deepEqual(meta.ids, headers, "ids in FASTA order");
        assert.deepEqual(meta.truncated, residues(fasta).map((n) => n > 1022), "truncation flags");
        assert.equal(vectors.length, meta.n * 320);
        assert.equal(coords.length, meta.n * 2);
        assert.ok(vectors.every(Number.isFinite) && coords.every(Number.isFinite), "finite embeddings and coordinates");
        const progress = messages.filter((m) => m && "embedProgress" in m);
        assert.ok(progress.every((m, i) => i === 0 || m.done >= progress[i - 1].done), "progress only grows");
        assert.ok(progress.at(-1).done === meta.n && progress.at(-1).total === meta.n, "progress ends at n/n");
        return { meta, vectors, coords };
    } finally {
        batch.free();
    }
}

function fidelity(vectors, reference) {
    let worstCosine = 1;
    let maxAbs = 0;
    for (let row = 0; row < reference.length / 320; row++) {
        let dot = 0, a = 0, b = 0;
        for (let i = row * 320; i < (row + 1) * 320; i++) {
            dot += vectors[i] * reference[i];
            a += vectors[i] ** 2;
            b += reference[i] ** 2;
            maxAbs = Math.max(maxAbs, Math.abs(vectors[i] - reference[i]));
        }
        worstCosine = Math.min(worstCosine, dot / Math.sqrt(a * b));
    }
    return { worstCosine: Number(worstCosine.toFixed(6)), maxAbs: Number(maxAbs.toExponential(2)) };
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    // Embedder.ts on the CPU: no GPU, so its power preference and task limit are unused.
    const embedder = await api.EsmEmbedder.create(asset("esm2_t6_8M_UR50D.bpkz"), asset("esm2_umap_encoder.bpkz"), false, 0, 0);
    assert.equal(embedder.dim(), 320);
    assert.equal(embedder.backend(), "cpu");
    const observed = {};
    for (const [fasta, npy] of [["seqs.faa", "reference.npy"], ["proteins20.faa", "reference_official.npy"]]) {
        const bytes = readFileSync(join(tests, fasta));
        const first = await embed(embedder, bytes);
        // Determinism on the small file only: proteins20.faa's long proteins make a second pass slow.
        if (fasta === "seqs.faa") {
            const again = await embed(embedder, bytes);
            assert.deepEqual(again.vectors, first.vectors, `${fasta}: the same embeddings twice`);
        }
        observed[fasta] = fidelity(first.vectors, readNpy(join(tests, npy)));
        assert.ok(observed[fasta].worstCosine >= 0.999, `${fasta}: worst cosine ${observed[fasta].worstCosine} to ${npy}`);
    }
    const proteome = await embed(embedder, gunzipSync(readFileSync(join(DATA, "UP000000807_50.fasta.gz"))));
    observed.UP000000807 = {
        proteins: proteome.meta.n,
        // Rounded, so that only a real change of the model shows.
        vectors: sha256(Array.from(proteome.vectors, (v) => v.toFixed(4)).join(",")),
    };
    golden(`esm-${target}`, observed, values["update-golden"]);
    await assert.rejects(embedder.embed_fasta(new TextEncoder().encode("MKV\n")));
    embedder.free();
    report(`PASS esm ${target} (${root})`);
}
