#!/usr/bin/env node
// sparrowhawk-asm (Assembly tab): the worker's calls at the page's defaults, on the assembler's E. coli
// regression fixture.
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";
import { LocalFile, REPO, canonical, fastaRecords, golden, loadPackage, messages, options, report, reverseComplement,
    sha256 } from "./browser.mjs";

const values = options({ "native-single": { type: "string" }, "native-multik": { type: "string" } });

// The E. coli fixture, and the criteria of rust/sparrowhawk-asm/tests/assembly_regression.rs.
const fixture = join(REPO, "rust", "sparrowhawk-asm", "tests", "data");
const FLANK = 1_000, WINDOW = 20_000;
const reference = fastaRecords(readFileSync(join(fixture, "ecoli_window.fasta"), "utf8"))[0];
const gzipped = ["ecoli_window_1.fastq.gz", "ecoli_window_2.fastq.gz"].map((name) => readFileSync(join(fixture, name)));
const sortedCanonical = (contigs) => contigs.map(canonical).sort();

// The page's sequence of calls, at the CLI's defaults: Bloom counting, a fitted cutoff, automatic
// quality floors (no floor given), contigs of 500 nt or more.
function assemble(api, files, { ks = null, minQual = undefined } = {}) {
    messages.length = 0;
    const helper = api.AssemblyHelper.new(31, false, 5, minQual, true, true, 500);
    try {
        if (ks) helper.set_multik(new Uint32Array(ks));
        helper.preprocess(...files);
        const info = JSON.parse(helper.get_preprocessing_info());
        helper.assemble();
        const states = messages.filter((m) => m?.assemblyState).map((m) => m.assemblyState);
        return { info, result: JSON.parse(helper.get_assembly()), states };
    } finally {
        helper.free();
    }
}

/** What holds for any assembly the page shows: FASTA lines of at most 80 bases, as many records as reported. */
function checkFasta(label, result) {
    const lines = result.outfasta.split("\n").filter(Boolean);
    assert.ok(lines.every((line) => line.startsWith(">") || line.length <= 80), `${label}: a FASTA line over 80 bases`);
    const contigs = fastaRecords(result.outfasta);
    assert.equal(contigs.length, result.ncontigs, `${label}: ncontigs`);
    return contigs;
}

function check(label, { info, result, states }, multik) {
    assert.ok(info.nkmers > 0 && info.used_min_count >= 2 && info.histo.length <= 500, `${label}: preprocessing info`);
    assert.equal(info.multik, multik, `${label}: preprocessing info's multik`);
    assert.equal(info.k, multik ? 21 : 31, `${label}: the spectrum's k`);
    const contigs = checkFasta(label, result);
    const covered = new Uint8Array(reference.length);
    contigs.forEach((contig, i) => {
        let start = reference.indexOf(contig);
        if (start < 0) start = reference.indexOf(reverseComplement(contig));
        assert.ok(start >= 0, `${label}: contig ${i} (${contig.length} bp) is not an exact piece of the reference`);
        covered.fill(1, start, start + contig.length);
    });
    const windowCovered = covered.subarray(FLANK, FLANK + WINDOW).reduce((n, c) => n + c, 0);
    assert.ok(windowCovered >= 0.99 * WINDOW, `${label}: contigs cover ${windowCovered} of ${WINDOW} bp`);
    assert.ok(Math.max(...contigs.map((c) => c.length)) >= 19_000, `${label}: longest contig under 19,000 bp`);
    for (const state of ["preprocess:bulk:fitting", "preprocess:end", "assembly:correct_graph", "assembly:end"]) {
        assert.ok(states.includes(state), `${label}: no "${state}" progress`);
    }
    if (multik) {
        // A ladder stopped by a low peak writes no graph, so only its steps are checked.
        assert.ok(states.some((s) => s.startsWith("multik:ladder:21,31,")), `${label}: the ladder was not posted`);
        assert.ok(states.includes("multik:step:2:31"), `${label}: the ladder never reached its second k`);
    } else {
        assert.ok(result.outgfa.startsWith("H\tVN:Z:1.0") && result.outdot.length > 0, `${label}: graph text`);
    }
    return sortedCanonical(contigs);
}

// The fixture is under one batch of reads, so no progress with a percentage is posted for it. Ten
// copies of it, as concatenated gzip members, are: the percentage must rise and end near 100.
function checkProgress(label, api) {
    const copies = gzipped.map((bytes, i) =>
        new LocalFile(`copies_${i + 1}.fastq.gz`, Buffer.concat(Array(10).fill(bytes))));
    const percents = assemble(api, copies).states
        .map((s) => s.split(":"))
        .filter((parts) => parts.length === 5 && parts[2] === "loop")
        .map((parts) => Number(parts[4]));
    assert.ok(percents.length >= 3, `${label}: ${percents.length} progress states with a percentage`);
    assert.ok(percents.every((p, i) => i === 0 || p >= percents[i - 1]), `${label}: the percentage went down`);
    assert.ok(percents.at(-1) >= 80 && percents.at(-1) <= 100, `${label}: the last percentage is ${percents.at(-1)}`);
}

// Hash-map order differs between targets (32-bit usize on wasm32, another SIMD group width), so the
// graph's node order and its tie-breaks can differ: equal contigs are reported, not required.
function compare(label, ours, theirs) {
    if (!theirs) return;
    const same = ours.length === theirs.length && ours.every((s, i) => s === theirs[i]);
    report(`${label}: ${same ? "identical contigs" : `contigs differ (${ours.length} against ${theirs.length})`}`);
}
const summary = (contigs) => ({ contigs: contigs.length, sha256: sha256(contigs.join("\n")) });

const native = (path) => path && sortedCanonical(fastaRecords(readFileSync(path, "utf8")));
const nativeSingle = native(values["native-single"]);
const nativeMultik = native(values["native-multik"]);
let previous = null;
for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    const gzFiles = gzipped.map((bytes, i) => new LocalFile(`ecoli_window_${i + 1}.fastq.gz`, bytes));
    const plainFiles = gzipped.map((bytes, i) => new LocalFile(`ecoli_window_${i + 1}.fastq`, gunzipSync(bytes)));
    const single = check(`${target} k=31`, assemble(api, gzFiles), false);
    assert.deepEqual(check(`${target} k=31 uncompressed`, assemble(api, plainFiles), false), single,
        `${target}: gzipped and uncompressed input assemble differently`);
    const multik = check(`${target} multi-k`, assemble(api, gzFiles, { ks: [] }), true);
    // A given floor is the only one: it assembles, and is never recounted at a looser floor.
    const fixed = assemble(api, gzFiles, { minQual: 25 });
    assert.ok(fixed.result.ncontigs > 0 && !fixed.states.some((s) => s.startsWith("preprocess:recount:")),
        `${target}: a fixed base-quality floor of 25`);
    assert.throws(() => assemble(api, gzFiles, { ks: [20] }), /odd number/, `${target}: an even k was accepted`);
    checkProgress(target, api);
    golden(`assembler-${target}`, { single: summary(single), multik: summary(multik) }, values["update-golden"]);
    compare(`${target} k=31 against native`, single, nativeSingle);
    compare(`${target} multi-k against native`, multik, nativeMultik);
    if (previous) {
        compare(`${target} k=31 against ${previous.target}`, single, previous.single);
        compare(`${target} multi-k against ${previous.target}`, multik, previous.multik);
    }
    previous = { target, single, multik };
    report(`PASS assembler ${target} (${root})`);
}
