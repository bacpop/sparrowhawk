#!/usr/bin/env node
// ska.rust (Mapping and Alignment tabs): the workers' calls at the page's defaults, on ska's own test
// genomes and on the assembler's E. coli fixture reads.
import assert from "node:assert/strict";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";
import { REPO, fastaRecords, fileFrom, golden, loadPackage, messages, options, report, sha256 } from "./browser.mjs";

const values = options();
const files = join(REPO, "rust", "ska.rust", "tests", "test_files_in");
const ecoli = join(REPO, "rust", "sparrowhawk-asm", "tests", "data");

// Mapper.ts: one SkaData per reference, then map() per sample. min_count is 0 when every file is FASTA,
// 5 when any is FASTQ.
function mapping(api, reference, samples, minCount) {
    const ska = api.SkaData.new(reference, 31, true, false, false);
    try {
        const chromosomes = ska.get_reference().split("\n");
        const mapped = samples.map(([name, file1, file2]) => {
            const out = JSON.parse(ska.map(file1, file2, 1, minCount, 20, 2, name));
            out["Mapped sequences"].forEach((seq, i) =>
                assert.equal(seq.length, chromosomes[i].length, `${name}: chromosome ${i}'s mapped length`));
            assert.ok(out.Coverage >= 0 && out.Coverage <= 1, `${name}: coverage ${out.Coverage}`);
            return {
                name,
                variants: out["Number of variants"],
                coverage: Number(out.Coverage.toFixed(6)),
                mapped: sha256(out["Mapped sequences"].join("\n")),
                vcf: sha256(out.VCF),
            };
        });
        return { reference: sha256(chromosomes.join("\n")), mapped };
    } finally {
        ska.free();
    }
}

// The alignment, as actions.ts and AlignmentExtractor.worker.ts drive it.
function alignment(api, inputs) {
    const runId = 1;
    messages.length = 0;
    const align = api.AlignData.new(31, runId);
    try {
        const plan = JSON.parse(align.plan_samples(inputs.map((f) => f.name), runId));
        for (const { sampleIndex, name, firstFileIndex, secondFileIndex } of plan) {
            const second = secondFileIndex == null ? null : inputs[secondFileIndex];
            const words = api.extract_alignment_sample(inputs[firstFileIndex], second, 31, sampleIndex, name,
                true, 1, 0, 20, 2);
            assert.equal(words.length % 2, 0, `${name}: packed words come in pairs at k=31`);
            align.add_sample_words(sampleIndex, name, words);
        }
        const finished = JSON.parse(align.finish_alignment());
        const csv = gunzipSync(Buffer.from(align.get_distances_csv_gzip())).toString();
        // A header row and a name column around the square matrix of distances.
        const matrix = csv.trim().split("\n").slice(1).map((row) => row.split(",").slice(1).map(Number));
        assert.equal(matrix.length, plan.length, "one distance row per sample");
        matrix.forEach((row, i) => row.forEach((d, j) => assert.equal(d, matrix[j][i], `distance ${i},${j}`)));
        assert.ok(matrix.every((row, i) => row[i] === 0), "zero diagonal");
        align.begin_alignment_export(plan.length);
        const chunks = [];
        while (!align.alignment_export_complete()) chunks.push(gunzipSync(Buffer.from(align.next_alignment_export_chunk())));
        align.finish_alignment_export();
        const fasta = Buffer.concat(chunks).toString();
        const seqs = fastaRecords(fasta);
        assert.equal(seqs.length, plan.length, "one aligned sequence per sample");
        assert.ok(seqs.every((s) => s.length === seqs[0].length), "aligned sequences of one length");
        const stages = messages.filter((m) => m?.alignmentProgress?.runId === runId).map((m) => m.alignmentProgress.stage);
        for (const stage of ["preparing", "distances"]) assert.ok(stages.includes(stage), `no "${stage}" progress`);
        return {
            samples: plan.length,
            tips: [...finished.treeTipNames].sort(),
            length: seqs[0].length,
            alignment: sha256(fasta),
            distances: sha256(csv),
        };
    } finally {
        align.free();
    }
}

for (const directory of values.package) {
    const { api, root, target } = await loadPackage(directory);
    const reference = fileFrom(join(files, "test_ref.fa"));
    golden(`ska-${target}`, {
        fasta: mapping(api, reference, [["test_1", fileFrom(join(files, "test_1.fa")), null],
            ["test_2", fileFrom(join(files, "test_2.fa")), null]], 0),
        // ska's own FASTQ files hold 10 bp reads, shorter than k=31, so the assembler's E. coli fixture
        // stands in: 30x reads of a 20 kb window, mapped to that window.
        fastq: mapping(api, fileFrom(join(ecoli, "ecoli_window.fasta")), [["ecoli_window",
            fileFrom(join(ecoli, "ecoli_window_1.fastq.gz")), fileFrom(join(ecoli, "ecoli_window_2.fastq.gz"))]], 5),
        alignment: alignment(api, [1, 2, 3, 4, 5].map((i) => fileFrom(join(files, "par_test", `test_${i}.fa.gz`)))),
    }, values["update-golden"]);
    report(`PASS ska ${target} (${root})`);
}
