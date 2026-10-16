// Shared by the WASM checks: stand-ins for the browser APIs the site's packages use, package loading,
// and the golden files. Only Node built-ins are used.
import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { existsSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { basename, dirname, join, resolve, sep } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { parseArgs } from "node:util";

/** The sparrowhawk repository, and the checks' own data. */
export const REPO = resolve(dirname(fileURLToPath(import.meta.url)), "..", "..");
export const DATA = join(REPO, "tests", "wasm", "data");

/** A browser File: the glue code reads it through slice() and FileReaderSync, the workers through arrayBuffer(). */
export class LocalFile {
    constructor(name, bytes) { this.name = name; this.bytes = Buffer.from(bytes); this.size = this.bytes.length; }
    slice(start = 0, end = this.size) { return new LocalFile(this.name, this.bytes.subarray(start, end)); }
    async arrayBuffer() { return Uint8Array.from(this.bytes).buffer; }
}
globalThis.FileReaderSync = class {
    readAsArrayBuffer(blob) { return Uint8Array.from(blob.bytes).buffer; }
};
/** Everything the packages post, as a worker's postMessage would carry it to the page. */
export const messages = [];
globalThis.postMessage = (message) => { messages.push(message); };

export const fileFrom = (path, name = basename(path)) => new LocalFile(name, readFileSync(path));
export const sha256 = (data) => createHash("sha256").update(data).digest("hex");
export const report = (text) => process.stdout.write(`${text}\n`);
export const fastaRecords = (text) => text.split(">").slice(1).map((record) => record.split("\n").slice(1).join(""));

const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A" };
export const reverseComplement = (s) => [...s].reverse().map((b) => COMPLEMENT[b] ?? b).join("");
export const canonical = (s) => { const r = reverseComplement(s); return s < r ? s : r; };

/** Common options, plus a check's own. The packages' console output is dropped unless --verbose. */
export function options(extra = {}) {
    const { values } = parseArgs({
        options: {
            package: { type: "string", multiple: true },
            "update-golden": { type: "boolean", default: false },
            verbose: { type: "boolean", default: false },
            ...extra,
        },
    });
    assert.ok(values.package?.length, "Pass --package DIR at least once (a web-target package).");
    if (!values.verbose) for (const level of ["log", "info", "warn", "debug"]) console[level] = () => {};
    return values;
}

/** A web-target package, its memory, and whether it is build_wasm64.py's wasm64 build. */
export async function loadPackage(directory) {
    const root = resolve(directory);
    const api = await import(pathToFileURL(join(root, "index.js")).href);
    assert.equal(typeof api.default, "function", `${root}: use a web-target package, not the bundler one.`);
    const exports = await api.default({ module_or_path: readFileSync(join(root, "index_bg.wasm")) });
    const target = root.includes(`${sep}pkg_wasm64${sep}`) ? "wasm64" : "wasm32";
    return { api, root, target, memory: exports.memory };
}

/** Compares `observed` with tests/wasm/golden/<name>.json, or records it there on a first run or with --update-golden. */
export function golden(name, observed, update) {
    const path = join(REPO, "tests", "wasm", "golden", `${name}.json`);
    if (update || !existsSync(path)) {
        mkdirSync(dirname(path), { recursive: true });
        writeFileSync(path, `${JSON.stringify(observed, null, 2)}\n`);
        report(`RECORDED ${path}`);
        return;
    }
    assert.deepEqual(observed, JSON.parse(readFileSync(path, "utf8")),
        `${name} differs from ${path}; rerun with --update-golden if the change is intended`);
}
