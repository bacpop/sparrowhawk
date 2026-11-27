---
name: rust-wasm-vue-tool
description: Building a browser tool from a Rust algorithm — the wasm interface, worker and Vue page. Use when adapting an existing method to WebAssembly, wrapping an existing Rust crate for the browser, or writing a new method directly for it.
---

# SKILLS.md — building a Rust → WebAssembly → Vue tool

A complete playbook for a computational tool that runs entirely in the browser: a
Rust algorithm compiled to WebAssembly, driven from a Web Worker, presented as a
Vue page.

**This document is self-contained.** Everything needed to act is here. Parts
I–III work for any project; Part IV adds what the sparrowhawk repository
expects when the tool goes into it rather than into a site of its own.

Four parts. **Part I** gets a correct Rust library, natively, with no wasm and no
Vue. **Part II** puts it in a browser. **Part III** is reference — long code
blocks and tables, consulted on demand, kept out of the way of the steps.
**Part IV** is the checklist for adding the tool to sparrowhawk itself.

**§2 is the skill.** Once you can tell a right answer from a wrong one, everything
downstream is mechanical. Spend disproportionate effort there.

Two worked examples, both real:

- [combine-web](https://github.com/vrbouza/combine-web) — a re-implementation of the maximum-likelihood fit from
  [CMS Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/),
  a C++/ROOT tool. Correctness meant agreeing numerically, to ~1e-6.
- [rammap-web](https://github.com/vrbouza/rammap-web) — a bridge around [rammap](https://github.com/jwanglab/rammap),
  a Rust minimap2-compatible mapper. Correctness meant byte-identical output
  against its own CLI.

The first is particle physics, the second bioinformatics. Nothing about this
architecture is domain-specific.

---

# Part I — get a correct Rust library

Native only. Writing wasm bindings or Vue components before this part is finished
means debugging through a browser round trip instead of a single command, unable
to tell whether a wrong number came from the algorithm or the plumbing.

## 1. Scope and feasibility

**Scope honestly, and refuse the rest by name.** Combine is ~41k lines of C++ on
ROOT and RooFit. Re-implementing it is not a project; re-implementing *the binned
likelihood fit for a single parameter of interest* is. Write down what is in and
what is out, and make the tool **reject out-of-scope input with a message naming
the construct and where it appeared**, rather than silently mis-computing.
combine-web's parser reads 166 real datacards and rejects 46, each by name.

**Check it can reach the browser at all**, before promising anything:

```bash
cargo tree --target wasm32-unknown-unknown 2>&1 | grep -E "error|could not"
```

Crates needing file I/O, threads, system clocks, native dynamic linking, or C
libraries via `cc`/`cmake` will not compile. `needletail` is a known example.
Three ways out: swap for a pure-Rust alternative, re-implement the subset you
need, or gate it out of the wasm build (Part III, *Cargo.toml*).

Some things cannot go at all: ROOT cannot be compiled to wasm — its `libCling`
JIT is architecturally incompatible with the sandbox, and `libCore` depends on
it. That ruled out wrapping Combine and forced a re-implementation.

**When upstream genuinely will not compile:** say so, and present the choice —
re-implement the subset, or find a different upstream. Proceeding on the hope
that the dependency tree resolves itself wastes the whole of Part I.

**Done when:**

- [ ] The in-scope set is written down, and every out-of-scope construct has a
      rejection path that names it.
- [ ] `cargo tree --target wasm32-unknown-unknown` is clean, or every offender
      has a decided remedy.

## 2. Where correctness comes from

**This is the skill.** It is the one question whose answer changes what you do,
and it decides what "done" means for everything after it. Answer it before
writing code.

The answer is always the same shape: **an oracle** — one command you can run that
says whether the output is right.

### If a reference implementation exists

It is your oracle. Install it and diff against it — byte for byte where the
output is text, number for number where it is numeric.

**Control for the legitimate divergences first, or you will chase ghosts:**

- **Thread count.** Run the reference single-threaded (`-t 1`); parallel
  reductions reorder floating-point work.
- **Lines that differ by construction.** A SAM `@PG CL:` line records the command
  line, and a browser has no argv. Exclude it rather than trying to match it.
- **Defaults the wasm path hardcodes differently** from the CLI. rammap's demo
  capped `index_max_occ` at 50 000 where the CLI leaves it uncapped.
- **Optional-field flags.** minimap2's `-y` comment copying adds a field to every
  record. Emitting it unconditionally is a real bug — and diffing found it.

**Generate the reference for the exact input you are running**, never the nearest
published number. In combine-web, comparing against a published result for a
*similar* input looked like a 14% discrepancy and cost real time; the two inputs
differed by one line.

Make the oracle **tight** — fast and deterministic enough to run after every
change — and script it, so a regression turns it **red** without anyone deciding
to look.

#### If upstream is already Rust

Two things change.

**Prefer a bridge crate to a fork.** Depend on upstream by git `rev` in
`Cargo.toml`, leave it untouched, and skip submodules entirely — there is then
nothing to keep in step. (sparrowhawk does both: its bridges depend on upstream
by git branch, pinned by a committed `Cargo.lock`, while forks it maintains, such
as ska and sketchlib, are submodules. A `rev` is the clearer pin.)

**Check first whether upstream already compiles to wasm.** rammap already had a
working wasm session and a demo page. The real job was not "port it" but "close
the gaps that make it unusable": no way to save output, invalid SAM with no
`@HD`/`@SQ`/`@PG`, and quality strings dropped. Assuming a greenfield port
rebuilds what exists and misses what matters.

**Be ready to drive upstream's internals rather than its public API.** In rammap,
every high-level constructor took either a path (there is no filesystem in a
browser) or a fully materialised `Vec<(String, Vec<u8>)>` — which would hold the
reference as ASCII alongside the packed index, exactly what the streaming index
builder exists to avoid. The one method that bridged the gap was private. Going a
level down was the only route, and also the only way to get SAM output at all.

### If no reference exists

Construct the oracle. In descending order of value:

1. **A slow oracle you write yourself.** An obviously-correct implementation,
   however inefficient — O(n²) is fine — differential-tested against the real one
   on small inputs. This is the single most effective technique for a new method.
2. **Closed-form cases**, where the answer follows from the mathematics.
   combine-web's sharpest test is exactly this: profiling a nuisance the data
   cannot constrain must reproduce its own prior — a parabola crossing 68% at
   exactly ±1. It generalises to any method with an analytically tractable limit.
3. **Invariants and metamorphic relations.** Properties that hold for every
   input: conservation, monotonicity, symmetry, idempotence — or knowing that
   `f(2x) = 2f(x)` even when `f(x)` is unknown. `proptest` fits here.

**An oracle derived from your own implementation proves nothing.** Writing the
test by reading the code encodes the same misunderstanding twice. Derive the
analytic cases from the specification.

**When no oracle can be constructed:** stop and say so. List what you tried, and
ask for a specification, a worked example with expected output, or a domain
expert who can adjudicate. Building on an undefined notion of correct produces
something nobody can ever validate.

**Done when:**

- [ ] You can **name one command** that decides whether the output is right —
      a diff script, or the oracle suite — and have **already run it** at least
      once, showing its invocation and output.
- [ ] It is **tight**: seconds, deterministic, runnable unattended.
- [ ] It goes **red** on a wrong answer, not merely green on a clean exit.

**No command that decides correct, no implementation.**

## 3. Build it headless

**Phase one is a library plus a CLI — or, for a bridge, a Node harness — with no
UI at all.** Ship that phase alone rather than bundling a bridge crate, four
features and a front end into one push. It answers the questions that can
invalidate everything downstream; if it fails, nothing else is wasted.

Give the binary **one subcommand per capability** — combine-web has `parse`,
`yields`, `fit`, `scan`, `impacts`. That is what makes the oracle scriptable and
gives each capability a natural test.

**Repo layout.** For a new project, the root holds `rust/` and `www/` and little
else but config (sparrowhawk itself also has `docs/`, `electron/`, `scripts/` and
`tests/`; Part IV):

```
rust/<tool>-bridge/          all #[wasm_bindgen] lives here
├── src/
├── scripts/                 build, compare-native, test harnesses
└── tests/data/              fixtures, or generators for any too large to commit
www/                         the Vue app
```

Scripts and test data for one crate belong **inside the crate**, mirroring
`rust/esm-bridge/`. Checks that span every crate — sparrowhawk's `scripts/` and
`tests/wasm/` — are the exception. Commit lockfiles: `rust/*/Cargo.lock` and
`www/package-lock.json`.

**When the headless phase cannot be made to reproduce:** stop there. A wasm build
will not fix a wrong answer, and every later layer makes the wrong answer harder
to see. Report which oracle case fails and by how much.

**Done when:**

- [ ] The CLI or Node harness runs the oracle from §2 and it comes back green.
- [ ] `git ls-files` shows no wasm bindings and no Vue anywhere in the tree.
- [ ] Every capability you intend to expose has its own subcommand.

## 4. Design for the target now

Anything you are *writing* — a re-implementation as much as a new method — can be
shaped by the browser's constraints from the first commit. Retrofitting them is
far more work:

- **Take bytes, never paths.** There is no filesystem.
- **Bound the working set**, and know what scales with input size versus
  reference size. The ceiling is 4 GB on wasm32 (16 GiB on wasm64, at a cost in
  build complexity; §6).
- **Design for chunk-level parallelism across workers**, which is what a pool
  gives you; `rayon` does not run here.
- **Structure work as a loop over records or chunks**, not one monolithic call.
  That is what makes progress reporting and cancellation possible at all — a
  single long call cannot be preempted (§8).
- **Emit output incrementally** rather than accumulating it.

**A bridge is the exception.** It inherits assumptions made years ago and must
work around them, which is what driving upstream's internals is about.

## 5. Record where the precision goes

Keep a running note — a `PRECISION.md` alongside the code works well — of every
place precision is lost or the reference disagrees with itself. Two findings from
combine-web, and why the habit earns its keep:

- Combine quantises a shape-derived constant through a `"%f"` format and
  downcasts templates to single precision during rebinning. **Reproducing those
  two losses deliberately** took agreement from ~1e-3 to ~1e-6. Sometimes
  matching the reference means copying its flaws.
- On one input the reference's own crossing search collapsed and it reported a
  spurious zero, while our answer was right. An oracle must detect that the
  *reference* is degenerate and refuse to treat it as truth.

**Measure performance claims; state only what you have run.** A predicted 3–5×
speedup from skipping an O(n²) step measured at **1%**. The arithmetic was right
and the conclusion was wrong.

**No green headless run, no wasm and no Vue.**

---

# Part II — put it in the browser

From here nothing depends on how Part I was reached.

## 6. The wasm interface

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack --version 0.15.0 --locked
```

**Pin wasm-pack**, locally and in CI. Older versions bundle a `wasm-opt` that
rejects current wasm-bindgen output ("Only 1 table definition allowed in MVP"),
and some CI actions silently fall back to 0.9.1 when left unpinned.

`wasm-pack` calls cargo, runs `wasm-bindgen` to generate the JS/TS glue,
optionally runs `wasm-opt`, and emits a `pkg/` with the `.wasm`, a JS module and
a `.d.ts`.

Use the **default `bundler` target** inside webpack: with
`experiments.asyncWebAssembly`, a dynamic `import()` instantiates the module, so
**there is no `init()` call**. A `--target web` package is different: it exports
a default `init` that must be awaited, which is how Node harnesses and wasm64
packages are loaded. `--out-dir` works fine (wasm-pack 0.15), despite older
documentation warning it is forwarded to cargo as the nightly-only
`--artifact-dir` — useful when you want two packages, a `web` or `nodejs` build
in a folder of its own for a Node harness, and the bundler build for the app.
Never point the harness build at the app's own `pkg*` folder.

**Start with `wasm32-unknown-unknown`, and plan to fit in 4 GB.** wasm64
(Memory64) works and lifts the ceiling to 16 GiB, but it is a second build path:
wasm-bindgen 0.2.120 or newer, nightly Rust with `rust-src` and
`-Z build-std=std,panic_abort`, `wasm-bindgen-cli` instead of wasm-pack, a
`--target web` package, and a loader that falls back to wasm32 where wasm64
fails. Dependencies can break on it: under simd128, `fixedbitset` imports
`core::arch::wasm32`, which does not exist on wasm64. sparrowhawk does all this
for its assembler and ska, in `scripts/build_wasm64.py` and
`www/src/wasm/loader.ts`; copy those rather than starting again.

Gate wasm-only code on **`target_family = "wasm"`** — the predicate this platform
uses, not a cargo feature, and not `target_arch = "wasm32"`, which silently drops
the code from a wasm64 build:

```rust
#[cfg(target_family = "wasm")]
#[wasm_bindgen]
pub fn init_panic_hook() { console_error_panic_hook::set_once(); }
```

See Part III, ***Cargo.toml and `.cargo/config.toml`*** — including the trap that
a dependency's `.cargo/config.toml` is never inherited, so SIMD flags silently
vanish and the stack stays at its 1 MB default. It also applies only when cargo
runs from the crate's folder, as wasm-pack does; a build driven from elsewhere,
such as a wasm64 script, must pass the flags itself.

### The interface pattern

Expose **one stateful handle**: a struct constructed once, mutated by work
methods, drained by getters. Real constructor, `&[u8]` in, `Result<_, JsValue>`
out. Full example in Part III, ***The wasm handle***.

**Return `Result`, not `panic!`** — and the reason is stronger than tidiness: **a
Rust panic poisons the wasm module.** The instance cannot be reused, so the
worker must be discarded and rebuilt rather than retried. Keep
`console_error_panic_hook` for genuine bugs; return `Result` for everything you
anticipate.

**Guardrail: exactly one `#[wasm_bindgen(start)]` per module.** If upstream has
wasm-bindgen exports gated on the target with no feature to opt out, its exports
link in alongside yours — harmless, and mildly useful since its panic hook
installs itself via `__wbindgen_start()`. Declaring a second start function fails
the build with `cannot specify two 'start' functions`. Keep every exported name
distinct from upstream's too, or you get two `export class X` in one ES module.

Three idioms from the existing bridges: **`Arc<Inner>`** so several sessions
share one large index; **`std::mem::take`** in `take_*()` methods, so handing
bytes to JS frees the wasm-side copy immediately; and **`js_sys::Object` +
`Reflect::set`** for small structured returns.

**Design the interface so work can be sliced.** combine-web's expensive stage
takes an anchor plus a list of indices, so any subset computes independently.
That one decision made a worker pool possible later with no restructuring.

### Getting data across

`&[u8]`, `&str` and numbers are copied at the boundary. That copy is cheap
relative to any real computation, so design for the clearest interface and let
the copy happen.

**Return `Vec<u8>` — a JS `Uint8Array` — for bulk output.** The demo rammap-web
replaced did `lastOutput += chunk` into a JS string: quadratic, and a hard wall
at V8's ~512 MB cap. Returning bytes makes it natural for the caller to stream
them into a sink. rammap-web produced 84 MB of SAM in-browser this way, and
253 MB through the Node path with peak RSS of 219 MB — more output than memory,
which is the point.

For inputs too large to hold, stream them:

```rust
pub fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<u8>, JsValue> { … }
pub fn finish(&mut self)                   -> Result<Vec<u8>, JsValue> { … }
```

**For structured results**, `serde-wasm-bindgen` with `#[derive(Serialize)]` and
`#[serde(rename_all = "camelCase")]` gives native JS objects mirroring your
TypeScript interfaces one-for-one. No sparrowhawk crate does this yet, so there
is no in-repo example to copy. Its older crates return JSON strings that JS
re-parses, which works but is untyped on both sides; the newer ones return typed
arrays (`Vec<f32>` arrives as a `Float32Array`) or small objects built with
`js_sys`.

**Optimise for speed and measure the result.** A profile of `opt-level = 3`,
`lto = "fat"`, `strip = "debuginfo"` produced **0.37 MB** in rammap-web, against
the 2.5–4 MB that had been predicted. For compute-bound tools the speed is free.

### Progress, reported from Rust

Bind `postMessage` and call it at named lifecycle points from inside the worker —
`post_state("preprocess:start")`, `"preprocess:end"`, `"assembly:saving"`. Full
example in Part III, ***Progress from Rust***.

**Done when:**

- [ ] `cargo build --lib --target wasm32-unknown-unknown` is green.
- [ ] One exported call round-trips through a Node harness and returns the same
      answer the CLI does.
- [ ] Every anticipated failure returns `Err`, and the panic hook is reserved for
      bugs.

## 7. Wiring the build

One `WasmPackPlugin` instance per crate, `experiments.asyncWebAssembly` on, and
`worker-loader` handling `*.worker.{js,ts}`. The complete block is in Part III,
***vue.config.js***; two settings in it are traps worth naming here.

**Decide where the site is served from.** If the app loads anything by a
root-relative URL — large data files in `public/`, a wasm64 package — it can only
run from a site root: use `publicPath: "/"` and skip the next two traps, as
sparrowhawk does. Otherwise **`publicPath` must differ between dev and
production.** A relative `"./"` makes a production `dist/` work from any
subdirectory, but Vue CLI passes the value verbatim as the dev server's URL
pathname, so `npm run serve` yields malformed asset URLs. Use
`process.env.NODE_ENV === "production" ? "./" : "/"`.

**A relative `publicPath` with the default `js` assets directory makes workers
404 on `js/js/…`.** The document resolves chunks against `dist/`, but a worker
resolves them against its own location in `dist/js/`, and the chunk filename
already starts with `js/`. Cure with flat output: `assetsDir: ""` plus explicit
`output.filename`/`chunkFilename`.

**If the same worker is built twice** — once against a wasm32 package and once
against a wasm64 one — give `worker-loader` hashed names
(`filename: "js/[name].[contenthash].worker.js"`), or the two files collide.

**There is no hot reload, deliberately.** Restart `npm run serve` after any
change, front end included. First start is slow — every crate builds in release
mode.

**Ignore the generated packages with a glob that works.** `www/src/pkg` ignores
only `pkg` and silently misses `pkg_ska`, `pkg_deacon` and the rest. Use
`www/src/pkg*/` — the trailing slash keeps it to directories, so a source file
merely starting with `pkg` survives.

**Done when:**

- [ ] `npm run build` produces a `dist/` that loads where it will be deployed —
      a subdirectory, if you chose a relative `publicPath`.
- [ ] `npm run serve` boots and the app runs — the command developers actually
      type, and the one usually left untested.
- [ ] `git status` is clean after a build: every generated `pkg*/` is ignored.

## 8. The worker, and a pool

Two files per tool: a thin `*.worker.ts` dispatcher keyed on boolean flags, and a
driver class holding the lazily-imported wasm module and the handle. Both in
Part III, ***The worker, driver and pool***.

### Three bugs that will not announce themselves

All three were found building rammap-web; the snippets are from that project.

**An `async` `onmessage` handler does not serialise messages.** Two
`postMessage` calls back to back start concurrently, so a page posting
`{buildIndex}` then `{map}` runs the map before an index exists. Fix with an
explicit promise queue, and let `cancel` bypass it so it reaches the job already
running:

```js
let queue = Promise.resolve();
ctx.onmessage = (evt) => {
    if (evt.data.cancel) { aligner.cancel(); return; }   // deliberately outside
    queue = queue.then(() => handle(evt.data)).catch(reportError);
};
```

**`for await` over a blob-backed stream starves the macrotask queue.** Its reads
resolve through microtasks, which run to exhaustion before the next task — so a
`{cancel: true}` message is not delivered until the whole file has already been
processed. Cancellation appears to work while doing nothing. Yield explicitly on
each chunk:

```js
await new Promise((r) => setTimeout(r, 0));
```

**Name progress fields `completed`.** Consumers test completion with
`if (msg.done)`, so a progress message carrying `done: <count>` reads as truthy
the moment any work happens and progress is mistaken for a finished run.
(sparrowhawk's `esm-bridge` predates this rule and still sends `done` as a
count; do not copy it.)

### Cancellation, and its real limit

The constraint is structural: one wasm instance lives in one worker, its memory
is a plain `ArrayBuffer` rather than shared, and every export is synchronous — so
the main thread cannot signal into a running call, and wasm has no preemption.
Three tiers are achievable:

1. **Between chunks** — a flag checked before each call. Adapt chunk size at
   runtime towards ~250 ms so cancel latency is bounded on any input; fixed large
   chunks mean seconds on short records.
2. **Within a chunk** — a Rust-side `Cell<bool>` checked per record, so a chunk
   already dispatched exits at its next record boundary.
3. **Force stop** — `worker.terminate()`, the only true mid-chunk abort. It
   destroys the instance, so any cached index is lost and must be rebuilt. Offer
   it once a soft cancel has been pending for a few seconds.

**Preempting a single long call is not possible.** One pathological record can
hold the worker for seconds, and only `terminate()` touches it.

### Pools

If the work is embarrassingly parallel, spawn a pool sized from a slider in the
page, terminating the old one first (Part III). Hand out **small chunks and
dispatch the next as each reply lands** — work stealing — so one slow item cannot
stall a whole static slice.

**Ordinary workers, each with its own wasm instance — not wasm threads.** Shared
memory needs `SharedArrayBuffer` and therefore COOP/COEP headers, which break
plain static hosting. Each worker recompiles the module and holds its own copy of
the data, so cap the pool: sparrowhawk's `clampWorkerCount()` caps it at
`navigator.hardwareConcurrency` (12 when unknown), and memory is often the
tighter limit.

**Expect sub-linear scaling and measure it.** combine-web's parallel stage:
5.3 s at one worker, 3.3 s at two, 2.8 s at four, 1.9 s at eight — 2.8×, not 8×,
because a laptop's "12 cores" are mostly efficiency cores. Measure before
attributing a result to chunk size: tuning it there changed nothing.

**Done when:**

- [ ] A cancel issued mid-run is honoured within roughly 250 ms.
- [ ] Progress messages name the current phase, and the field is `completed`.
- [ ] Two jobs posted back to back run in order.
- [ ] Pool scaling is measured at one worker and at the cap, and the numbers are
      written down.

## 9. The page

**Parameters in a 350 px left rail; drop zone and results to the right.** The
layout, the exact class strings and the full conventions table are in Part III,
***The page***. The rules that decide whether it looks right:

- Every control is a label row with an `Info` icon at exactly
  `w-3.5 h-3.5 text-gray-400 cursor-help`.
- A **threads slider** wherever work is parallelised, wired to the pool size.
- Export buttons in the **action bar at the top** — results run to several
  screens and nobody scrolls back.
- Results **side by side** on wide screens; widen the page if a split would
  shrink a plot.
- Progress that names what is happening (`Impacts 37/72 across 8 workers`).
- Warn **in the interface**, not the console, when a result is approximate. A
  value that could not be computed must never read as a real zero.

**Saving large results.** Bytes must never accumulate in JavaScript. Keep a
bounded preview — 2 MB or 5000 lines — for display, and write the rest to a sink.
Two backends: the **File System Access API**, where `showSaveFilePicker()` runs
on the main thread under user activation and the handle is passed into the
worker, so `await writable.write()` applies real backpressure to disk (Chromium
desktop only); and **blob parts**, coalesced into ~8 MB pieces and assembled at
the end, which Firefox and Safari take, and which escapes the string cap because
blob parts spill to disk.

In sparrowhawk, all of this sits behind `www/src/platform/files.ts`
(`openBinaryFileSink`, `saveTextFile`, `saveBinaryFile`), which opens the sink on
the main thread and also routes saves through native dialogs in the desktop
program. Use it rather than the browser APIs directly.

**Plotting** is hand-rolled declarative SVG for anything numeric and precise —
`computed` scales, module-level layout constants, no imperative drawing. Reach
for plotly or d3 only for genuinely interactive views. Two traps:

- **Check the *computed* font size, not the number in the `viewBox`.** An SVG
  with `width: 100%` and a `viewBox` wider than its column is scaled down: 9 px
  labels in a 920-unit plot rendered into a 656 px column reach the eye at
  **6.4 px**. Measure with `getComputedStyle`, in a browser, at the real width.
- **A scoped `svg { … }` rule compiles to `svg[data-v-…]`**, specificity 0,1,1 —
  which **outranks Tailwind's `w-3.5`** at 0,1,0. Put a class on your chart
  (`svg.chart { … }`) or the next icon inside that component inherits the chart's
  full-width bordered styling and renders enormous.

**Done when:**

- [ ] The smallest plot label measures **≥10 CSS px** via `getComputedStyle`, in
      a browser, at the real column width.
- [ ] Every parameter has an `Info` tooltip; parallel work has a threads slider.
- [ ] Exports sit above the results, and an approximate value is marked in the
      interface.

## 10. Verifying you are done

- [ ] The oracle from §2 is green on every input you claim support for, and the
      agreement level is written down.
- [ ] `cargo test` covers the algorithm, not just the parser.
- [ ] `cargo build --lib --target wasm32-unknown-unknown` is green. Note `--lib`:
      the CLI binary usually depends on crates that are not wasm-compatible.
- [ ] `npm run build` **and** `npm run serve` both work.
- [ ] The browser produces the same numbers as the CLI, checked in a browser.

**Build a Node harness that mirrors the worker** and diff it against the native
CLI byte for byte. It is cheap, and it is what catches real bugs.
sparrowhawk's is `tests/wasm/` (Part IV): a small shim stands in for `File`,
`FileReaderSync` and `postMessage`, and outputs are compared with committed
golden files.

**Headless browser testing needs no puppeteer.** Node 22+ ships a `WebSocket`
client, so about sixty lines drives Chrome over the DevTools Protocol. Drive the
app **through its own store and actions** rather than a parallel test page — the
tests then exercise what a user drives, and there is no second front end to keep
in sync.

Three ways that harness bites, all of which read as code bugs and are not:

- **Chrome's blob quota is sized from the filesystem backing its profile.** On a
  root filesystem reporting full, fetching a ~77 MB fixture into a Blob fails
  with a bare `net::ERR_FAILED` while smaller fetches succeed. Put the profile on
  tmpfs (`/dev/shm`). Passing `--disable-dev-shm-usage` redirects the blob store
  onto the full disk and causes the same failure — use it only when `/dev/shm` is
  genuinely tiny.
- **Killing a dev server started by `npm run serve` needs a process-group kill**
  — `detached: true` at spawn, then `process.kill(-pid)` — because npm spawns a
  shell which spawns the real server.
- **Pre-flight the port.** `vue-cli-service` silently falls back to the next free
  one, so a leaked server means the test polls one process while driving another,
  and the symptom is an unexplained boot timeout. Fail loudly rather than
  attaching to a server you did not start.

**Guardrail: identify a process before killing it.** `readlink /proc/<pid>/cwd`,
as a *separate step*, then kill that pid. Combining discovery and `pkill` in one
command has killed unrelated dev servers belonging to other projects.

---

# Part III — Reference

Consulted on demand. Nothing here is a step.

## Cargo.toml and `.cargo/config.toml`

```toml
[lib]
crate-type = ["cdylib", "rlib"]     # cdylib for wasm, rlib so the CLI can link it

[dependencies]
wasm-bindgen = "0.2"
js-sys = "0.3"
web-sys = { version = "0.3", features = ["console"] }
console_error_panic_hook = "0.1"

# Native-only: keeps rayon, clap and indicatif out of the wasm build entirely.
[target.'cfg(not(target_family = "wasm"))'.dependencies]
clap = { version = "4", features = ["derive"] }
rayon = "1"
indicatif = "0.17"
```

**`.cargo/config.toml` is never inherited from a dependency.** Cargo walks up
from the invocation directory, so upstream's flags apply only inside *its*
checkout. If upstream sets `-C target-feature=+simd128` and you do not, the v128
kernels compile out and the browser runs scalar fallbacks — no warning, no error,
just slow. Keep your own copy:

```toml
# .cargo/config.toml — scoped to the target so native `cargo test` is unaffected
[target.wasm32-unknown-unknown]
rustflags = [
    "-C", "target-feature=+simd128",
    "-C", "link-args=-zstack-size=8000000",   # the 1 MB default overflows on deep code
]
```

## The wasm handle

From `rust/sparrowhawk-amr/src/lib.rs`:

```rust
#[wasm_bindgen]
pub struct AmrDetector { index: AmrIndex }

#[wasm_bindgen]
impl AmrDetector {
    #[wasm_bindgen(constructor)]
    pub fn new(index_bytes: &[u8]) -> Result<AmrDetector, JsValue> {
        init_panic_hook();
        let index = load_index_from_bytes(index_bytes).map_err(js_error)?;
        if index.alphabet != IndexAlphabet::Dna {
            return Err(JsValue::from_str("AMR direct mode requires a DNA AMR index"));
        }
        Ok(Self { index })
    }

    pub fn info(&self) -> String { self.index.stats_string() }
}
```

## Progress from Rust

From `rust/sparrowhawk-asm/src/lib.rs`:

```rust
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_name = postMessage)]
    fn post_message(data: &JsValue);
}

#[cfg(target_family = "wasm")]
pub(crate) fn post_state(state: &str) {
    let obj = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&obj, &JsValue::from_str("assemblyState"),
                                 &JsValue::from_str(state));
    post_message(&obj.into());
}

// Natively a no-op, so shared code can call it on every target.
#[cfg(not(target_family = "wasm"))]
pub(crate) fn post_state(_state: &str) {}
```

## vue.config.js

```js
const path = require("path");
const { defineConfig } = require("@vue/cli-service");
const WasmPackPlugin = require("@wasm-tool/wasm-pack-plugin");

module.exports = defineConfig({
    // Relative for production so dist/ works from any subdirectory; absolute for
    // dev, because Vue CLI uses this verbatim as the dev server's URL pathname.
    publicPath: process.env.NODE_ENV === "production" ? "./" : "/",
    // Flat output, so a worker in dist/js/ does not resolve chunks to js/js/….
    assetsDir: "",
    configureWebpack: {
        experiments: { asyncWebAssembly: true },
        devServer: {
            // Otherwise every keystroke retriggers a release wasm build.
            hot: false, liveReload: false, watchFiles: { paths: [] },
        },
    },
    chainWebpack: (config) => {
        config.output
            .filename("[name].[contenthash:8].js")
            .chunkFilename("[name].[contenthash:8].js");

        config
            .plugin("wasm-pack_mytool")
            .use(WasmPackPlugin)
            .init((Plugin) => new Plugin({
                crateDirectory: path.resolve(__dirname, "../rust/mytool"),
                outDir:         path.resolve(__dirname, "./src/pkg_mytool"),
                forceMode:      "production",
            }))
            .end();

        // Worker entrypoints bypass the default rules.
        config.module.rule("js").exclude.add(/\.worker\.(js|ts)$/);
        config.module.rule("ts").exclude.add(/\.worker\.ts$/);
        config.module
            .rule("worker")
            .test(/\.worker\.(js|ts)$/)
            .use("worker-loader").loader("worker-loader").end()
            .use("ts-loader").loader("ts-loader").options({ transpileOnly: true }).end();
    },
});
```

This block is a **starting point** for a new project. It is not a copy of
sparrowhawk's `www/vue.config.js`, which differs: `publicPath: "/"` (it loads
assets from the site root), the default `js/` assets folder with hashed worker
names, an `outputDir` switch for the desktop build, and one plugin per crate,
some with `extraArgs`.

## The worker, driver and pool

`Depleter.worker.ts` — a thin dispatcher on boolean flags:

```ts
import { Depleter } from './Depleter';

interface LoadIndexMessage { loadIndex: true; file: File; }
interface FilterMessage    { filter: true; file: File; deplete: boolean; }
type WorkerMessage = LoadIndexMessage | FilterMessage | { reset: true };

const ctx: Worker = self as unknown as Worker;
const depleter = new Depleter(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (!(evt.data instanceof Object)) return;
    if ('loadIndex' in evt.data && evt.data.loadIndex) {
        depleter.loadIndex((evt.data as LoadIndexMessage).file);
    } else if ('filter' in evt.data && evt.data.filter) {
        const m = evt.data as FilterMessage;
        depleter.filterReads(m.file, m.deplete);
    } else if ('reset' in evt.data && evt.data.reset) {
        depleter.resetAll();
    }
};
```

`Depleter.ts` — the driver, importing the wasm lazily and keeping the handle:

```ts
import type { WasmIndex } from "@/pkg_deacon";

export class Depleter {
    worker: Worker;
    index: WasmIndex | null = null;

    constructor(worker: Worker) { this.worker = worker; }

    async loadIndex(file: File): Promise<void> {
        try {
            const wasm = await import("@/pkg_deacon");
            const buf = await file.arrayBuffer();
            this.index = new wasm.WasmIndex(new Uint8Array(buf));
            this.worker.postMessage({ indexLoaded: true, info: this.index.info() });
        } catch {
            this.worker.postMessage({ error: true, message: 'index' });
        }
    }
}
```

sparrowhawk's drivers also time each run (`elapsedMs`) and read the instance's
memory (`wasmMemoryBytes`, from an import of `@/pkg_<tool>/index_bg.wasm`); the
pages show both.

The pool, from `www/src/store/actions.ts`:

```ts
async initSketchlibWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
    const {commit, state} = context;
    const workerCount = clampWorkerCount(numWorkers);
    for (const worker of state.workerState.workers_sketchlib) {
        worker.terminate();
    }
    commit("clearIdentifyWorkerMemory");
    const pool: Worker[] = [];
    for (let i = 0; i < workerCount; i++) {
        pool.push(new WorkerSketcher());
    }
    commit("SET_WORKERS_SKETCHLIB", pool);
},
```

## The page

```html
<div class="flex flex-col gap-6 md:flex-row md:gap-0">
  <div class="w-full md:w-[350px] md:shrink-0">
    <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
      <ScanFace class="w-6 h-6" /> Taxonomic ID
    </h1>

    <TaxonomicIDHelpCollapsible />

    <TooltipProvider>
      <div class="flex flex-col gap-4">
        <div>
          <p class="flex items-center gap-1">
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">What this parameter controls.</p>
              </TooltipContent>
            </Tooltip>
            Min read quality
          </p>
          <div class="flex flex-row items-center w-full gap-2">
            <VueSlider class="flex-grow" v-model="min_qual" :lazy="true"
                       :min="0" :max="33" :interval="1" :disabled="isRunning" />
            <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
              {{ min_qual }}
            </span>
          </div>
        </div>
      </div>
    </TooltipProvider>
  </div>

  <div class="min-w-0 w-full flex-1 overflow-visible pt-0 md:overflow-hidden md:pt-12">
    <!-- drop zone, status rows, results -->
  </div>
</div>
```

| Convention | Why |
|---|---|
| Parameters in the 350 px rail; drop zone and results right | Every page does this; deviating reads as broken |
| A **threads slider** wherever work is parallelised, wired to the pool size | The platform's pools work this way and users look for it |
| Export buttons in the **action bar at the top** | Results run to several screens and nobody scrolls back |
| Results **side by side** on wide screens | And widen the page if a split would shrink a plot |
| A help collapsible per tool, with `Overview` and `Parameters` tabs | Say plainly in the overview if it is a re-implementation rather than a port |
| Progress naming what is happening | A spinner alone is indistinguishable from a hang |
| Warn in the interface when a result is approximate | A value you could not compute must never read as a real zero |

## Limits to design around

| Limit | Consequence |
|---|---|
| **4 GB linear memory** (wasm32) | Stream or chunk; cap input size in the UI; wasm64 lifts it to 16 GiB (§6) |
| **Shallow call stack** | Rewrite deep recursion iteratively, or raise it via `.cargo/config.toml` |
| **No filesystem** | All input arrives from JS as bytes or strings |
| **No threads by default** | `rayon` does not run; parallelise with workers |
| **SIMD128 ≠ AVX/NEON** | x86 intrinsics will not compile; use `std::simd` or `wide` |
| **Incompatible crates** | Audit with `cargo tree --target wasm32-unknown-unknown` |

Expect 60–90% of native performance for compute-bound code, and worse for
allocation-heavy code or anything crossing the JS boundary often.

## House style

- **British English, Oxford comma**, sentence case for headings.
- **Comments: two or three lines at most per block, explaining why**, not what.
- **Stable Rust, explicit error handling, and a small dependency set.** Check
  what a dependency costs in wasm before adding it — `regex` is ~1.5 MB.
- **Pin a dependency; write down the trigger that would justify a fork** and wait
  for it.
- **Prefer the simplest correct thing plus a warning** over a clever conversion.
- **`uv` or `pixi`** for Python tooling.

---

# Part IV — adding the tool to sparrowhawk itself

Only when the tool goes into the [sparrowhawk](https://github.com/bacpop/sparrowhawk)
repository rather than a site of its own. Parts I–III still apply in full; this
is what the repository additionally expects, roughly in the order it comes up.
Copy from the newest tools — AMR detection (`sparrowhawk-amr`), Embed proteins
(`esm-bridge`) and Host depletion (`deacon-bridge`) — and read `AGENTS.md` for a
map of the repository.

## The crate

- [ ] **Where it lives.** A crate maintained in its own repository goes in as a
      git submodule under `rust/` (`.gitmodules`, cloned with
      `--recurse-submodules`). A bridge around an upstream library goes in-tree
      as `rust/<tool>-bridge/`, with `rust/<tool>-bridge/target/*` added to
      `.gitignore`. Commit its `Cargo.lock`.
- [ ] **`.cargo/config.toml`** with `+simd128` for `wasm32-unknown-unknown`,
      plus the stack size if the code recurses deeply (Part III).
- [ ] **It also builds natively**, unless it is wasm-only by construction like
      `deacon-bridge`: CI lints and tests in-tree crates natively.

## The build

- [ ] **`www/vue.config.js`**: one more `WasmPackPlugin`, with
      `outDir: path.resolve(__dirname, "./src/pkg_<tool>")`,
      `forceMode: "production"`, and `extraArgs` only if the crate needs
      features.
- [ ] **`www/.eslintignore`**: add `src/pkg_<tool>`. It lists every generated
      package by name, while `.gitignore`'s `www/src/pkg*/` glob already covers
      it.
- [ ] **Large data** (indexes, weights) goes in `www/public/` and is loaded with
      `loadAssetBytes()` or `loadAssetBlob()` from `www/src/platform/files.ts`,
      by root-relative name, so the same code works on the website and in the
      desktop program.

## The front end

- [ ] **Worker and driver** in `www/src/workers/` (`<Tool>.worker.ts` and
      `<Tool>.ts`), plus a `declare module '@/workers/<Tool>.worker'` entry in
      `www/src/workers/workers.d.ts`; without it the import fails type-checking.
- [ ] **Store.** In `www/src/types.ts`, the `AllResults<Tool>` type and the
      worker's field in `WorkerState`. In `RootState` (`store/state.ts`),
      `allResults_<tool>`, **and** its initial value in `emptyState()` in
      `www/src/utils.ts`. Mutations for the worker (`SET_WORKER(S)_<TOOL>`), the
      results and their reset; the run action and, for a pool, an
      `init<Tool>Workers` action modelled on `initSketchlibWorkers`; getters if
      the page needs derived values.
- [ ] **Page** in `www/src/components/pages/<Tool>Page.vue` (Part III, *The
      page*), and a help collapsible in `www/src/components/help/`.
- [ ] **Navigation in `www/src/App.vue`.** Import and register the page and its
      `lucide` icon, add a `v-else-if` branch for it, and add an entry to
      `tabs` (`id`, `sidebar_label`, `label`, `icon`); the `id` is also the URL
      hash. Start its worker or default pool in `mounted()` if it must be ready
      when the tab opens.
- [ ] **Saving** goes through `www/src/platform/files.ts` (`saveTextFile`,
      `saveBinaryFile`, `openBinaryFileSink`), never a raw download link, so the
      desktop program gets native save dialogs.
- [ ] **Run metrics:** the driver returns `elapsedMs` and `wasmMemoryBytes`, as
      the other pages show.

## wasm64, only if 4 GB is not enough

- [ ] `scripts/build_wasm64.py --crate rust/<crate> --module <tool>` must build
      it: one `cdylib`, and wasm-bindgen 0.2.120 or newer. SIMD or stack flags go
      in the script, which ignores the crate's `.cargo/config.toml`.
- [ ] `www/src/wasm/packages/<tool>.ts`, modelled on `ska.ts`. The driver loads
      it with `loadWasmPackage()`, and the worker is wrapped with
      `observeWasmRuntime()` where it is created.
- [ ] The same `build_wasm64.py` line in `cf_deploy.yaml` and
      `electron_release.yaml`, and the module in `electron_release.yaml`'s check
      that the packages reached the renderer.

## Tests and CI

- [ ] **A wasm check.** An entry in `TESTS` in `scripts/test_wasm.py` (crate,
      check, the same `wasm_pack_args` as `vue.config.js`, `wasm64`, `uses`); a
      `tests/wasm/<tool>.mjs` built on `browser.mjs`; its inputs in
      `tests/wasm/data/`, with their origin in `data/README.txt`; and the golden
      file its first run writes, reviewed and committed.
- [ ] **CI.** A row in the `wasm` matrix of `.github/workflows/ci.yaml`, whose
      `module` must match the `TESTS` key. An in-tree crate also goes in the
      `fmt`, `clippy` and `test` jobs and their `rust-cache` workspaces; its
      coverage then counts towards `.octocov.yml`'s threshold.

## Documentation

- [ ] A methods page, `docs/src/<tool>.md`, and its line in
      `docs/src/SUMMARY.md`.
- [ ] A "What can I do in the <Tool> tab?" entry in `FaqPage.vue`, and a line in
      the README's feature list.

**Done when:** everything in §10 holds, and `npm run lint`,
`npm run test:wasm -- --only <tool>` and `mdbook build` (in `docs/`) are green.

---

## Further reading — not required

- [`docs/src/wasm_guide.md`](docs/src/wasm_guide.md) — the same mechanics written
  for a human reader, with a narrative worked example, notes on wasm64 toolchain
  progress and a section on WebGPU via `wgpu`.
- [AGENTS.md](AGENTS.md) — this repository specifically: its crates, store,
  design tokens, wasm64, the desktop program, CI, tests and gotchas.
