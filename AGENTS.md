# AGENTS.md — sparrowhawk

Orientation for an AI agent working in this repository. Everything below was read
from the tree on `main`; where something could not be verified it says so. File
paths are the source of truth: where this file and the code disagree, trust the
code, and fix this file.

For *how to add a new tool*, read [SKILLS.md](SKILLS.md): Parts I–III are a
self-contained playbook for any Rust → WebAssembly → Vue tool, and Part IV is the
checklist for adding one to this repository. [`docs/src/wasm_guide.md`](docs/src/wasm_guide.md)
is the same WebAssembly material written for a human reader. `rust/mandrake` has
its own `AGENTS.md`.

---

## What this is

A bioinformatics platform whose Rust crates compile to WebAssembly and run
**entirely client-side** — no server, nothing uploaded. Each tool is a tab in a
single Vue 3 SPA, served as a website (https://sparrowhawk.bacpop.org) and also
packaged as a desktop program with Electron (`electron/`).

The repository was renamed: this is `bacpop/sparrowhawk` (formerly
`sparrowhawk-web`), and the assembler crate is now `bacpop/sparrowhawk-asm`.
Older checkouts and documents may still use the old names.

The tabs, in the order of `tabs` in `www/src/App.vue`, with their pages under
`www/src/components/pages/`:

| Tab | Page | Crate |
|---|---|---|
| Assembly | `AssemblyPage.vue` | `sparrowhawk-asm` |
| Mapping, Alignment | `MappingAlignmentPage.vue` (one page, two tabs) | `ska.rust` |
| Transmission | `TransmissionPage.vue` | `ska.rust` (`ska_cluster`, run in the ska worker) |
| Clustering | `MandrakePage.vue` | `mandrake` |
| Taxonomic ID | `TaxonomicIDPage.vue` | `sketchlib.rust` |
| Gene calling | `GeneCallingPage.vue` | `orphos-bridge` (its driver also loads `sparrowhawk-amr`) |
| Host depletion | `HostDepletionPage.vue` | `deacon-bridge` |
| AMR detection | `AMRDetectionPage.vue` | `sparrowhawk-amr` |
| Embed proteins | `ProteinEmbeddingsPage.vue` | `esm-bridge` |

Plus `FaqPage.vue`, reached as the `faq` tab.

---

## Getting a working checkout

```sh
git clone --recurse-submodules https://github.com/bacpop/sparrowhawk.git
cd sparrowhawk/www && npm install && npm run serve
```

**Without `--recurse-submodules` the `rust/` submodule directories are empty and
every wasm build fails.** This is the most common way to arrive at a broken tree.
Recover with `git submodule update --init --recursive`; `--recursive` matters,
because `rust/mandrake` has a submodule of its own.

Requires:

- the Rust toolchain plus `wasm32-unknown-unknown`;
- `wasm-pack` **0.15.0**, the version every workflow pins. Older ones — such as
  the 0.9.1 that `jetli/wasm-pack-action` falls back to when left unpinned — fail
  on the assembler and ska with "Only 1 table definition allowed in MVP";
- for the optional wasm64 packages only: nightly Rust with `rust-src`, and
  Python 3.9+ (see *wasm64*).

### Two facts that will otherwise cost you an hour

1. **There is no hot reload.** `www/vue.config.js` sets `hot: false`,
   `liveReload: false` and `watchFiles: { paths: [] }` deliberately — otherwise
   every file touch retriggers every crate's `wasm-pack` release build. **Restart
   `npm run serve` to see any change**, front end included.
2. **First start is slow.** Every crate compiles in release mode before webpack
   runs.

---

## Layout

```
sparrowhawk/
├── .github/workflows/  CI, the website and docs deploys, the desktop releases
├── docs/               mdBook: user docs + the WebAssembly developer guide
├── electron/           the desktop program: Electron shell and packaging
├── rust/               the crates: git submodules and in-tree bridges
├── scripts/            build_wasm64.py (wasm64 packages), test_wasm.py (wasm checks)
├── tests/wasm/         Node checks of each crate's wasm package, with goldens
├── www/                the Vue 3 SPA
├── CITATION.cff        citation metadata (GitHub's "Cite this repository")
└── .octocov.yml        coverage thresholds for CI
```

Builds also leave ignored local folders: `.wasm64-build/` and `.wasm64-tools/`
(wasm64 builds and their `wasm-bindgen-cli`), `.wasm-test/` (the checks'
packages) and `www/public/pkg_wasm64/`.

### `rust/`

**Git submodules**, per `.gitmodules`. Each is pinned to the commit recorded in
this repository; the branch is only what `git submodule update --remote` follows.

| Crate | Repo | Branch |
|---|---|---|
| `sparrowhawk-asm` | `bacpop/sparrowhawk-asm` | `master` |
| `ska.rust` | `bacpop/ska.rust` | `sparrowhawk-dev` |
| `sketchlib.rust` | `bacpop/sketchlib.rust` | `sparrowhawk-dev` |
| `sparrowhawk-amr` | `bacpop/sparrowhawk-amr` | `master` |
| `mandrake` | `bacpop/mandrake.rust` | `sparrowhawk-dev`, with a nested `sketchlib.rust` submodule |

**In-tree bridge crates**, which wrap an upstream library for the browser:
`deacon-bridge`, `orphos-bridge` and `esm-bridge`. The first two depend on
upstream by git branch, pinned in practice by their committed `Cargo.lock`.
`esm-bridge` also carries `build.rs`, `model/`, `scripts/` and `tests/`.

### `www/src/`

| Path | Contents |
|---|---|
| `App.vue` | the shell: sidebar, the `tabs` list, the page chain, worker start-up |
| `assets/app.css` | the entire design system (see below) |
| `components/pages/` | one page per tool, plus `amr-detection/`, `protein-embeddings/`, `taxonomic-id/` subfolders |
| `components/ui/` | shadcn-vue primitives over `reka-ui` |
| `components/help/` | one help collapsible per tool |
| `components/` | result displays and shared widgets: `MSAViewer/`, `SequenceViewer/`, `gene-calling/`, plots, download buttons, `WasmRuntimeWarnings.vue`, … |
| `lib/utils.ts` | `cn()` — `twMerge(clsx(...))` |
| `platform/` | `files.ts` (saving and asset loading, web or Electron), `gpu.ts` (WebGPU detection), `electron.d.ts` |
| `wasm/` | `loader.ts` and `packages/`: load a wasm64 package, falling back to wasm32; `runtime.ts` reports which one ran |
| `store/` | Vuex: `index.ts`, `state.ts`, `actions.ts`, `mutations.ts`, `getters.ts` |
| `workers/` | the workers and their drivers (see *Workers*) |
| `types.ts`, `utils.ts` | shared types; `emptyState()` (the store's initial state) and `clampWorkerCount()` |
| `pkg*/` | wasm-pack output, generated, not committed |

---

## The design system

All of it lives in `www/src/assets/app.css`:

```css
@import "tailwindcss";
@import "tw-animate-css";
@import "@fontsource/dm-sans/400.css";   /* …500, 600, 700 */

@custom-variant dark (&:is(.dark *));

:root {
    --background: oklch(1 0 0);
    /* …the full shadcn token set, plus --sidebar-* … */
    --radius: 0.625rem;
}
.dark { /* … */ }
```

**Tailwind 4**, configured through `postcss.config.js`
(`@tailwindcss/postcss` + `autoprefixer`). `tailwind.config.js` exists but is
effectively vestigial — Tailwind 4 does not need its `content` array. DM Sans is
the interface font; `App.vue` also loads IBM Plex.

`www/src/components/ui/` holds the primitives `button`, `collapsible`, `dialog`,
`input`, `separator`, `sheet`, `sidebar`, `skeleton`, `slider`, `table`, `tabs`
and `tooltip`.

Note there are **two UI libraries** in `package.json`: `reka-ui` (which the
`ui/` primitives wrap) and `primevue` + `tailwindcss-primeui`. `react` and
`react-dom` are present only because `taxonium-component` needs them.

### The page layout every tool follows

From `TaxonomicIDPage.vue` — copy this shape rather than inventing one:

```html
<div class="flex flex-col gap-6 md:flex-row md:gap-0">
  <div class="w-full md:w-[350px] md:shrink-0">
    <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
      <ScanFace class="w-6 h-6" /> Taxonomic ID
    </h1>
    <TaxonomicIDHelpCollapsible />
    <TooltipProvider>
      <div class="flex flex-col gap-4">
        <!-- one block per parameter -->
      </div>
    </TooltipProvider>
  </div>
  <!-- results column -->
</div>
```

- Parameters live in a **350 px left rail**; the drop zone and results go right.
- Every control is a label row carrying an `Info` icon at exactly
  `w-3.5 h-3.5 text-gray-400 cursor-help`, wrapped in `Tooltip`/`TooltipTrigger`
  (`as-child`) / `TooltipContent` with `<p class="max-w-xs">`.
- Sliders pair `VueSlider class="flex-grow"` with a `w-[40px]` readout box.
- Each help collapsible has `Overview` and `Parameters` tabs.

The app shell in `App.vue` is the shadcn `Sidebar` family
(`SidebarProvider` → `Sidebar` → `SidebarMenu`), with the main content in a
white rounded panel. **Navigation is a `tabName` string plus a `v-if`/`v-else-if`
chain — there is no router**, despite `@vue/cli-plugin-router` being a
devDependency. The tab is mirrored in the URL hash (`#Assembly`), so links and
back/forward work; a tab's `id` is its hash.

---

## Rust → wasm

Every crate exposes a **stateful handle**: a `#[wasm_bindgen] pub struct`
constructed once, mutated by work methods, drained by getters. Beyond that, the
styles differ, and the split is simply older code versus newer code:

- **Older**: `sparrowhawk-asm`, `sketchlib.rust`, `orphos-bridge`, and most of
  `ska.rust`. They take a `web_sys::File`, return JSON strings built with the
  `json` crate, and `panic!` on bad input (shown through
  `console_error_panic_hook`).
- **Newer**: `sparrowhawk-amr`, `deacon-bridge`, `esm-bridge`, `mandrake`, and
  ska's newer functions. They return `Result<_, JsValue>`, and take `&[u8]`
  (`mandrake` also takes `web_sys::File`). Outputs vary: typed arrays
  (`esm-bridge`, `mandrake`), small JS objects (`deacon-bridge`'s `stats()`), or
  JSON strings (`sparrowhawk-amr`).

**Prefer the newer style.** From `rust/sparrowhawk-amr/src/lib.rs`:

```rust
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
```

A real constructor, explicit errors, and `&[u8]` so large inputs can be streamed
in chunks — `deacon-bridge`'s `push_chunk(&mut self, chunk: &[u8]) ->
Result<Vec<u8>, JsValue>` / `finish()` pair is the model for anything that will
not fit in memory at once.

**No crate uses `serde-wasm-bindgen`.** For structured results, deriving
`Serialize` with `#[serde(rename_all = "camelCase")]` and returning through
`serde-wasm-bindgen` gives native JS objects that mirror your TypeScript
interfaces directly, instead of JSON strings that are untyped on both sides.

### Conditional compilation

Wasm-only code is gated on **`target_family = "wasm"`** — not `target_arch =
"wasm32"`, which silently drops the code from a wasm64 build — and mostly no
longer by a cargo feature (the commented-out `--features wasm` lines in
`vue.config.js` are leftovers):

```rust
#[cfg(target_family = "wasm")]
#[wasm_bindgen]
pub fn init_panic_hook() { console_error_panic_hook::set_once(); }
```

The exception is `mandrake`: its browser code needs the `wasm-sketchlib` feature
(`vue.config.js` passes `--no-default-features --features wasm-sketchlib`).

`Cargo.toml` uses the same predicate to keep native-only dependencies out, as
`[target.'cfg(not(target_family = "wasm"))'.dependencies]`; how much each crate
moves there varies (`sketchlib.rust` and `ska.rust` keep `rayon` and `clap` in
the wasm build too). `deacon-bridge` and `orphos-bridge` gate nothing at all —
they are wasm-only by construction.

**Each crate has its own `.cargo/config.toml`** adding `+simd128` for
`wasm32-unknown-unknown`; `orphos-bridge` and `esm-bridge` also raise the stack
to 8 MB. Cargo never inherits these from a dependency, and they apply only when
cargo runs from the crate's folder, as the wasm-pack plugin and `test_wasm.py`
do. The wasm64 build does not read them (see *wasm64*).

### Progress reporting

Three crates report progress **from Rust, mid-run**, by binding `postMessage`
and calling it from inside the worker:

- `sparrowhawk-asm` (`src/lib.rs`): `post_state()` posts `{assemblyState}` at
  named stages — `preprocess:start`, `preprocess:end`, `assembly:create_graph`,
  `assembly:saving`, `assembly:end`, and more;
- `esm-bridge` (`src/wasm.rs`): `post_progress()` posts
  `{embedProgress, done, total}`, and `report_event()` posts `{gpuEvent, message}`;
- `ska.rust` (`src/lib.rs`): `post_alignment_progress()` posts
  `{alignmentProgress: {…}}`.

The assembler's version, with a no-op twin so native code can call it too:

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

#[cfg(not(target_family = "wasm"))]
pub(crate) fn post_state(_state: &str) {}
```

`mandrake` works differently: its worker calls the handle's `advance()`
repeatedly, and each call returns the progress so far. The other crates are
silent until they finish. **Copy one of these patterns for anything that runs
longer than a second or two**, and call a count field `completed`, not `done`
(SKILLS.md §8): `esm-bridge`'s `done` predates that rule.

---

## Workers

Drivers and workers live in `www/src/workers/`. Most tools have a driver class
plus a `*.worker.ts`: `Assembler`, `Mapper` (mapping, alignment and
transmission), `Sketcher`, `Caller` (gene calling), `Depleter` (host
depletion), `AmrDetector` and `Embedder`. The exceptions:

- `AlignmentExtractor.worker.ts` has no driver; `processQueryAlign` spawns a
  pool of them.
- `Mandrake.ts` is a main-thread runner that starts a fresh `Mandrake.worker.ts`
  for each operation, outside the store.

Also there: `amrIndex.ts` and `esmWeights.ts` (asset loading), and
`workers.d.ts`, with one `declare module` per worker (a new worker needs one).

The worker file is a thin dispatcher keyed on boolean flags of the message
object (`Mandrake.worker.ts` switches on a `type` field instead); the driver
holds the lazily-imported wasm module and the stateful handle.

**Pools exist for the parallelisable tools.** `store/actions.ts` has
`initSketchlibWorkers`, `initCallerWorkers`, `initAmrWorkers` and
`initEmbedderWorkers`. Each terminates the previous pool, clamps the size with
`clampWorkerCount()` (in `www/src/utils.ts`: at most
`navigator.hardwareConcurrency`), spawns the new workers and commits them;
`initSketchlibWorkers` is the model to copy.

`App.vue` starts the default pools at mount (gene calling and taxonomic ID with
four workers, AMR with one); Taxonomic ID and Gene calling have a threads slider
that respawns their pool, and Embed proteins starts its own. Mapping and Host
depletion are single workers created in `App.vue`. The assembler's worker comes
from the `spawnAssemblerWorker` action, which also re-spawns it after a crash.

Workers that may load a wasm64 package report which build ran: they are wrapped
with `observeWasmRuntime()` (`www/src/wasm/runtime.ts`), which commits
`recordWasmRuntimeStatus`, and `WasmRuntimeWarnings.vue` shows a banner when
wasm64 fell back to wasm32.

These are ordinary workers, each with its own wasm instance — **not** wasm
threads. Shared memory would need `SharedArrayBuffer` and therefore COOP/COEP
headers, which break static hosting.

---

## Store

Vuex 4, one flat root store, no modules. The fields are declared in `RootState`
(`store/state.ts`), and their initial values are set in `emptyState()`, which is
in **`www/src/utils.ts`**, not in `store/`. Types live in `www/src/types.ts`.

Per-tool results are `allResults_<tool>`; assembly's is the unsuffixed
`allResults`, for historical reasons. Other fields worth knowing: `workerState`
(every worker and pool), `wasmRuntimeStatus` (which build each module loaded),
`errors` and `errorMessage` (the assembler's), `gpuAdapters`, and
`processingState`.

Mutations mix `SCREAMING_SNAKE` (worker handles, `SET_WORKER*`) with
`camelCase` (results). Actions `postMessage` and mostly assign
`worker.onmessage` per dispatch; a few add and remove `message` listeners
instead.

---

## Build, deploy, quality

```sh
npm run serve            # dev server on :8080, no hot reload
npm run build            # dist/
npm run build:electron   # dist-electron/, sets SPARROWHAWK_TARGET=electron
npm run lint             # fixes files in place; CI runs it with --no-fix
npm run test:wasm        # scripts/test_wasm.py: build and check every wasm package
```

The three builds run node with `--max-old-space-size=8192`; they need it.

`vue.config.js` registers one `WasmPackPlugin` per crate, all
`forceMode: "production"`, each writing to `www/src/pkg_<name>` (the
assembler's is plain `src/pkg`). Two pass `extraArgs`: `deacon-bridge`
(`--no-default-features`) and `mandrake`
(`--no-default-features --features wasm-sketchlib`). `publicPath` is `'/'`: the
app loads its assets by root-relative URLs, so it cannot be served from a
subdirectory. Workers go through `worker-loader` (matching both `.js` and `.ts`,
with `ts-loader` in `transpileOnly` mode), named
`js/[name].[contenthash].worker.js` so that a worker's wasm32 and wasm64
variants do not collide; `experiments.asyncWebAssembly` is on.

### CI and deploys

| Workflow | Runs on | Does |
|---|---|---|
| `ci.yaml` | pull requests, pushes to `main` | `fmt` and `clippy -D warnings` on the in-tree bridges; `build` (deacon-bridge for wasm32); `test` (orphos- and esm-bridge, with coverage reported by octocov against `.octocov.yml`); `lint`; `docs` (`mdbook build`); and `wasm`, one job per module running `scripts/test_wasm.py --only <module>` |
| `cf_deploy.yaml` | every push, on any branch | builds the wasm64 packages and the site, and deploys to Cloudflare Pages (`sparrowhawk`) |
| `docs_deploy.yaml` | pushes to `main` | builds the mdBook and deploys it as `sparrowhawk-docs` |
| `electron_release.yaml` | a published release, or by hand | builds the renderer once (with wasm64), packages an AppImage, a Windows installer and an Apple-silicon `.dmg`, and attaches them to the release |

Every workflow that runs wasm-pack pins it to 0.15.0 and checks out submodules
recursively.

### Tests

- **`tests/wasm/`**: one Node 24 script per crate (`<tool>.mjs`). It loads the
  crate's web-target wasm package through `browser.mjs` (which stands in for
  `File`, `FileReaderSync` and `postMessage`), replays the worker's calls, and
  compares deterministic outputs with `golden/<tool>-<target>.json`. Inputs are
  committed in `tests/wasm/data/`. Run it with `npm run test:wasm`
  (`-- --only <module>`, `--no-build`, `--update-golden`); `tests/wasm/README.md`
  explains it. Host depletion is not covered.
- **Rust unit tests** in the in-tree bridges; `esm-bridge`'s `tests/fidelity.rs`
  needs a GPU and is ignored by default.
- **No front-end tests**: no jest, vitest or `*.spec.*` in `www/`. Back any claim
  about the interface with something you ran yourself.

---

## wasm64

`sparrowhawk-asm` and `ska.rust` are also built for wasm64 (Memory64), which
lifts the 4 GB ceiling to 16 GiB:

- `scripts/build_wasm64.py --crate <path> --module <name>` builds with nightly
  Rust (`-Z build-std`) and `wasm-bindgen-cli`, not wasm-pack, into
  `www/public/pkg_wasm64/<module>/`, a web-target package served as a static
  file. It needs exactly one `cdylib` and wasm-bindgen 0.2.120 or newer, and it
  ignores the crate's `.cargo/config.toml`: its flags are set in the script
  (simd128 for ska only, since it breaks `fixedbitset` in the assembler).
- `www/src/wasm/packages/<module>.ts` describes each package. The driver loads
  it with `loadWasmPackage()` (`www/src/wasm/loader.ts`), which tries wasm64 and
  **falls back to wasm32 silently** if wasm64 fails or is missing; only a banner
  says so.
- `cf_deploy.yaml` and `electron_release.yaml` build both packages before the
  site; `electron_release.yaml` fails if they did not reach the renderer.
- A local `npm run serve` has only wasm32 unless you ran the script first.

---

## The desktop program

`electron/` packages the same renderer as a desktop program:

- `npm run build:electron` (in `www/`) builds the renderer into
  `www/dist-electron/`. `electron/package.json` holds electron-builder's
  configuration (AppImage x64, NSIS x64, dmg arm64, all unsigned) and the
  `dev` and `build:*` scripts.
- `main.js` serves the renderer from a privileged `app://bundle/` scheme, so
  root-relative URLs work as on the website; on Linux it enables Vulkan for
  WebGPU.
- `preload.js` exposes `window.electronAPI`: save dialogs, streamed alignment
  saves, and reading the taxonomic index. The front end reaches it only through
  `www/src/platform/files.ts`; save and load files with its functions
  (`saveTextFile`, `saveBinaryFile`, `openBinaryFileSink`, `loadAssetBytes`, …)
  so a tool works in both.
- The version is in `www/package.json` and `electron/package.json`; release
  builds take it from the tag instead.

---

## Documentation

`docs/` is an mdBook (`docs/src/SUMMARY.md`):

- **Methods** — user-facing pages for assembly, mapping, alignment, clustering,
  taxonomic ID, gene calling and host depletion. AMR detection, protein
  embeddings and transmission have none yet.
- **Developer guide** — `wasm_guide.md`, covering targets (with a note on
  wasm64), toolchain, designing the Rust interface, `#[wasm_bindgen]`, passing
  data across the boundary, conditional compilation, `vue.config.js` wiring, the
  worker pattern, wasm limitations, WebGPU, and a complete minimum working
  example. **Read it before writing any Rust that has to reach the browser.**

---

## Constraints worth knowing

- **4 GB memory ceiling** on wasm32: every crate except the two wasm64 ones, and
  those too when they fall back. Large inputs need chunking or streaming;
  `deacon-bridge`'s `push_chunk` is the pattern.
- **No filesystem, no threads by default, limited SIMD** — see the guide's
  "Limitations" section.
- Crates that assume `std::fs`, spawn threads, or link C libraries generally will
  not compile to wasm without a bridge crate in between.
- **Large data files** (the AMR index, the ESM weights, the taxonomic-ID index)
  are committed in `www/public/` and loaded from the site root through
  `platform/files.ts`.
- **Every `www/src/pkg*/` directory is generated by `wasm-pack` and ignored.**
  `.gitignore` uses the glob `www/src/pkg*/`; the narrower `www/src/pkg` it
  replaced covered only the assembler's output. `www/.eslintignore`, however,
  lists each directory by name, so a new crate needs its own line there.
