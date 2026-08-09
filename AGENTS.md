# AGENTS.md — sparrowhawk

Orientation for an AI agent working in this repository. Everything below was read
from the tree on `main`; where something could not be verified it says so.

For *how to add a new tool*, read [SKILLS.md](SKILLS.md) — it is self-contained
and needs neither this file nor the guide. [`docs/src/wasm_guide.md`](docs/src/wasm_guide.md)
is the same WebAssembly material written for a human reader.

---

## What this is

A browser-based bioinformatics platform. Rust crates compile to WebAssembly and
run **entirely client-side** — no server, nothing uploaded. Each tool is a tab in
a single Vue 3 SPA.

The repository was renamed: this is `bacpop/sparrowhawk` (formerly
`sparrowhawk-web`), and the assembler crate is now `bacpop/sparrowhawk-asm`.
Older checkouts and documents may still use the old names.

Eight tools, one page each under `www/src/components/pages/`:

| Tool | Page | Crate |
|---|---|---|
| Assembly | `AssemblyPage.vue` | `sparrowhawk-asm` |
| Mapping / Alignment | `MappingAlignmentPage.vue` | `ska.rust` |
| Taxonomic ID | `TaxonomicIDPage.vue` | `sketchlib.rust` |
| Gene calling | `GeneCallingPage.vue` | `orphos-bridge` |
| Host depletion | `HostDepletionPage.vue` | `deacon-bridge` |
| AMR detection | `AMRDetectionPage.vue` | `sparrowhawk-amr` |
| Protein embeddings | `ProteinEmbeddingsPage.vue` | `esm-bridge` |
| Transmission | `TransmissionPage.vue` | (front-end analysis) |

Plus `FaqPage.vue`.

---

## Getting a working checkout

```sh
git clone --recurse-submodules https://github.com/bacpop/sparrowhawk.git
cd sparrowhawk/www && npm install && npm run serve
```

**Without `--recurse-submodules` the `rust/` submodule directories are empty and
every wasm build fails.** This is the most common way to arrive at a broken tree.
Recover with `git submodule update --init --recursive`.

Requires the Rust toolchain plus `wasm32-unknown-unknown`, and `wasm-pack`.

### Two facts that will otherwise cost you an hour

1. **There is no hot reload.** `www/vue.config.js` sets `hot: false`,
   `liveReload: false` and `watchFiles: { paths: [] }` deliberately — otherwise
   every file touch retriggers seven `wasm-pack` release builds. **Restart
   `npm run serve` to see any change**, front end included.
2. **First start is slow.** Seven crates compile in release mode before webpack
   runs.

---

## Layout

```
sparrowhawk/
├── docs/            mdBook: user docs + the WebAssembly developer guide
├── electron/        desktop shell
├── rust/            7 crates — 4 submodules, 3 in-tree
├── www/             the Vue 3 SPA
├── .github/workflows/   Cloudflare Pages deploys (app + docs)
└── netlify.toml     older deploy path, still present
```

### `rust/`

Four **git submodules**, per `.gitmodules`:

| Crate | Repo | Branch |
|---|---|---|
| `sparrowhawk-asm` | `bacpop/sparrowhawk-asm` | `master` |
| `ska.rust` | `bacpop/ska.rust` | `sparrowhawk-dev` |
| `sketchlib.rust` | `bacpop/sketchlib.rust` | `sparrowhawk-dev` |
| `sparrowhawk-amr` | `bacpop/sparrowhawk-amr` | `master` |

Three **in-tree bridge crates**, which wrap an upstream library for the browser:
`deacon-bridge`, `orphos-bridge`, `esm-bridge`. `esm-bridge` also carries
`build.rs`, `model/`, `scripts/` and `tests/`.

### `www/src/`

| Path | Contents |
|---|---|
| `assets/app.css` | the entire design system (see below) |
| `components/pages/` | one page per tool, plus `amr-detection/`, `protein-embeddings/`, `taxonomic-id/` subfolders |
| `components/ui/` | shadcn-vue primitives over `reka-ui` |
| `components/help/` | per-tool help collapsibles |
| `components/` | result displays, `MSAViewer/`, `SequenceViewer/`, `MinimisedSequenceViewer/`, `gene-calling/`, transmission views |
| `lib/utils.ts` | `cn()` — `twMerge(clsx(...))` |
| `platform/` | `files.ts`, `gpu.ts`, `electron.d.ts` — web-vs-Electron abstraction and WebGPU detection |
| `store/` | Vuex: `index.ts`, `state.ts`, `actions.ts`, `mutations.ts`, `getters.ts` |
| `workers/` | driver + `*.worker.ts` pair per tool |
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
effectively vestigial — Tailwind 4 does not need its `content` array.

`www/src/components/ui/` holds twelve primitives: `button`, `collapsible`,
`dialog`, `input`, `separator`, `sheet`, `sidebar`, `skeleton`, `slider`,
`table`, `tabs`, `tooltip`.

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

The app shell in `App.vue` is the shadcn `Sidebar` family
(`SidebarProvider` → `Sidebar` → `SidebarMenu`), with the main content in a
white rounded panel. **Navigation is a `tabName` string plus a `v-if`/`v-else-if`
chain — there is no router**, despite `@vue/cli-plugin-router` being a
devDependency.

---

## Rust → wasm

All seven crates expose a **stateful handle**: a `#[wasm_bindgen] pub struct`
constructed once, mutated by work methods, drained by getters. Beyond that there
are **two styles**, and the split is *not* submodule versus in-tree — it is
simply older code and newer code:

| Style | Crates | Input | Output | Errors |
|---|---|---|---|---|
| Older | `sparrowhawk-asm`, `ska.rust`, `sketchlib.rust`, `orphos-bridge` | `web_sys::File` | JSON `String` via the `json` crate | `panic!` through `console_error_panic_hook` |
| Newer | `sparrowhawk-amr`, `deacon-bridge`, `esm-bridge` | `&[u8]` | typed values, `JsValue` | `Result<_, JsValue>` |

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

**No crate uses `serde-wasm-bindgen`.** Results cross the boundary as JSON
strings that JavaScript re-parses. That works, but it is untyped on both sides;
if you add a tool, deriving `Serialize` with
`#[serde(rename_all = "camelCase")]` and returning through `serde-wasm-bindgen`
gives native JS objects that mirror your TypeScript interfaces directly.

### Conditional compilation

Wasm-only code is gated on **`target_family = "wasm"`** — not `target_arch =
"wasm32"`, and no longer by a cargo feature, which is why the `--features wasm`
argument is commented out in `vue.config.js`:

```rust
#[cfg(target_family = "wasm")]
#[wasm_bindgen]
pub fn init_panic_hook() { console_error_panic_hook::set_once(); }
```

`Cargo.toml` uses the same predicate to swap native-only dependencies out, e.g.
`[target.'cfg(not(target_family = "wasm"))'.dependencies]` in `sketchlib.rust`,
which is how rayon, `indicatif` and `clap` stay out of the wasm build.
`deacon-bridge` and `orphos-bridge` gate nothing at all — they are wasm-only by
construction.

### Progress reporting

Two crates report progress **from Rust, mid-run**, by binding `postMessage`
directly and calling it from inside the worker —
`sparrowhawk-asm/src/lib.rs` and `esm-bridge/src/wasm.rs`:

```rust
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_name = postMessage)]
    fn post_message(data: &JsValue);
}

#[cfg(target_family = "wasm")]
pub fn post_state(state: &str) {
    let obj = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&obj, &JsValue::from_str("assemblyState"),
                                 &JsValue::from_str(state));
    post_message(&obj.into());
}
```

called at named lifecycle points — `post_state("preprocess:start")`,
`"preprocess:saving"`, `"preprocess:end"`, `"assembly:saving"`. The other five
crates are silent until they finish, and their pages fall back to a spinner
derived from store getters. **Copy the `post_state` pattern for anything that
runs longer than a second or two.**

---

## Workers

One driver class plus one `*.worker.ts` per tool, in `www/src/workers/`:
`Assembler`, `Mapper`, `Sketcher`, `Caller` (gene calling), `Depleter` (host
depletion), `AmrDetector`, `Embedder`. Plus `amrIndex.ts`, `esmWeights.ts` and
`workers.d.ts`.

The worker file is a thin dispatcher keyed on boolean flags of the message
object; the driver holds the lazily-imported wasm module and the stateful handle.

**Pools exist for the parallelisable tools.** `store/actions.ts` has
`initSketchlibWorkers`, `initCallerWorkers`, `initAmrDetectorWorkers` and
`initEmbedderWorkers`, each of which terminates the previous pool and spawns
`numWorkers` fresh ones:

```ts
async initSketchlibWorkers(context, numWorkers: number) {
    for (const worker of state.workerState.workers_sketchlib) worker.terminate();
    const pool: Worker[] = [];
    for (let i = 0; i < numWorkers; i++) pool.push(new WorkerSketcher());
    commit("SET_WORKERS_SKETCHLIB", pool);
}
```

`numWorkers` comes from a **slider in the page** (`const numWorkers = ref(4)`),
watched so that moving it respawns the pool. Assembly, Mapping and Host depletion
are single workers created once in `App.vue`.

These are ordinary workers, each with its own wasm instance — **not** wasm
threads. Shared memory would need `SharedArrayBuffer` and therefore COOP/COEP
headers, which break static hosting.

---

## Store

Vuex 4, one flat root store, no modules. `state.ts`:

```ts
export interface RootState {
    workerState: WorkerState,
    readsFileNames: string | null,
    readsPreprocessing: ReadsPreprocessing,
    allResults: AllResults,              // assembly — unsuffixed, for historical reasons
    refSet: string | null,
    allResults_ska: AllResultsSka,
    errors: string,
    min_count: number,
    allResults_sketchlib: AllResultsSketchlib,
    allResults_orphos: AllResultsOrphos,
    allResults_deacon: AllResultsDeacon,
    allResults_amr: AllResultsAmr,
    allResults_esm: AllResultsEsm,
    esmRetry: EsmRetry,
    gpuAdapters: GpuAdapterInfo[],
    transmissionStandalone: TransmissionStandaloneResults,
    processingState: ProcessingState,
}
```

Per-tool results are `allResults_<tool>`; assembly is the unsuffixed
`allResults`. Types live in `www/src/types.ts`. Mutations mix
`SCREAMING_SNAKE` (worker handles) with `camelCase` (results). Actions
`postMessage` and assign `worker.onmessage` per dispatch rather than registering
a persistent listener.

---

## Build, deploy, quality

```sh
npm run serve            # dev server on :8080, no hot reload
npm run build            # dist/
npm run build:electron   # dist-electron/, sets SPARROWHAWK_TARGET=electron
npm run lint
```

`serve` and `build` run node with `--max-old-space-size=8192`; the build needs it.

`vue.config.js` registers **seven `WasmPackPlugin` instances**, one per crate,
all `forceMode: "production"`, writing to `src/pkg`, `src/pkg_ska`,
`src/pkg_sketchlib`, `src/pkg_orphos-bridge`, `src/pkg_deacon`, `src/pkg_amr`,
`src/pkg_esm`. `deacon-bridge` alone passes `--no-default-features`. Workers go
through `worker-loader` (matching both `.js` and `.ts`, with `ts-loader` in
`transpileOnly` mode), and `experiments.asyncWebAssembly` is on.

**CI** (`.github/workflows/`): `cf_deploy.yaml` builds and deploys the app to
Cloudflare Pages on push to `main` — `submodules: true`, Node 24, Rust stable
with `wasm32-unknown-unknown`, `Swatinem/rust-cache`. `docs_deploy.yaml` builds
the mdBook and deploys it as `sparrowhawk-docs`. `netlify.toml` is an older
deploy path that is still present.

**There is no test suite.** No `test` script, no jest or vitest, no `*.spec.*`.
`npm run lint` is the only automated gate, and CI does not run it. Any correctness
claim you make must be backed by something you ran yourself.

---

## Documentation

`docs/` is an mdBook (`docs/src/SUMMARY.md`):

- **Methods** — user-facing pages for assembly, mapping, alignment, taxonomic ID,
  gene calling, host depletion.
- **Developer guide** — `wasm_guide.md`, 1100 lines covering targets, toolchain,
  designing the Rust interface, `#[wasm_bindgen]`, passing data across the
  boundary, conditional compilation, `vue.config.js` wiring, the worker pattern,
  wasm limitations, and a complete minimum working example. **Read it before
  writing any Rust that has to reach the browser.**

---

## Constraints worth knowing

- **4 GB memory ceiling** from 32-bit wasm addressing. Large inputs need chunking
  or streaming; `deacon-bridge`'s `push_chunk` is the pattern.
- **No filesystem, no threads by default, limited SIMD** — see the guide's
  "Limitations" section.
- Crates that assume `std::fs`, spawn threads, or link C libraries generally will
  not compile to wasm without a bridge crate in between.
