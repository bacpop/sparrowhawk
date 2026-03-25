# Guide to compile methods or tools to WebAssembly

This guide is aimed at bioinformaticians who are comfortable writing Rust and want to bring their tools to the web via WebAssembly (wasm). It covers the essential concepts, toolchain, workflow, and known limitations, using the same approach as this repository: Rust compiled with `wasm-pack`, and a plain-JavaScript front end that calls into the resulting module.

---

## Table of Contents

1. [What is WebAssembly?](#what-is-webassembly)
2. [Available targets](#available-targets)
3. [Toolchain](#toolchain)
4. [Workflow](#workflow)
5. [Designing the Rust interface](#designing-the-rust-interface)
6. [Building and bundling](#building-and-bundling)
7. [Web integration](#web-integration)
8. [Limitations](#limitations)
9. [A note on WebGPU](#a-note-on-webgpu)
10. [Minimum working example (MWE)](#minimum-working-example-mwe)

---

## What is WebAssembly?

WebAssembly (wasm) is a binary instruction format designed as a portable compilation target for high-level languages, though can be loosenly defined as **a somewhat-universal compilation target**. It can run inside the browser (and increasingly outside of it, via runtimes such as wasmtime or WASI), with near-native performance. This in particular allows to run the same code (or close to it) as one would run locally in the terminal for doing bioinformatics analyses, in a website under whichever platform. This also allows to use languages such as Rust to bear the burden of the most computing-demanding part of the algorithms, instead of using JavaScript, that might be slower for those.

Another main feature from wasm is that allows for running methods fully locally. With this compilation target, as the exact same code could be potentially run with it, no upload of input data (e.g. a genome) to some remote server is needed, as many other websites offer, the same goes for downloading the results. This makes it ideally for resource-constrained environments (in particular, bad Internet connection), as well as for working with sensitive data, as it does not need to leave your computer.

From the browser's perspective, a wasm module is essentially a well-defined binary blob that exports functions and, optionally, memory. JavaScript loads the module, instantiates it, and calls its exports as ordinary function calls. There is no subprocess, no file system access by default, and no shared state beyond what the module explicitly exposes.

---

## Available Targets

Rust expresses compilation targets using the standard `<arch>-<vendor>-<os>` triple. The most relevant wasm targets for us are:

| Target                       | Status in Rust | Notes                                                                                              |
| ---------------------------- | -------------- | -------------------------------------------------------------------------------------------------- |
| `wasm32-unknown-unknown`     | Tier 2         | No OS abstraction; the standard target for browser-facing wasm. Used by this project.             |
| `wasm64-unknown-unknown`     | Tier 3         | 64-bit address space (Memory64 proposal); support in `wasm-bindgen` and `wasm-pack` is in active development but not yet released (see below). |
<!-- | `wasm32-wasip1`              | Tier 2         | Targets the WASI (WebAssembly System Interface) preview 1; adds a POSIX-like layer.               | -->
<!-- | `wasm32-wasip2`              | Tier 2         | Targets WASI preview 2, based on the Component Model.                                              | -->

For browser-facing Rust code combined with `wasm-bindgen` and `wasm-pack`, **`wasm32-unknown-unknown` is the only practical choice**. It provides no operating system layer whatsoever — no file I/O, no threading primitives, no environment variables — which is both its strength (minimal, portable binary) and a source of constraints discussed in the [Limitations](#limitations) section.

### A note on wasm64

The `wasm64-unknown-unknown` target, corresponding to the WebAssembly Memory64 proposal, would allow modules to address more than 4 GiB of linear memory. This is potentially very relevant for bioinformatics tools operating on large reference genomes or graphs. However, as of the time of writing, toolchain support is still in progress and has not yet landed in a stable release:

- `wasm-bindgen`: PR [#5004](https://github.com/wasm-bindgen/wasm-bindgen/pull/5004) — "Add wasm64/memory64 support for `wasm64-unknown-unknown` target" — is open and pending review. It was filed in February 2026 and closes the two underlying tracking issues ([#4436](https://github.com/wasm-bindgen/wasm-bindgen/issues/4436) and [#4499](https://github.com/wasm-bindgen/wasm-bindgen/issues/4499)).
- `wasm-pack`: open issue [#1464](https://github.com/drager/wasm-pack/issues/1464), with a companion PR [#1553](https://github.com/drager/wasm-pack/pull/1553) linked from the `wasm-bindgen` PR above.

Until both PRs are merged and released, `wasm32-unknown-unknown` remains the only usable target for browser-facing Rust with `wasm-pack`. Memory-intensive computations must be structured to fit within the 4 GiB linear memory limit, and the interface between Rust and JavaScript should be designed accordingly.

---

## Toolchain

Install the following tools before starting:

```bash
# Install the wasm target for the Rust standard library
rustup target add wasm32-unknown-unknown

# Install wasm-pack — handles compilation, wasm-bindgen, and JS/TS glue generation
cargo install wasm-pack
```

The central tool is [`wasm-pack`](https://rustwasm.github.io/wasm-pack/). It orchestrates, so that you don't need to take care of this:

1. Calling `cargo build --target wasm32-unknown-unknown`.
2. Running `wasm-bindgen` to generate the JavaScript/TypeScript glue layer.
3. Optionally running `wasm-opt` for further binary optimisation.
4. Producing a `pkg/` directory containing the `.wasm` file, the generated JS module, and a TypeScript declaration file.

[`wasm-bindgen`](https://github.com/wasm-bindgen/wasm-bindgen) is the lower-level library that makes it possible to annotate Rust types and functions so that they can be called from JavaScript with idiomatic types (strings, typed arrays, JS objects) rather than raw integer pointers.



---

## Workflow

The overall final way of running Rust code with WebAssembly in your website will look like this:

```
┌─────────────────────────────────────────────────────────────────┐
│  Rust library crate (crate-type = ["cdylib"])                   │
│                                                                 │
│  1. Annotate public API with #[wasm_bindgen]                    │
│  2. Ensure that it compiles for your wasm target                │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  Web front end (JavaScript / TypeScript / Vue / …)              │
│                                                                 │
│  3. Use wasm-pack to automatically compile, bind (wasm-bindgen),│
│     optimise (wasm-opt), and put the binaries in a folder in    │
│     your web repository                                         │
│  4. Initialise the wasm module (async, one-time)                │
│  5. Call exported Rust functions or construct exported structs  │
│  6. Receive results as JS-native types                          │
└─────────────────────────────────────────────────────────────────┘
```

---

## Designing the Rust interface
### `Cargo.toml` setup

The crate must be declared as a `cdylib` (C-compatible dynamic library), which is the form wasm modules take:

```toml
[lib]
crate-type = ["cdylib"]

[dependencies]
### Whatever version you want/need, here 0.2 is an example
wasm-bindgen = "0.2"

# For returning errors to JS:
# wasm-bindgen = { version = "0.2", features = ["serde-serialize"] }
```

If you want, you can also declare it as `rlib` in the same crate-type variable in the toml file.


### Annotating with `#[wasm_bindgen]`

The `#[wasm_bindgen]` attribute is the bridge between Rust and JavaScript. It can annotate free functions or `impl` blocks on structs.

**Option A — Free functions** (simplest, good for stateless operations):

```rust
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn gc_content(sequence: &str) -> f64 {
    let len = sequence.len() as f64;
    if len == 0.0 {
        return 0.0;
    }
    let gc = sequence
        .bytes()
        .filter(|&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
        .count() as f64;
    gc / len
}
```

**Option B — Exported struct with methods** (better for stateful objects, e.g., a k-mer index or an assembler graph that must persist across multiple JS calls):

```rust
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct KmerCounter {
    counts: std::collections::HashMap<String, usize>,
    k: usize,
}

// In this example, all these methods will become accessible to JavaScript code
#[wasm_bindgen]
impl KmerCounter {
    // Here we re-declare this particular method with #[wasm_bindgen(constructor)],
    // as this will allow us to later, in JavaScript, just create an instance of this
    // object, exactly as if it was a normal instance in JavaScript.
    #[wasm_bindgen(constructor)]
    pub fn new(k: usize) -> Self {
        Self {
            counts: std::collections::HashMap::new(),
            k,
        }
    }

    pub fn add_sequence(&mut self, sequence: &str) {
        for i in 0..sequence.len().saturating_sub(self.k - 1) {
            let kmer = &sequence[i..i + self.k];
            *self.counts.entry(kmer.to_string()).or_insert(0) += 1;
        }
    }

    pub fn count(&self, kmer: &str) -> usize {
        *self.counts.get(kmer).unwrap_or(&0)
    }

    pub fn total_kmers(&self) -> usize {
        self.counts.values().sum()
    }
}
```


### Passing data across the boundary

The types that `wasm-bindgen` handles transparently include: `bool`, `i32`/`u32`/`i64`/`u64`/`f32`/`f64`, `String`/`&str`, `Vec<u8>`/`&[u8]` (as `Uint8Array` in JS), and `JsValue`. For returning richer structured data (e.g., multiple fields), consider:

- Returning a `JsValue` built with `serde_wasm_bindgen` (requires the `serde` feature).
- Using multiple getter methods on an exported struct.
- Encoding the result as JSON and returning it as a `String`.

Avoid passing very large strings or byte slices across the boundary on every call; instead, pass them once to a stateful object and keep results on the Rust side until they are needed.

### Conditional compilation

It is good practice to gate `wasm-bindgen` imports behind a `cfg` attribute so that the crate still compiles and can be tested natively:

```rust
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
```

This is also a way to have, in the same crate/repository, code that compiles to `x86_64`, and also to WebAssembly, as you can essentially create alternative "versions" of those methods, functions, or crates that need to be updated.

---

## Building and bundling, automatically

In sparrowhawk-web, we use the `@wasm-tool/wasm-pack-plugin` to integrate WASM compilation directly into the Vue build process. This approach automatically handles compilation, binding generation, and module integration.

### Vue Configuration with wasm-pack-plugin

The `vue.config.js` file configures the wasm-pack-plugin for each Rust crate:

```javascript
const WasmPackPlugin = require("@wasm-tool/wasm-pack-plugin");
const path = require("path");

module.exports = {
    configureWebpack: {
        experiments: {
            asyncWebAssembly: true,
        },
    },
    chainWebpack: (config) => {
        // Configure wasm-pack plugin for your Rust crate
        config
            .plugin("wasm-pack_your_crate")
            .use(WasmPackPlugin)
            .init(
                (Plugin) =>
                    new Plugin({
                        crateDirectory: path.resolve(__dirname, "../rust/your-crate"),
                        outDir: path.resolve(__dirname, "./src/pkg"),
                        forceMode: "production",
                    })
            )
            .end()
    },
}
```

### Key Configuration Options

- `crateDirectory`: Path to your Rust crate directory
- `outDir`: Output directory for generated files (relative to Vue project)
- `forceMode`: Set to `"production"` for optimized builds

### Build Process

When you run `npm run build` or `npm run serve`, the wasm-pack-plugin automatically:

1. Compiles your Rust crate to WASM using `wasm32-unknown-unknown` target
2. Runs `wasm-bindgen` to generate JavaScript/TypeScript glue code
3. Optimizes the WASM binary with `wasm-opt`
4. Places the output files in the specified `outDir`

The output in `src/pkg/` will contain:
- `<crate_name>_bg.wasm` — the compiled WASM binary
- `<crate_name>.js` — the JavaScript glue module
- `<crate_name>.d.ts` — TypeScript declarations
- `package.json` — npm package metadata

### Manual Build (Alternative)

If you need to build manually outside the Vue build process:

```bash
# Navigate to your Rust crate directory
cd ../rust/your-crate

# Build with wasm-pack
wasm-pack build --target web --release --no-default-features --features wasm

# Copy the pkg directory to your Vue project
cp -r pkg ../web/www/src/pkg
```

> **Note**: The manual approach requires you to manually copy files and may not integrate as smoothly with the Vue development server.

<!--## Calling from JavaScript/TypeScript

<a name="calling-from-javascripttypescript"></a>

In sparrowhawk-web, we use a worker pattern to handle WASM operations, which provides better performance and avoids blocking the main thread. Here's how to integrate WASM with Vue and TypeScript:

### Vue Component with Worker Pattern

```html
<script lang="ts">
import { defineComponent, ref } from 'vue';

// Worker that handles WASM operations
export default defineComponent({
  name: 'WasmExample',
  setup() {
    const result = ref<string>('');
    const isLoading = ref<boolean>(false);
    
    const calculateGCContent = async (sequence: string) => {
      isLoading.value = true;
      try {
        // Create worker and handle WASM operations
        const worker = new Worker(new URL('./workers/WasmWorker.ts', import.meta.url));
        
        worker.postMessage({ sequence });
        
        worker.onmessage = (event) => {
          if (event.data.result) {
            result.value = `GC Content: ${(event.data.result * 100).toFixed(2)}%`;
          }
          worker.terminate();
          isLoading.value = false;
        };
        
        worker.onerror = (error) => {
          console.error('WASM worker error:', error);
          result.value = 'Error calculating GC content';
          worker.terminate();
          isLoading.value = false;
        };
      } catch (error) {
        console.error('Failed to create worker:', error);
        result.value = 'Failed to initialize WASM';
        isLoading.value = false;
      }
    };
    
    return { result, isLoading, calculateGCContent };
  }
});
</script>

<template>
  <div>
    <h2>GC Content Calculator</h2>
    <textarea v-model="sequence" placeholder="Enter DNA sequence..."></textarea>
    <button @click="calculateGCContent(sequence)" :disabled="isLoading">
      {{ isLoading ? 'Calculating...' : 'Calculate GC Content' }}
    </button>
    <div>{{ result }}</div>
  </div>
</template>
```

### Worker Implementation

The worker handles the WASM module loading and function calls:

```typescript
// src/workers/WasmWorker.ts
interface WasmModule {
  gc_content: (sequence: string) => number;
  // Add other exported functions here
}

type WasmModuleAny = any;

self.onmessage = async (event) => {
  try {
    // Dynamically import the WASM module
    const wasm: WasmModuleAny = await import('@/pkg');
    
    // Call the WASM function
    const gc = wasm.gc_content(event.data.sequence);
    
    // Send result back to main thread
    self.postMessage({ result: gc });
  } catch (error) {
    console.error('WASM error:', error);
    self.postMessage({ error: error.message });
  }
};
```

### Direct Import Approach (Alternative)

For simpler cases, you can import WASM directly in components:

```typescript
// In a Vue component or utility function
let wasmModule: any = null;

async function loadWasm() {
  if (!wasmModule) {
    wasmModule = await import('@/pkg');
  }
  return wasmModule;
}

async function calculateGC(sequence: string): Promise<number> {
  const wasm = await loadWasm();
  return wasm.gc_content(sequence);
}
```

### TypeScript Interfaces

For better type safety, define interfaces for your WASM module:

```typescript
// types/wasm.d.ts
declare module '@/pkg' {
  export function gc_content(sequence: string): number;
  
  export class KmerCounter {
    constructor(k: number);
    add_sequence(sequence: string): void;
    count(kmer: string): number;
    total_kmers(): number;
    free(): void;
  }
  
  export function init(): Promise<void>;
}
```

> **Memory management**: When using exported structs, you must call `.free()` to release WASM-allocated memory. In the worker pattern, this is typically handled within the worker when it terminates.

> **Performance tip**: The worker pattern is recommended for computationally intensive operations to avoid blocking the UI thread. WASM modules are cached, so subsequent calls are faster.-->

## Web integration

<a name="web-integration"></a>

One simple way of dealing with your Rust code in your website is through **web workers**: these are, in essence, independent threads that you can focus on simply talking to your WebAssembly interface.

### Architecture overview

```
┌─────────────────────────────────────────────────────────────────┐
│  Web (Vue) components (main thread)                             │
│                                                                 │
│  1. User interacts with UI                                      │
│  2. Component creates Worker                                    │
│  3. Sends data to Worker via postMessage                        │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  Web Worker (background thread)                                 │
│                                                                 │
│  4. Worker dynamically imports WASM module                      │
│  5. Calls WASM functions                                        │
│  6. Processes results                                           │
│  7. Sends results back to main thread                           │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  WASM module (Rust compiled to WebAssembly)                     │
│                                                                 │
│  8. Executes bioinformatics algorithms                          │
│  9. Returns computed results to the worker                      │
└─────────────────────────────────────────────────────────────────┘
```

### Worker class implementation

Here's a more complete example showing how you can write your web workers. You will create usually two files, one with the actual class definition, and other with how to deal with the messages that you will exchange between it and the main thread that runs your website. Let's go with the first file (class declaration):

```typescript
// src/workers/BioWorker.ts
interface WasmModule {
  // Define your WASM module interface
  process_sequence: (sequence: string, params: any) => string;
  // Add other exported functions and classes
}

type WasmModuleAny = any;

export class BioWorker {
  private worker: Worker;
  private wasm: WasmModuleAny | null = null;
  private wasmPromise: Promise<WasmModuleAny>;
  
  constructor() {
    this.worker = new Worker(new URL('./bio-worker.impl.ts', import.meta.url));
    this.wasm = null;
    
    this.wasmPromise = new Promise((resolve) => {
      this.worker.onmessage = (event) => {
        if (event.data.type === 'wasm_loaded') {
          this.wasm = event.data.wasm;
          resolve(this.wasm);
        }
      };
    });
  }
  
  private async waitForWasm(): Promise<WasmModuleAny> {
    return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
  }
  
  public async processData(sequence: string, params: any): Promise<string> {
    await this.waitForWasm();
    
    return new Promise((resolve, reject) => {
      const messageId = Math.random().toString(36).substring(2);
      
      this.worker.postMessage({
        type: 'process',
        messageId,
        sequence,
        params
      });
      
      const handler = (event: MessageEvent) => {
        if (event.data.messageId === messageId) {
          this.worker.removeEventListener('message', handler);
          if (event.data.error) {
            reject(event.data.error);
          } else {
            resolve(event.data.result);
          }
        }
      };
      
      this.worker.addEventListener('message', handler);
    });
  }
  
  public terminate() {
    this.worker.terminate();
  }
}
```

And now, the second file:

```typescript
// src/workers/bio-worker.impl.ts
interface WasmModule {
  process_sequence: (sequence: string, params: any) => string;
  init?: () => Promise<void>;
}

type WasmModuleAny = any;

let wasm: WasmModuleAny | null = null;

self.onmessage = async (event) => {
  try {
    switch (event.data.type) {
      case 'init':
        // Load WASM module
        wasm = await import('@/pkg');
        if (wasm.init) {
          await wasm.init();
        }
        self.postMessage({ type: 'wasm_loaded', wasm });
        break;
        
      case 'process':
        if (!wasm) {
          throw new Error('WASM not initialized');
        }
        
        // Process the sequence with WASM
        const result = wasm.process_sequence(event.data.sequence, event.data.params);
        
        self.postMessage({
          messageId: event.data.messageId,
          result: result
        });
        break;
    }
  } catch (error) {
    console.error('Worker error:', error);
    self.postMessage({
      messageId: event.data?.messageId,
      error: error.message
    });
  }
};
```

### Using the Worker in Vue Components

And here you can see and example of how to call the previous worker example in a web (a Vue "component").

```vue
<script lang="ts">
import { defineComponent, ref, onMounted, onUnmounted } from 'vue';
import { BioWorker } from '@/workers/BioWorker';

export default defineComponent({
  name: 'BioProcessing',
  setup() {
    const result = ref('');
    const isProcessing = ref(false);
    const bioWorker = ref<BioWorker | null>(null);
    
    onMounted(() => {
      // Initialize worker when component mounts
      bioWorker.value = new BioWorker();
    });
    
    onUnmounted(() => {
      // Clean up worker when component unmounts
      bioWorker.value?.terminate();
    });
    
    const processSequence = async (sequence: string) => {
      if (!bioWorker.value) return;
      
      isProcessing.value = true;
      try {
        const processed = await bioWorker.value.processData(sequence, {
          // Add processing parameters
        });
        result.value = processed;
      } catch (error) {
        result.value = `Error: ${error.message}`;
      } finally {
        isProcessing.value = false;
      }
    };
    
    return { result, isProcessing, processSequence };
  }
});
</script>
```

### Benefits of the worker pattern

1. **Non-blocking UI**: Computationally intensive operations don't freeze the user interface
2. **Better memory management**: Workers can be terminated to clean up memory
3. **Reusability**: Worker logic can be shared across multiple components
4. **Error isolation**: Errors in workers don't crash the main application
5. **Performance**: Multiple workers can run in parallel (within browser limits)


---

## Limitations

WebAssembly in the browser is a capable environment, but it imposes real constraints that are particularly relevant for bioinformatics tools.

### Memory cap

`wasm32` modules have a maximum linear memory of **4 GiB**. This is a hard constraint for tools that load large reference genomes or full-resolution pangenome graphs into memory. Strategies include streaming and chunked processing, or restricting input sizes in the web interface. The forthcoming `wasm64` target would lift this restriction, and toolchain support is actively being developed (see [Available Targets](#available-targets)), but it has not yet landed in a stable release of the Rust/Web interfacing libraries we are considering in this guide.

### Call stack depth

The wasm call stack is much shallower than a typical native stack. Deeply recursive algorithms may overflow the stack at runtime. Rewrite critical paths iteratively or increase the stack size at link time (with a `.cargo/config.toml` entry such as `[target.wasm32-unknown-unknown] rustflags = ["-C", "link-args=-z stack-size=<bytes>"]`), accepting the trade-off of increased module size.

### Incompatible crates

Not every Rust crate compiles to `wasm32-unknown-unknown`. Crates that depend on OS-level facilities — file I/O, threads, system clocks, native dynamic linking, or C libraries via `cc`/`cmake` — will fail to compile. A notable example in bioinformatics is [`needletail`](https://crates.io/crates/needletail), which cannot currently be compiled to this target. When porting an existing tool, audit your dependency tree with:

```bash
cargo tree --target wasm32-unknown-unknown 2>&1 | grep -E "error|could not"
```

Common workarounds: replace the offending crate with a pure-Rust alternative; re-implement only the required subset; or gate the incompatible dependency behind `#[cfg(not(target_arch = "wasm32"))]` and providing a stub for the wasm build.

### SIMD and advanced CPU instructions

The WebAssembly SIMD proposal (`wasm32` SIMD128) is now widely supported in major browsers, but it is **not the same instruction set as x86 SSE/AVX or ARM NEON**. Code that relies on platform-specific SIMD intrinsics (e.g., via the `packed_simd` or `std::arch` crates targeting x86) will not compile. You can, however, use portable SIMD via the `std::simd` API (nightly) or crates such as `wide`, which will lower to SIMD128 when targeting wasm. Performance from SIMD128 is typically lower than from AVX2/AVX-512 on a desktop CPU.

### No multithreading (by default)

True threading in wasm requires `SharedArrayBuffer`, which is behind cross-origin isolation headers and not universally straightforward to configure. The `rayon` crate does not work out of the box in `wasm32-unknown-unknown`. For parallel workloads, consider Web Workers at the JavaScript level, passing disjoint chunks of data to separate wasm instances. Sharing data at this level is possible with shared memory objects, although it will apply the overall limit of 4MB.

### Performance expectations

WebAssembly can achieve 60–90 % of native performance for well-optimised, compute-bound code. However, memory-intensive code (with many small allocations or frequent cache misses) and code with a high wasm/JS boundary-crossing frequency can perform significantly worse. Profile the bottlenecks before assuming a native algorithm translates directly.

### No filesystem

`wasm32-unknown-unknown` has no access to the host filesystem. Input data must be passed from JavaScript as byte slices or strings (e.g., from a `<input type="file">` element read by the browser's `FileReader`/`File` API). This is a deliberate design constraint of the target, not a bug.

---

## A note on WebGPU

For computations that would benefit from GPU parallelism, [WebGPU](https://www.w3.org/TR/webgpu/) is the browser's emerging standard for general-purpose GPU computation (GPGPU). It is accessible from JavaScript and can be combined with a wasm module: the wasm code handles data preparation and result parsing, while dispatch and shader execution are managed through the WebGPU JavaScript API.

From Rust, the [`wgpu`](https://crates.io/crates/wgpu) crate provides a cross-platform WebGPU abstraction that compiles to `wasm32-unknown-unknown` (targeting the browser's WebGPU backend). This requires writing WGSL shaders and managing GPU buffer lifetimes, but it opens the door to genuinely GPU-accelerated methods in the browser. WebGPU is currently supported in Chrome and (behind a flag) in Firefox; Safari support is partial. Browser compatibility should be verified before committing to a GPU-dependent architecture.

---

## Minimum working example (MWE)

<a name="minimum-working-example-mwe"></a>

The following self-contained example computes the **GC content** of a DNA sequence entered by the user. It mirrors the structure used in this repository: a small Rust library compiled with `wasm-pack`, called from a single HTML file with plain JavaScript.

### Project layout

```
gc-wasm/
├── Cargo.toml
├── src/
│   └── lib.rs
└── www/
    ├── public/
    │   └── index.html
    ├── src/
    │   ├── components/
    │   │   └── GCContent.vue
    │   ├── workers/
    │   │   ├── gc.worker.ts
    │   │   └── GCWorker.ts
    │   ├── App.vue
    │   ├── main.ts
    │   └── shims-vue.d.ts
    ├── vue.config.js
    ├── package.json
    └── tsconfig.json
```

### `Cargo.toml`

```toml
[package]
name = "gc-wasm"
version = "0.1.0"
edition = "2021"

[lib]
crate-type = ["cdylib"]

[dependencies]
wasm-bindgen = "0.2"
```

### `src/lib.rs`

```rust
use wasm_bindgen::prelude::*;

/// Computes the GC content of a DNA sequence as a value in [0, 1].
/// Non-ACGT characters are counted in the denominator but not the numerator,
/// so ambiguous bases dilute the GC fraction rather than being ignored.
#[wasm_bindgen]
pub fn gc_content(sequence: &str) -> f64 {
    let total = sequence.len();
    if total == 0 {
        return 0.0;
    }
    let gc = sequence
        .bytes()
        .filter(|&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
        .count();
    gc as f64 / total as f64
}

/// Returns the length of the sequence after stripping whitespace and
/// newlines — useful when the user pastes a FASTA body.
#[wasm_bindgen]
pub fn clean_length(sequence: &str) -> usize {
    sequence
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .count()
}
```

### Vue integration

Instead of the manual build approach, integrate the WASM module with Vue using the wasm-pack-plugin:

#### Vue configuration (`vue.config.js`)

```javascript
const WasmPackPlugin = require("@wasm-tool/wasm-pack-plugin");
const path = require("path");

module.exports = {
    configureWebpack: {
        experiments: {
            asyncWebAssembly: true,
        },
    },
    chainWebpack: (config) => {
        config
            .plugin("wasm-pack_gc")
            .use(WasmPackPlugin)
            .init(
                (Plugin) =>
                    new Plugin({
                        crateDirectory: path.resolve(__dirname, "../gc-wasm"),
                        outDir: path.resolve(__dirname, "./src/pkg"),
                        forceMode: "production",
                    })
            )
            .end()
    },
}
```

#### Vue component (`src/components/GCContent.vue`)

```vue
<script lang="ts">
import { defineComponent, ref, onMounted, onUnmounted } from 'vue';
import { GCWorker } from '@/workers/GCWorker';

export default defineComponent({
  name: 'GCContent',
  setup() {
    const sequence = ref('');
    const result = ref('');
    const isLoading = ref(false);
    const gcWorker = ref<GCWorker | null>(null);
    
    // Initialize worker when component mounts
    onMounted(() => {
      gcWorker.value = new GCWorker();
    });
    
    // Clean up worker when component unmounts
    onUnmounted(() => {
      gcWorker.value?.terminate();
    });
    
    const calculate = async () => {
      if (!sequence.value.trim()) {
        result.value = 'Please enter a sequence';
        return;
      }
      
      if (!gcWorker.value) {
        result.value = 'Worker not initialized';
        return;
      }
      
      isLoading.value = true;
      try {
        // Clean sequence (remove whitespace and FASTA headers)
        const cleanSeq = sequence.value
          .split('\n')
          .filter(line => !line.startsWith('>'))
          .join('')
          .replace(/\s+/g, '');
        
        if (cleanSeq.length === 0) {
          result.value = 'No valid sequence found';
          return;
        }
        
        // Call WASM functions through worker
        const gc = await gcWorker.value.calculateGC(cleanSeq);
        const length = cleanSeq.length;
        
        result.value =
          `Length: ${length} bp\n` +
          `GC: ${(gc * 100).toFixed(2)}%`;
          
      } catch (error) {
        console.error('WASM worker error:', error);
        result.value = `Error: ${error.message}`;
      } finally {
        isLoading.value = false;
      }
    };
    
    return { sequence, result, isLoading, calculate };
  }
});
</script>

<template>
  <div class="gc-calculator">
    <h2>GC Content Calculator</h2>
    <p>Paste a DNA sequence (FASTA body or raw bases):</p>
    <textarea 
      v-model="sequence" 
      placeholder="ATGCGCATGCTTAAGC... or FASTA format"
      rows="6"
      class="w-full p-2 border rounded"
    ></textarea>
    <button 
      @click="calculate" 
      :disabled="isLoading"
      class="mt-2 px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600 disabled:bg-gray-300"
    >
      {{ isLoading ? 'Calculating...' : 'Calculate' }}
    </button>
    <div class="mt-4 p-2 bg-gray-100 rounded font-mono whitespace-pre">
      {{ result }}
    </div>
  </div>
</template>

<style scoped>
.gc-calculator {
  max-width: 600px;
  margin: 0 auto;
  font-family: sans-serif;
}
textarea {
  font-family: monospace;
  width: 100%;
}
</style>
```

#### Worker files

##### Worker implementation (`src/workers/gc.worker.ts`)

```typescript
// Handle messages from main thread
self.onmessage = async (event) => {
  try {
    // Dynamically import the WASM module
    const wasm = await import('@/pkg');
    
    // Call the WASM function
    const result = wasm.gc_content(event.data.sequence);
    
    // Send result back to main thread
    self.postMessage({ result });
  } catch (error) {
    console.error('WASM error in worker:', error);
    self.postMessage({ error: error.message });
  }
};
```

##### Worker class (`src/workers/GCWorker.ts`)

```typescript
// Worker class that manages the web worker lifecycle
export class GCWorker {
  private worker: Worker;
  
  constructor() {
    // Create worker instance
    this.worker = new Worker(new URL('./gc.worker.ts', import.meta.url));
  }
  
  public calculateGC(sequence: string): Promise<number> {
    return new Promise((resolve, reject) => {
      // Generate unique message ID for this request
      const messageId = Math.random().toString(36).substring(2);
      
      // Message handler for this specific request
      const handler = (event: MessageEvent) => {
        if (event.data.messageId === messageId) {
          this.worker.removeEventListener('message', handler);
          if (event.data.result !== undefined) {
            resolve(event.data.result);
          } else if (event.data.error) {
            reject(new Error(event.data.error));
          }
        }
      };
      
      this.worker.addEventListener('message', handler);
      
      // Send sequence to worker
      this.worker.postMessage({
        messageId,
        sequence
      });
    });
  }
  
  public terminate() {
    // Clean up worker when no longer needed
    this.worker.terminate();
  }
}
```

### `www/index.html`
=======

Place this file inside the `www/` directory, and serve both `www/` and `pkg/` from the same origin (e.g., with `python3 -m http.server` from the `gc-wasm/` directory, adjusting the import path accordingly, or by copying `pkg/` into `www/pkg/`).

```html
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>GC Content Calculator</title>
  <style>
    body      { font-family: monospace; max-width: 640px; margin: 2rem auto; }
    textarea  { width: 100%; height: 8rem; font-family: monospace; }
    output    { display: block; margin-top: 1rem; font-size: 1.2rem; }
    .error    { color: #c00; }
  </style>
</head>
<body>
  <h1>GC Content Calculator</h1>
  <p>Paste a DNA sequence (FASTA body or raw bases):</p>
  <textarea id="seq" placeholder="ATGCGCATGCTTAAGC..."></textarea>
  <button id="calc">Calculate</button>
  <output id="result"></output>

  <script type="module">
    // Adjust the path if pkg/ lives elsewhere relative to this file.
    import init, { gc_content, clean_length } from "../pkg/gc_wasm.js";

    // Initialise the wasm module once, before any calls.
    await init();

    document.getElementById("calc").addEventListener("click", () => {
      const raw   = document.getElementById("seq").value;
      const out   = document.getElementById("result");

      if (!raw.trim()) {
        out.textContent = "";
        return;
      }

      // Strip FASTA header lines if present.
      const body = raw
        .split("\n")
        .filter(line => !line.startsWith(">"))
        .join("");

      const len = clean_length(body);
      const gc  = gc_content(body);

      out.textContent =
        `Length : ${len} bp\n` +
        `GC     : ${(gc * 100).toFixed(2)} %`;
    });
  </script>
</body>
</html>
```

### Serving locally

```bash
# Navigate to the Vue project directory
cd www

# Install dependencies
npm install

# Start development server
npm run serve

# The Vue CLI development server will start and show the URL
# Typically http://localhost:8080
# Open this URL in your browser
```

> **Note**: The Vue CLI development server handles all the necessary configuration for WASM module loading. Unlike the Python HTTP server approach, it provides hot module replacement, proper asset handling, and automatic reloading during development.

### Building for production

```bash
# Create optimized production build
npm run build

# The built files will be in www/dist/ directory
# You can serve these with any static file server
```

### `www/package.json`

```json
{
  "name": "gc-wasm-vue",
  "version": "0.1.0",
  "private": true,
  "scripts": {
    "serve": "vue-cli-service serve",
    "build": "vue-cli-service build",
    "lint": "vue-cli-service lint"
  },
  "dependencies": {
    "vue": "^3.2.13",
    "core-js": "^3.8.3"
  },
  "devDependencies": {
    "@vue/cli-service": "^5.0.0",
    "@wasm-tool/wasm-pack-plugin": "^1.0.0",
    "typescript": "^4.5.0",
    "@vue/compiler-sfc": "^3.2.13"
  }
}
```

---

## Further Reading

- [The `wasm-bindgen` Guide](https://rustwasm.github.io/docs/wasm-bindgen/)
- [Rust and WebAssembly Book](https://rustwasm.github.io/docs/book/)
- [`wasm-pack` documentation](https://rustwasm.github.io/wasm-pack/book/)
- [WebAssembly specification](https://webassembly.github.io/spec/)
- [WebGPU specification](https://www.w3.org/TR/webgpu/)
- [`wgpu` crate](https://crates.io/crates/wgpu)
