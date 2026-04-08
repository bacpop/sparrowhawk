# Sparrowhawk-Web Repository Documentation

## Repository Overview

**Sparrowhawk-Web** is a web-based bioinformatics platform that provides genome assembly, sequence analysis, and taxonomic identification capabilities directly in the browser using WebAssembly technology.

### Purpose and Scope
- Web interface for the Sparrowhawk genome assembler
- Browser-based bioinformatics tools using Rust compiled to WebAssembly
- Local processing without server requirements (all computation happens in-browser)
- Support for various bioinformatics workflows: assembly, mapping, alignment, taxonomic ID, gene calling, and host depletion

### Key Features
- **Genome Assembly**: WebAssembly-compiled Sparrowhawk assembler
- **Sequence Mapping/Alignment**: Within-species comparison using ska.rust
- **Taxonomic Identification**: Species-level classification using sketchlib.rust
- **Gene Calling**: Prodigal-based gene prediction via Orphos
- **Host Depletion**: Read filtering using Deacon
- **Interactive Visualization**: MSA viewer, phylogenetic trees, genome browser

### Target Audience
- Microbiologists and bioinformaticians
- Researchers working with bacterial genomes
- Users needing local, private data analysis
- Low-resource contexts (works on smartphones/tablets)

## Architecture Overview

```
┌───────────────────────────────────────────────────────┐
│                 Sparrowhawk-Web Architecture            │
├─────────────────┬─────────────────┬───────────────────┤
│   Frontend     │   WebAssembly  │    Rust Crates    │
│  (Vue.js)       │    Modules      │   (Compiled)      │
├─────────────────┼─────────────────┼───────────────────┤
│ - Vue 3         │ - @/pkg         │ - sparrowhawk    │
│ - Vuex          │ - @/pkg_ska     │ - ska.rust       │
│ - Workers       │ - @/pkg_orphos │ - sketchlib.rust │
│ - Components    │ - @/pkg_sketchlib│ - orphos-bridge │
│ - UI Framework  │ - @/pkg_deacon  │ - deacon-bridge  │
└─────────────────┴─────────────────┴───────────────────┘
```

### Frontend-Backend Separation
The application follows a clear separation:
- **Frontend**: Vue.js components, state management, UI
- **Backend**: Rust crates compiled to WebAssembly for computational tasks
- **Bridge**: Worker-based communication between frontend and WASM modules

### WebAssembly Integration Points
- Workers handle heavy computation in separate threads
- Vuex stores manage application state
- Components trigger workers and display results
- WASM modules expose JavaScript-compatible interfaces

## Frontend Structure

### Vue.js Component Hierarchy

```
App.vue
├── Sidebar (Navigation)
└── Main Content Area
    ├── AssemblyPage
    ├── MappingAlignmentPage
    ├── TaxonomicIDPage
    ├── GeneCallingPage
    ├── HostDepletionPage
    └── FaqPage
```

### Key Components

#### Core Components
- **App.vue**: Main application container with sidebar navigation
- **AssemblyPage.vue**: Genome assembly interface
- **MappingAlignmentPage.vue**: Mapping and alignment tools
- **TaxonomicIDPage.vue**: Taxonomic identification interface
- **GeneCallingPage.vue**: Gene prediction interface
- **HostDepletionPage.vue**: Host read filtering interface

#### UI Components
- **DownloadButton.vue**: Result download functionality
- **KmerHistogram.vue**: k-mer frequency visualization
- **SequenceViewer/**: MSA visualization components
- **MSAViewer.vue**: Multiple sequence alignment viewer
- **help/**: Contextual help components for each page

#### Workers
- **Assembler.worker.ts**: Genome assembly worker
- **Mapper.worker.ts**: Mapping worker (ska.rust)
- **Sketcher.worker.ts**: Taxonomic ID worker (sketchlib)
- **Caller.worker.ts**: Gene calling worker (orphos)
- **Depleter.worker.ts**: Host depletion worker (deacon)

### State Management with Vuex

The application uses Vuex for centralized state management:

**Store Modules:**
- **state.ts**: Root state definition
- **actions.ts**: Async operations and worker coordination
- **mutations.ts**: State modifications
- **getters.ts**: Computed state properties

**Key State Elements:**
- `workerState`: Manages WebWorker instances
- `readsPreprocessing`: Stores preprocessing results
- `allResults`: Assembly outputs and visualizations
- `allResults_ska`: Mapping/alignment results
- `allResults_sketchlib`: Taxonomic ID results
- `allResults_orphos`: Gene calling results
- `allResults_deacon`: Host depletion results

### Worker-Based Architecture

Each computational module runs in a dedicated WebWorker:

```
Frontend → Worker → WASM Module → Results → Frontend
```

**Worker Communication Pattern:**
1. Frontend component dispatches action
2. Action sends message to appropriate worker
3. Worker loads WASM module and executes computation
4. Worker posts results back to frontend
5. Frontend updates state and UI

## Rust/WASM Components

### Sparrowhawk Assembler

**Location**: `rust/sparrowhawk/`

**Purpose**: Genome assembly using de Bruijn graphs

**Key Features:**
- k-mer based assembly algorithm
- Preprocessing with optional Bloom filters
- Graph construction and simplification
- Contig generation in FASTA format
- Graph export (DOT, GFAv1, GFAv2 formats)

**WASM Integration:**
- Compiled with `wasm-pack`
- Exposed as `@/pkg` in frontend
- Provides `AssemblyHelper` class for preprocessing and assembly

### ska.rust

**Location**: `rust/ska.rust/`

**Purpose**: Sequence mapping and alignment

**Key Features:**
- k-mer based sequence comparison
- Within-species distance calculation
- Multiple sequence alignment
- Phylogenetic tree construction
- Reference-based mapping

**WASM Integration:**
- Compiled as `@/pkg_ska`
- Used for both mapping and alignment workflows
- Provides reference indexing and query processing

### sketchlib.rust

**Location**: `rust/sketchlib.rust/`

**Purpose**: Taxonomic identification

**Key Features:**
- k-mer sketching for genomic sequences
- Reference database comparison
- Species-level classification
- Probability-based identification
- Metadata integration (Gemsparcl, GTDB)

**WASM Integration:**
- Compiled as `@/pkg_sketchlib`
- Parallel worker pool for multi-file processing
- Returns taxonomic probabilities and metadata

### orphos (Prodigal Rust Port)

**Location**: `rust/orphos-bridge/`

**Purpose**: Gene calling and annotation

**Key Features:**
- Prodigal algorithm implementation
- Open reading frame detection
- Start codon identification
- Gene feature prediction
- GFF output generation
- Genome visualization support

**WASM Integration:**
- Compiled as `@/pkg_orphos-bridge`
- Parallel processing for multiple genomes
- Generates indexed GFF and FASTA outputs

### deacon

**Location**: `rust/deacon-bridge/`

**Purpose**: Host read depletion/filtering

**Key Features:**
- k-mer based read classification
- Index-based filtering
- Absolute and relative thresholding
- Host contamination removal
- Read preservation or depletion modes

**WASM Integration:**
- Compiled as `@/pkg_deacon`
- Single worker for sequential processing
- Handles paired-end read preservation

### Compilation and Integration Process

1. **Rust Toolchain Setup**:
   ```bash
   rustup target add wasm32-unknown-unknown
   cargo install wasm-pack
   ```

2. **WASM Compilation**:
   - Each crate has its own `Cargo.toml` with WASM features
   - Compiled using `wasm-pack build --target web`
   - Outputs JavaScript/TypeScript bindings

3. **Frontend Integration**:
   - Import statements in worker files
   - Type definitions in `workers.d.ts`
   - Dynamic imports in main application

4. **Build Process**:
   - `npm run serve` triggers WASM compilation
   - Webpack handles WASM module bundling
   - Workers load WASM modules on demand

## Build System

### NPM Dependencies and Scripts

**Key Dependencies:**
- `@vue/cli-service`: Vue CLI for development
- `vuex`: State management
- `vue3-dropzone`: File upload handling
- `fflate`: Compression utilities
- `plotly.js-dist`: Data visualization
- `d3`: SVG manipulation
- `taxonium-component`: Phylogenetic tree viewer
- `mgnify-jbrowse`: Genome browser
- `reka-ui`: UI component library

**Scripts:**
- `serve`: Development server with hot reload
- `build`: Production build
- `lint`: Code linting

### Rust Toolchain Requirements

- Rust 1.70+
- `wasm32-unknown-unknown` target
- `wasm-pack` for WASM compilation
- `wasm-bindgen` for JS interop

### WASM Compilation Process

Each Rust crate follows this compilation pattern:

```bash
# In each rust crate directory
wasm-pack build --target web --out-dir ../www/pkg_crate_name
```

The frontend imports these compiled modules:
- `@/pkg` - sparrowhawk
- `@/pkg_ska` - ska.rust
- `@/pkg_sketchlib` - sketchlib
- `@/pkg_orphos-bridge` - orphos
- `@/pkg_deacon` - deacon

### Development vs Production Builds

**Development:**
- Source maps enabled
- Hot module replacement
- Debug logging
- Larger bundle sizes

**Production:**
- Optimized WASM modules
- Tree shaking
- Minification
- Smaller bundle sizes

## Key Workflows

### Assembly Workflow

```
User Upload → Preprocessing → Assembly → Results Display
```

1. **File Upload**: FASTQ files (single or paired-end)
2. **Preprocessing**:
   - k-mer extraction
   - Quality filtering
   - Count thresholding
   - Optional Bloom filter
3. **Assembly**:
   - Graph construction
   - Bubble collapse
   - Dead-end removal
   - Contig generation
4. **Results**:
   - FASTA download
   - Graph visualization (DOT/GFA)
   - k-mer histogram

### Mapping/Alignment Workflow

```
Reference Upload → Query Upload → Processing → Visualization
```

1. **Reference Indexing**:
   - FASTA file upload
   - k-mer indexing
   - Canonical k-mer generation
2. **Query Processing**:
   - Read mapping
   - Variant detection
   - Coverage calculation
3. **Alignment**:
   - Multiple sequence alignment
   - Phylogenetic tree construction
4. **Visualization**:
   - MSA viewer
   - Taxonium tree viewer
   - Alignment downloads

### Taxonomic Identification Workflow

```
Sample Upload → Sketching → Comparison → Results
```

1. **File Processing**:
   - FASTQ/FASTA upload
   - Subsampling (if needed)
   - Quality filtering
2. **Sketching**:
   - k-mer extraction
   - MinHash sketch generation
3. **Comparison**:
   - Reference database matching
   - Probability calculation
4. **Results**:
   - Species identification
   - Probability scores
   - Metadata display
   - TSV export

### Gene Calling Workflow

```
Genome Upload → ORF Detection → Annotation → Visualization
```

1. **File Upload**:
   - FASTA genome files
   - Multiple file support
2. **ORF Detection**:
   - Start codon identification
   - Reading frame analysis
   - Gene boundary detection
3. **Annotation**:
   - GFF file generation
   - Feature indexing
4. **Visualization**:
   - JBrowse genome browser
   - Interactive feature exploration
   - Bulk download options

### Host Depletion Workflow

```
Index Upload → Reads Upload → Filtering → Results
```

1. **Index Loading**:
   - Deacon index file (.idx)
   - Host reference k-mers
2. **Reads Processing**:
   - FASTQ file upload
   - k-mer matching
   - Threshold application
3. **Filtering**:
   - Host read removal/preservation
   - Paired-end handling
4. **Results**:
   - Filtered FASTQ downloads
   - Statistics display
   - Bulk export options

## Testing and Deployment

### Available Test Suites

The repository includes:
- **Vue component tests**: Unit tests for UI components
- **Worker tests**: Integration tests for worker-WASM communication
- **Rust tests**: Unit tests for Rust crates

**Running Tests:**
```bash
# Frontend tests
npm test

# Rust tests (in each crate directory)
cargo test
```

### Deployment Considerations

**Memory Constraints:**
- 4GB RAM limit due to 32-bit WASM memory addressing
- Chunking options for large datasets
- Bloom filter alternatives for memory efficiency
- Progress indicators for long-running operations

**Browser Compatibility:**
- Modern browsers with WebAssembly support
- WebWorker support required
- Minimum Chrome 57+, Firefox 52+, Safari 11+, Edge 16+

**Performance Optimization:**
- Worker pooling for parallel processing
- Virtual scrolling for large datasets
- Lazy loading of components
- WASM module caching

### Memory Constraints and Limitations

**Current Limitations:**
- 4GB total memory (32-bit addressing)
- Large read files may fail without chunking
- Browser tab crashes possible with excessive memory usage

**Mitigation Strategies:**
- Chunked preprocessing (`csize` parameter)
- Bloom filter mode (`do_bloom` parameter)
- Automatic k-mer count fitting (`do_fit` parameter)
- Progress monitoring and error handling

## Common Patterns

### Worker Communication Patterns

**Request-Response Pattern:**
```javascript
// Frontend → Worker
worker.postMessage({ action: 'process', data: files });

// Worker → Frontend
worker.onmessage = (msg) => {
  if (msg.data.results) {
    // Update state
  }
  if (msg.data.error) {
    // Handle error
  }
};
```

**Progress Reporting:**
```javascript
// Worker posts progress updates
worker.postMessage({ 
  progress: current/total, 
  status: 'Processing reads...'
});
```

### State Management Patterns

**Centralized State Updates:**
```javascript
// In actions.ts
async processReads({ commit }, payload) {
  commit('setProcessingState', true);
  try {
    // Worker processing
    commit('setResults', results);
  } catch (error) {
    commit('setError', error);
  } finally {
    commit('setProcessingState', false);
  }
}
```

**Reactive UI Updates:**
```vue
<!-- In components -->
<template>
  <div v-if="isProcessing">Loading...</div>
  <div v-else-if="hasResults">Results</div>
</template>

<script>
export default {
  computed: {
    isProcessing() { return this.$store.getters.isProcessing; },
    hasResults() { return this.$store.getters.hasResults; }
  }
}
</script>
```

### File Handling and Processing

**File Validation:**
```javascript
// In utils.ts
const validExtensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz'];
const isValidFile = (file) => {
  return validExtensions.some(ext => file.name.endsWith(ext));
};
```

**Paired-End Handling:**
```javascript
// In utils.ts
export function getFilesToProcess(files) {
  // Logic to pair _1 and _2 files
  // Return array of file groups
}
```

### Error Handling Strategies

**Memory Error Handling:**
```javascript
// In Assembler.ts
try {
  await helper.preprocess(file1, file2);
} catch (error) {
  console.error('Memory issue:', error);
  worker.postMessage({ reset: true });
  return;
}
```

**User Feedback:**
```vue
<!-- Error display -->
<div v-if="errorType === 'memory'" class="error-message">
  Memory error! Try chunking or Bloom filter.
</div>
```

## Development Setup

### Prerequisites

1. **Node.js and NPM**:
   ```bash
   # Install Node.js 16+ and NPM
   ```

2. **Rust Toolchain**:
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   rustup target add wasm32-unknown-unknown
   cargo install wasm-pack
   ```

3. **Project Setup**:
   ```bash
   git clone --recurse-submodules https://github.com/bacpop/sparrowhawk-web.git
   cd sparrowhawk-web/www
   npm install
   ```

### Running the Development Server

```bash
cd www
npm run serve
```

This will:
1. Compile all Rust crates to WASM
2. Start Vue development server
3. Enable hot module replacement
4. Open browser at localhost:8080

### Building for Production

```bash
cd www
npm run build
```

Production build includes:
- Optimized WASM modules
- Minified JavaScript
- Tree-shaken dependencies
- Output in `dist/` directory

## Repository Structure

```
sparrowhawk-web/
├── .vibe/                  # VIBE documentation
├── rust/                  # Rust crates
│   ├── sparrowhawk/      # Genome assembler
│   ├── ska.rust/          # Mapping/alignment
│   ├── sketchlib.rust/    # Taxonomic ID
│   ├── orphos-bridge/     # Gene calling
│   └── deacon-bridge/     # Host depletion
├── www/                   # Frontend code
│   ├── public/            # Static assets
│   ├── src/               # Vue.js source
│   │   ├── components/    # UI components
│   │   ├── workers/       # Web Workers
│   │   ├── store/        # Vuex store
│   │   └── ...           # Other frontend code
│   ├── package.json      # NPM dependencies
│   └── vue.config.js     # Vue configuration
├── README.md              # Main documentation
└── LICENSE                # License information
```

## Key Files Reference

### Frontend Entry Points
- `www/src/main.ts`: Application entry point
- `www/src/App.vue`: Main application component
- `www/src/store/index.ts`: Vuex store setup

### Worker Files
- `www/src/workers/Assembler.worker.ts`: Assembly worker
- `www/src/workers/Mapper.worker.ts`: Mapping worker
- `www/src/workers/Sketcher.worker.ts`: Taxonomic ID worker
- `www/src/workers/Caller.worker.ts`: Gene calling worker
- `www/src/workers/Depleter.worker.ts`: Host depletion worker

### Component Files
- `www/src/components/pages/`: Main page components
- `www/src/components/ui/`: UI component library
- `www/src/components/SequenceViewer/`: MSA visualization
- `www/src/components/help/`: Help components

### Rust Crate Entry Points
- `rust/sparrowhawk/src/lib.rs`: Sparrowhawk WASM interface
- `rust/ska.rust/src/lib.rs`: ska.rust WASM interface
- `rust/sketchlib.rust/src/lib.rs`: sketchlib WASM interface
- `rust/orphos-bridge/src/lib.rs`: Orphos WASM interface
- `rust/deacon-bridge/src/lib.rs`: Deacon WASM interface

## Troubleshooting

### Common Issues

**WASM Compilation Errors:**
- Ensure Rust toolchain is properly installed
- Check `wasm32-unknown-unknown` target is available
- Verify `wasm-pack` version compatibility

**Memory Errors:**
- Reduce input file sizes
- Enable chunking (`csize` parameter)
- Use Bloom filter mode (`do_bloom`)
- Close other browser tabs

**Worker Communication Issues:**
- Check worker instantiation
- Verify message formats
- Ensure proper error handling

**Browser Compatibility:**
- Use Chrome/Firefox for best results
- Ensure WebAssembly support
- Check console for detailed errors

### Debugging Tips

**Frontend Debugging:**
```bash
# Run with debug logging
npm run serve
# Open browser dev tools
# Check console and network tabs
```

**Rust Debugging:**
```bash
# Build with debug symbols
wasm-pack build --debug
# Use console_error_panic_hook in Rust
```

**Worker Debugging:**
```javascript
// Add logging in workers
console.log('Worker received:', message);
// Check worker console output
```
