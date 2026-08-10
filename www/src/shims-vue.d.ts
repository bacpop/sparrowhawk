/* eslint-disable */
declare module '*.vue' {
    import type {DefineComponent} from 'vue'
    const component: DefineComponent<{}, {}, any>
    export default component
}

declare module 'taxonium-component';
declare module 'mgnify-jbrowse';

// wasm-bindgen module namespace (webpack asyncWebAssembly); used to read
// the linear-memory size for the peak-memory summaries.
declare module '*/index_bg.wasm' {
    export const memory: WebAssembly.Memory;
}
