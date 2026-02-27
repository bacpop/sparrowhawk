// Per-workgroup offset computation.
// For each digit d, computes the exclusive prefix sum of wg_hist[*][d]
// across all workgroups and writes the result to wg_prefix[wg * 256 + d].
//
// This is a column-wise scan of the [num_wg x 256] wg_hist matrix:
//   wg_prefix[wg * 256 + d] = sum_{wg' < wg} wg_hist[wg' * 256 + d]
//
// Dispatched as 1 workgroup of 256 threads.
// Thread d handles the entire column for digit value d sequentially.

@group(0) @binding(0) var<storage, read>       wg_hist:   array<u32>;
@group(0) @binding(1) var<storage, read_write> wg_prefix: array<u32>;
@group(0) @binding(2) var<uniform>             num_wg:    u32;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3<u32>) {
    let d = lid.x; // this thread is responsible for digit value d
    var running = 0u;
    for (var wg = 0u; wg < num_wg; wg++) {
        let count = wg_hist[wg * 256u + d];
        wg_prefix[wg * 256u + d] = running;
        running += count;
    }
}
