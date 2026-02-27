// For each workgroup wg, compute exclusive prefix sum of
// wg_hist[wg*256 .. wg*256+256] → wg_prefix[wg*256 .. wg*256+256].
// Dispatch num_wg workgroups of 256 threads.

@group(0) @binding(0) var<storage, read>       wg_hist:   array<u32>;
@group(0) @binding(1) var<storage, read_write> wg_prefix: array<u32>;

var<workgroup> scan_data: array<u32, 256>;

@compute @workgroup_size(256)
fn main(
    @builtin(local_invocation_id) lid:  vec3<u32>,
    @builtin(workgroup_id)        wgid: vec3<u32>,
) {
    let base = wgid.x * 256u;
    scan_data[lid.x] = wg_hist[base + lid.x];
    workgroupBarrier();

    // Up-sweep
    var offset = 1u;
    for (var d = 128u; d > 0u; d >>= 1u) {
        workgroupBarrier();
        if lid.x < d {
            let ai = offset * (2u * lid.x + 1u) - 1u;
            let bi = offset * (2u * lid.x + 2u) - 1u;
            scan_data[bi] += scan_data[ai];
        }
        offset <<= 1u;
    }

    if lid.x == 0u {
        scan_data[255] = 0u;
    }
    workgroupBarrier();

    // Down-sweep
    for (var d = 1u; d < 256u; d <<= 1u) {
        offset >>= 1u;
        workgroupBarrier();
        if lid.x < d {
            let ai = offset * (2u * lid.x + 1u) - 1u;
            let bi = offset * (2u * lid.x + 2u) - 1u;
            let t = scan_data[ai];
            scan_data[ai] = scan_data[bi];
            scan_data[bi] += t;
        }
    }
    workgroupBarrier();

    wg_prefix[base + lid.x] = scan_data[lid.x];
}
