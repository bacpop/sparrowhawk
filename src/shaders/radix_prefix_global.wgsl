// Exclusive Blelloch scan of global_hist[256] → global_prefix[256].
// Single workgroup, 256 threads.

@group(0) @binding(0) var<storage, read>       global_hist:   array<u32>;
@group(0) @binding(1) var<storage, read_write> global_prefix: array<u32>;

var<workgroup> scan_data: array<u32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3<u32>) {
    scan_data[lid.x] = global_hist[lid.x];
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

    global_prefix[lid.x] = scan_data[lid.x];
}
