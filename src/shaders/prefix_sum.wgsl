// Multi-pass prefix sum for arbitrary-length arrays.
// Three entry points for a 3-pass Blelloch scan:
//   local_scan:  each workgroup scans its 256-element tile; writes tile sum to block_sums[]
//   block_scan:  single workgroup scans block_sums[] (up to 256 blocks)
//   propagate:   add block_sums[workgroup_id] to each element in tile

@group(0) @binding(0) var<storage, read_write> data:       array<u32>;
@group(0) @binding(1) var<storage, read_write> block_sums: array<u32>;

var<workgroup> tile: array<u32, 256>;

// Pass A: local tile scan (exclusive), write tile total to block_sums[wg_id]
@compute @workgroup_size(256)
fn local_scan(
    @builtin(global_invocation_id) gid:  vec3<u32>,
    @builtin(local_invocation_id)  lid:  vec3<u32>,
    @builtin(workgroup_id)         wgid: vec3<u32>,
) {
    let n = arrayLength(&data);
    if gid.x < n {
        tile[lid.x] = data[gid.x];
    } else {
        tile[lid.x] = 0u;
    }
    workgroupBarrier();

    // Up-sweep
    var offset = 1u;
    for (var d = 128u; d > 0u; d >>= 1u) {
        workgroupBarrier();
        if lid.x < d {
            let ai = offset * (2u * lid.x + 1u) - 1u;
            let bi = offset * (2u * lid.x + 2u) - 1u;
            tile[bi] += tile[ai];
        }
        offset <<= 1u;
    }

    if lid.x == 0u {
        block_sums[wgid.x] = tile[255];
        tile[255] = 0u;
    }
    workgroupBarrier();

    // Down-sweep
    for (var d = 1u; d < 256u; d <<= 1u) {
        offset >>= 1u;
        workgroupBarrier();
        if lid.x < d {
            let ai = offset * (2u * lid.x + 1u) - 1u;
            let bi = offset * (2u * lid.x + 2u) - 1u;
            let t = tile[ai];
            tile[ai] = tile[bi];
            tile[bi] += t;
        }
    }
    workgroupBarrier();

    if gid.x < n {
        data[gid.x] = tile[lid.x];
    }
}

// Pass B: exclusive prefix scan of block_sums[] (single workgroup).
// Thread 0 does a sequential scan so correctness is not limited by tile size.
// The Blelloch approach above only handles nb ≤ 256; for large buckets
// (nb ≈ 14 000) we need the full sequential pass.
@compute @workgroup_size(256)
fn block_scan(
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    if lid.x == 0u {
        let nb = arrayLength(&block_sums);
        var running = 0u;
        for (var i = 0u; i < nb; i++) {
            let v = block_sums[i];
            block_sums[i] = running;
            running += v;
        }
    }
}

// Pass C: add block_sums[wg_id] to each element in tile
@compute @workgroup_size(256)
fn propagate(
    @builtin(global_invocation_id) gid:  vec3<u32>,
    @builtin(workgroup_id)         wgid: vec3<u32>,
) {
    let n = arrayLength(&data);
    if gid.x < n {
        data[gid.x] += block_sums[wgid.x];
    }
}
