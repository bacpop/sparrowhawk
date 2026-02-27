// Radix sort histogram kernel.
// Each workgroup accumulates a local histogram in workgroup memory, then:
// - Atomically adds to global_hist[256] (one slot per thread)
// - Writes its local sub-histogram to wg_hist[wg_id * 256 + d] (non-atomic, unique slots)

struct KmerEntry {
    canon_lo: u32,
    canon_hi: u32,
    nc_lo:    u32,
    nc_hi:    u32,
    bases:    u32,
    pad:      u32,
}

@group(0) @binding(0) var<storage, read>       elements:    array<KmerEntry>;
@group(0) @binding(1) var<storage, read_write> global_hist: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> wg_hist:     array<u32>;
@group(0) @binding(3) var<uniform>             pass_idx:    u32;

var<workgroup> local_hist: array<atomic<u32>, 256>;

fn digit(e: KmerEntry, p: u32) -> u32 {
    if p < 4u {
        return (e.canon_lo >> (p * 8u)) & 0xFFu;
    } else {
        return (e.canon_hi >> ((p - 4u) * 8u)) & 0xFFu;
    }
}

@compute @workgroup_size(256)
fn main(
    @builtin(global_invocation_id) gid:  vec3<u32>,
    @builtin(local_invocation_id)  lid:  vec3<u32>,
    @builtin(workgroup_id)         wgid: vec3<u32>,
) {
    // Each thread clears its own slot in the local histogram
    atomicStore(&local_hist[lid.x], 0u);
    workgroupBarrier();

    let n = arrayLength(&elements);
    if gid.x < n {
        let d = digit(elements[gid.x], pass_idx);
        atomicAdd(&local_hist[d], 1u);
    }

    workgroupBarrier();

    // Thread lid.x handles digit value lid.x:
    // - accumulate into global histogram (atomic, across workgroups)
    // - write this workgroup's count to wg_hist (non-atomic, unique slot)
    let local_count = atomicLoad(&local_hist[lid.x]);
    atomicAdd(&global_hist[lid.x], local_count);
    wg_hist[wgid.x * 256u + lid.x] = local_count;
}
