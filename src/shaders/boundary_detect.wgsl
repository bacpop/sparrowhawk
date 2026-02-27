// Boundary detection kernel.
// Thread i sets boundary[i] = 1 if it starts a new kmer run, else 0.

struct KmerEntry {
    canon_lo: u32,
    canon_hi: u32,
    nc_lo:    u32,
    nc_hi:    u32,
    bases:    u32,
    pad:      u32,
}

@group(0) @binding(0) var<storage, read>       elements: array<KmerEntry>;
@group(0) @binding(1) var<storage, read_write> boundary: array<u32>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = arrayLength(&elements);
    if i >= n { return; }

    if i == 0u {
        boundary[i] = 1u;
    } else {
        let cur  = elements[i];
        let prev = elements[i - 1u];
        if cur.canon_lo != prev.canon_lo || cur.canon_hi != prev.canon_hi {
            boundary[i] = 1u;
        } else {
            boundary[i] = 0u;
        }
    }
}
