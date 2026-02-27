// Scatter-output kernel.
// For each i where output_mask[i] == 1:
//   output[output_index[i]] = KmerOut { canon, nc, bases, count }

struct KmerEntry {
    canon_lo: u32,
    canon_hi: u32,
    nc_lo:    u32,
    nc_hi:    u32,
    bases:    u32,
    pad:      u32,
}

struct KmerOut {
    canon_lo: u32,
    canon_hi: u32,
    nc_lo:    u32,
    nc_hi:    u32,
    bases:    u32,
    count:    u32,
}

@group(0) @binding(0) var<storage, read>       elements:     array<KmerEntry>;
@group(0) @binding(1) var<storage, read>       output_mask:  array<u32>;
@group(0) @binding(2) var<storage, read>       output_index: array<u32>;
@group(0) @binding(3) var<storage, read>       segment_id:   array<u32>;
@group(0) @binding(4) var<storage, read>       seg_count:    array<u32>;
@group(0) @binding(5) var<storage, read_write> output:       array<KmerOut>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i >= arrayLength(&elements) { return; }

    if output_mask[i] == 1u {
        let e = elements[i];
        let out_pos = output_index[i];
        output[out_pos] = KmerOut(
            e.canon_lo,
            e.canon_hi,
            e.nc_lo,
            e.nc_hi,
            e.bases,
            seg_count[segment_id[i]],
        );
    }
}
