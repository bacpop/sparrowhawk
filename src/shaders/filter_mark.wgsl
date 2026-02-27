// Filter-mark kernel.
// output_mask[i] = 1 if boundary[i]==1 AND seg_count[segment_id[i]] >= min_count, else 0.

@group(0) @binding(0) var<storage, read>       boundary:    array<u32>;
@group(0) @binding(1) var<storage, read>       segment_id:  array<u32>;
@group(0) @binding(2) var<storage, read>       seg_count:   array<u32>;
@group(0) @binding(3) var<storage, read_write> output_mask: array<u32>;
@group(0) @binding(4) var<uniform>             min_count:   u32;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i >= arrayLength(&boundary) { return; }

    if boundary[i] == 1u && seg_count[segment_id[i]] >= min_count {
        output_mask[i] = 1u;
    } else {
        output_mask[i] = 0u;
    }
}
