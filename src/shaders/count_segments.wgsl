// Count-segments kernel.
// For each element i, atomicAdd to the correct 0-based segment index.
//
// segment_id[i] = exclusive prefix sum of boundary[].
// The actual segment index is segment_id[i] + boundary[i] - 1:
//   - When boundary[i]=1 (start of a new run): formula = segment_id[i] (correct).
//   - When boundary[i]=0 (continuation):       exclusive prefix at i overcounts by 1
//     because boundary[0]=1 always; subtracting 1 corrects this.
// This produces 0-based segment indices matching those expected by filter_mark.

@group(0) @binding(0) var<storage, read>       segment_id: array<u32>;
@group(0) @binding(1) var<storage, read_write> seg_count:  array<atomic<u32>>;
@group(0) @binding(2) var<storage, read>       boundary:   array<u32>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i >= arrayLength(&segment_id) { return; }
    let seg = segment_id[i] + boundary[i] - 1u;
    atomicAdd(&seg_count[seg], 1u);
}
