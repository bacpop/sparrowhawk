// Stable radix scatter kernel.
//
// Each workgroup builds a 256-bit bitmask per digit value in workgroup memory,
// recording which thread IDs have each digit.  Each thread then computes its
// stable local rank by counting set bits before its own thread ID in the mask
// for its digit.  This guarantees that within a workgroup, elements are
// scattered in thread-ID (= original input) order, making the scatter stable.
//
// Output position:
//   out_pos = global_prefix[digit]           (start for this digit globally)
//           + wg_prefix[wg_id * 256 + digit] (offset from all prior workgroups)
//           + local_rank                     (position within this workgroup)

struct KmerEntry {
    canon_lo: u32,
    canon_hi: u32,
    nc_lo:    u32,
    nc_hi:    u32,
    bases:    u32,
    pad:      u32,
}

@group(0) @binding(0) var<storage, read>       src:           array<KmerEntry>;
@group(0) @binding(1) var<storage, read_write> dst:           array<KmerEntry>;
@group(0) @binding(2) var<storage, read>       global_prefix: array<u32>;
@group(0) @binding(3) var<storage, read>       wg_prefix:     array<u32>;
@group(0) @binding(4) var<uniform>             pass_idx:      u32;

fn digit(e: KmerEntry, p: u32) -> u32 {
    if p < 4u {
        return (e.canon_lo >> (p * 8u)) & 0xFFu;
    } else {
        return (e.canon_hi >> ((p - 4u) * 8u)) & 0xFFu;
    }
}

// Bitmask: 256 threads × 8 u32 words (= 256 bits) per digit value, 256 digit values.
// Layout: digit_mask[d * 8 + word] covers threads [word*32 .. word*32+31] for digit d.
// Total: 256 * 8 = 2048 u32 = 8 KB workgroup memory.
var<workgroup> digit_mask: array<atomic<u32>, 2048>;

@compute @workgroup_size(256)
fn main(
    @builtin(global_invocation_id) gid:  vec3<u32>,
    @builtin(local_invocation_id)  lid:  vec3<u32>,
    @builtin(workgroup_id)         wgid: vec3<u32>,
) {
    // Clear the 8 words this thread is responsible for (lid.x * 8 .. lid.x * 8 + 7).
    // Because there are 2048 words and 256 threads, each thread clears exactly 8 words.
    for (var i = 0u; i < 8u; i++) {
        atomicStore(&digit_mask[lid.x * 8u + i], 0u);
    }
    workgroupBarrier();

    let n = arrayLength(&src);
    let valid = gid.x < n;
    var d = 0u;
    if valid {
        d = digit(src[gid.x], pass_idx);
    }

    // Set the bit for this thread in the bitmask for its digit.
    if valid {
        atomicOr(&digit_mask[d * 8u + lid.x / 32u], 1u << (lid.x % 32u));
    }
    workgroupBarrier();

    if valid {
        // Count set bits before thread lid.x in digit_mask[d]:
        //   - sum popcount of full words [0 .. my_word)
        //   - plus popcount of the partial word at my_word, masked to bits < my_bit
        let my_word = lid.x / 32u;
        let my_bit  = lid.x % 32u;
        var local_rank = 0u;
        for (var w = 0u; w < my_word; w++) {
            local_rank += countOneBits(atomicLoad(&digit_mask[d * 8u + w]));
        }
        local_rank += countOneBits(
            atomicLoad(&digit_mask[d * 8u + my_word]) & ((1u << my_bit) - 1u)
        );

        let out_pos = global_prefix[d]
                    + wg_prefix[wgid.x * 256u + d]
                    + local_rank;
        dst[out_pos] = src[gid.x];
    }
}
