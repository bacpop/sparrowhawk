//! GPU-accelerated radix sort + count + filter for kmer preprocessing.
//!
//! This module implements the bulk-mode kmer pipeline using wgpu compute shaders:
//! - 8-pass LSD radix sort (256 buckets per pass) on canonical hash
//! - Count runs of identical canonical hashes
//! - Filter by min_count threshold
//!
//! Data is partitioned into BUCKET_COUNT=64 buckets (top 6 bits of canonical hash)
//! at extraction time so each bucket fits within GPU buffer limits.

use std::{
    cell::RefCell,
    collections::HashMap,
    hash::BuildHasherDefault,
    path::PathBuf,
};

use bytemuck::{Pod, Zeroable};
use nohash_hasher::NoHashHasher;
use wgpu::util::DeviceExt;

use crate::HashInfoSimple;

/// Number of partition buckets (top 6 bits of canonical hash → 64 buckets).
pub const BUCKET_COUNT: usize = 64;

// ──────────────────────────────────────────────────────────────────────────────
// GPU-side data structures (must match WGSL structs exactly)
// ──────────────────────────────────────────────────────────────────────────────

/// One kmer occurrence uploaded to the GPU (24 bytes, aligned).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuKmerEntry {
    /// Canonical hash bits [31:0]
    pub canon_lo: u32,
    /// Canonical hash bits [63:32]
    pub canon_hi: u32,
    /// Non-canonical hash bits [31:0]
    pub nc_lo:    u32,
    /// Non-canonical hash bits [63:32]
    pub nc_hi:    u32,
    /// First/last bases (u8 padded to u32)
    pub bases:    u32,
    /// Padding to 24 bytes
    pub pad:      u32,
}

/// One filtered kmer result returned from the GPU (24 bytes, aligned).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GpuKmerOut {
    /// Canonical hash bits [31:0]
    pub canon_lo: u32,
    /// Canonical hash bits [63:32]
    pub canon_hi: u32,
    /// Non-canonical hash bits [31:0]
    pub nc_lo:    u32,
    /// Non-canonical hash bits [63:32]
    pub nc_hi:    u32,
    /// First/last bases
    pub bases:    u32,
    /// Occurrence count
    pub count:    u32,
}

// ──────────────────────────────────────────────────────────────────────────────
// Pipeline cache
// ──────────────────────────────────────────────────────────────────────────────

struct GpuPipelines {
    radix_count:      wgpu::ComputePipeline,
    radix_prefix_g:   wgpu::ComputePipeline, // global histogram exclusive scan
    radix_wg_offset:  wgpu::ComputePipeline, // column-wise wg prefix scan
    radix_scatter:    wgpu::ComputePipeline,
    prefix_local:    wgpu::ComputePipeline,
    prefix_block:    wgpu::ComputePipeline,
    prefix_prop:     wgpu::ComputePipeline,
    boundary:        wgpu::ComputePipeline,
    count_seg:       wgpu::ComputePipeline,
    filter_mark:     wgpu::ComputePipeline,
    scatter_out:     wgpu::ComputePipeline,
}

fn compile_pipelines(device: &wgpu::Device) -> GpuPipelines {
    let mk = |src: &str| device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: None,
        source: wgpu::ShaderSource::Wgsl(src.into()),
    });

    let sm_radix_count      = mk(include_str!("shaders/radix_count.wgsl"));
    let sm_radix_prefix_g   = mk(include_str!("shaders/radix_prefix_global.wgsl"));
    let sm_radix_wg_offset  = mk(include_str!("shaders/radix_wg_offset.wgsl"));
    let sm_radix_scatter    = mk(include_str!("shaders/radix_scatter.wgsl"));
    let sm_prefix_sum    = mk(include_str!("shaders/prefix_sum.wgsl"));
    let sm_boundary      = mk(include_str!("shaders/boundary_detect.wgsl"));
    let sm_count_seg     = mk(include_str!("shaders/count_segments.wgsl"));
    let sm_filter_mark   = mk(include_str!("shaders/filter_mark.wgsl"));
    let sm_scatter_out   = mk(include_str!("shaders/scatter_output.wgsl"));

    let cp = |sm: &wgpu::ShaderModule, ep: &str| {
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(ep),
            layout: None,
            module: sm,
            entry_point: Some(ep),
            compilation_options: Default::default(),
            cache: Default::default(),
        })
    };

    GpuPipelines {
        radix_count:      cp(&sm_radix_count,     "main"),
        radix_prefix_g:   cp(&sm_radix_prefix_g,  "main"),
        radix_wg_offset:  cp(&sm_radix_wg_offset, "main"),
        radix_scatter:    cp(&sm_radix_scatter,    "main"),
        prefix_local:    cp(&sm_prefix_sum,    "local_scan"),
        prefix_block:    cp(&sm_prefix_sum,    "block_scan"),
        prefix_prop:     cp(&sm_prefix_sum,    "propagate"),
        boundary:        cp(&sm_boundary,      "main"),
        count_seg:       cp(&sm_count_seg,     "main"),
        filter_mark:     cp(&sm_filter_mark,   "main"),
        scatter_out:     cp(&sm_scatter_out,   "main"),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Helper: run 3-pass prefix sum on an arbitrary u32 buffer
// ──────────────────────────────────────────────────────────────────────────────

fn run_prefix_sum(
    device:  &wgpu::Device,
    queue:   &wgpu::Queue,
    pipes:   &GpuPipelines,
    buf:     &wgpu::Buffer,  // read_write storage, length n
    n:       u32,
) {
    let num_wg = n.div_ceil(256) as u32;

    // block_sums buffer: one u32 per workgroup
    let block_sums_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("block_sums"),
        size:               (num_wg as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    // ---- Pass A: local_scan ----
    {
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  None,
            layout: &pipes.prefix_local.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: block_sums_buf.as_entire_binding() },
            ],
        });
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipes.prefix_local);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(num_wg, 1, 1);
        }
        queue.submit([enc.finish()]);
    }

    if num_wg > 1 {
        // ---- Pass B: block_scan (single workgroup over block_sums) ----
        // block_scan only accesses block_sums (binding 1); the reflected layout has 1 binding.
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.prefix_block.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 1, resource: block_sums_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                pass.set_pipeline(&pipes.prefix_block);
                pass.set_bind_group(0, &bg, &[]);
                pass.dispatch_workgroups(1, 1, 1);
            }
            queue.submit([enc.finish()]);
        }

        // ---- Pass C: propagate ----
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.prefix_prop.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: block_sums_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                pass.set_pipeline(&pipes.prefix_prop);
                pass.set_bind_group(0, &bg, &[]);
                pass.dispatch_workgroups(num_wg, 1, 1);
            }
            queue.submit([enc.finish()]);
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Per-bucket GPU sort + count + filter
// ──────────────────────────────────────────────────────────────────────────────

/// Convert CPU tuple vec to GPU entry vec.
fn to_gpu_entries(bucket: &[(u64, u64, u8)]) -> Vec<GpuKmerEntry> {
    bucket.iter().map(|&(hc, hnc, b)| GpuKmerEntry {
        canon_lo: hc as u32,
        canon_hi: (hc >> 32) as u32,
        nc_lo:    hnc as u32,
        nc_hi:    (hnc >> 32) as u32,
        bases:    b as u32,
        pad:      0,
    }).collect()
}

/// Run GPU radix sort on the elements buffer (ping-pong A↔B over 8 passes).
/// Returns the buffer that holds the sorted result.
fn gpu_radix_sort(
    device:  &wgpu::Device,
    queue:   &wgpu::Queue,
    pipes:   &GpuPipelines,
    entries: &[GpuKmerEntry],
) -> wgpu::Buffer {
    let n    = entries.len() as u32;
    let size = (entries.len() * std::mem::size_of::<GpuKmerEntry>()) as u64;
    let num_wg = n.div_ceil(256) as u32;

    // Two element buffers for ping-pong
    let buf_a = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label:    Some("sort_a"),
        contents: bytemuck::cast_slice(entries),
        usage:    wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
    });
    let buf_b = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("sort_b"),
        size,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // Histogram buffers: global (256 u32) and per-wg (num_wg * 256 u32)
    let global_hist_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("global_hist"),
        size:               256 * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let wg_hist_size = (num_wg as u64) * 256 * 4;
    let wg_hist_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("wg_hist"),
        size:               wg_hist_size,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let global_prefix_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("global_prefix"),
        size:               256 * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let wg_prefix_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("wg_prefix"),
        size:               wg_hist_size,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let zeros_256: Vec<u8> = vec![0u8; 256 * 4];
    let zeros_wg:  Vec<u8> = vec![0u8; wg_hist_size as usize];

    // pass_idx uniform (one u32 per pass → reuse same buffer, write new value each pass)
    let pass_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label:    Some("pass_idx"),
        contents: bytemuck::bytes_of(&0u32),
        usage:    wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // num_wg uniform — used by radix_wg_offset to know how many workgroups to scan
    let num_wg_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label:    Some("num_wg"),
        contents: bytemuck::bytes_of(&num_wg),
        usage:    wgpu::BufferUsages::UNIFORM,
    });

    let bufs = [&buf_a, &buf_b];

    for pass in 0u32..8u32 {
        let src_buf = bufs[(pass % 2) as usize];
        let dst_buf = bufs[((pass + 1) % 2) as usize];

        // Update pass index uniform
        queue.write_buffer(&pass_buf, 0, bytemuck::bytes_of(&pass));

        // Clear histograms
        queue.write_buffer(&global_hist_buf,   0, &zeros_256);
        queue.write_buffer(&wg_hist_buf,       0, &zeros_wg);
        queue.write_buffer(&global_prefix_buf, 0, &zeros_256);
        queue.write_buffer(&wg_prefix_buf,     0, &zeros_wg);

        // ---- Dispatch 1: radix_count ----
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.radix_count.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: src_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: global_hist_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: wg_hist_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 3, resource: pass_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass_enc = enc.begin_compute_pass(&Default::default());
                pass_enc.set_pipeline(&pipes.radix_count);
                pass_enc.set_bind_group(0, &bg, &[]);
                pass_enc.dispatch_workgroups(num_wg, 1, 1);
            }
            queue.submit([enc.finish()]);
        }

        // ---- Dispatch 2a: radix_prefix_global (global_hist → global_prefix) ----
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.radix_prefix_g.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: global_hist_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: global_prefix_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass_enc = enc.begin_compute_pass(&Default::default());
                pass_enc.set_pipeline(&pipes.radix_prefix_g);
                pass_enc.set_bind_group(0, &bg, &[]);
                pass_enc.dispatch_workgroups(1, 1, 1);
            }
            queue.submit([enc.finish()]);
        }

        // ---- Dispatch 2b: radix_wg_offset (column-wise scan: wg_hist → wg_prefix) ----
        // Thread d scans wg_hist[*][d] across all num_wg workgroups → exclusive column prefix.
        // Dispatched as a single workgroup of 256 threads (one thread per digit value).
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.radix_wg_offset.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: wg_hist_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: wg_prefix_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: num_wg_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass_enc = enc.begin_compute_pass(&Default::default());
                pass_enc.set_pipeline(&pipes.radix_wg_offset);
                pass_enc.set_bind_group(0, &bg, &[]);
                pass_enc.dispatch_workgroups(1, 1, 1);
            }
            queue.submit([enc.finish()]);
        }

        // ---- Dispatch 3: radix_scatter ----
        {
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label:  None,
                layout: &pipes.radix_scatter.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: src_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: dst_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: global_prefix_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 3, resource: wg_prefix_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 4, resource: pass_buf.as_entire_binding() },
                ],
            });
            let mut enc = device.create_command_encoder(&Default::default());
            {
                let mut pass_enc = enc.begin_compute_pass(&Default::default());
                pass_enc.set_pipeline(&pipes.radix_scatter);
                pass_enc.set_bind_group(0, &bg, &[]);
                pass_enc.dispatch_workgroups(num_wg, 1, 1);
            }
            queue.submit([enc.finish()]);
        }
    }

    // After 8 passes (even), result is in buf_a
    buf_a
}

/// Run GPU count+filter pipeline. Returns output buffer and output count.
fn gpu_count_filter(
    device:    &wgpu::Device,
    queue:     &wgpu::Queue,
    pipes:     &GpuPipelines,
    sorted:    &wgpu::Buffer,
    n:         u32,
    min_count: u32,
) -> (wgpu::Buffer, u32) {
    let num_wg = n.div_ceil(256) as u32;

    // boundary[n]
    let boundary_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("boundary"),
        size:               (n as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // ---- Dispatch 1: boundary_detect ----
    {
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  None,
            layout: &pipes.boundary.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: sorted.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: boundary_buf.as_entire_binding() },
            ],
        });
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipes.boundary);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(num_wg, 1, 1);
        }
        queue.submit([enc.finish()]);
    }

    // segment_id = prefix_sum(boundary) → then subtract 1 per element on GPU
    // We clone boundary into segment_id_buf, run prefix_sum, subtract 1
    let segment_id_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("segment_id"),
        size:               (n as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    // Copy boundary → segment_id_buf
    {
        let mut enc = device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(&boundary_buf, 0, &segment_id_buf, 0, (n as u64) * 4);
        queue.submit([enc.finish()]);
    }

    // Run prefix sum on segment_id_buf to get inclusive prefix (each element = its 1-based segment id).
    // But we need exclusive prefix for 0-based segment id, which prefix_sum.wgsl produces.
    // After the exclusive prefix sum, segment_id[i] = number of boundaries before i = 0-based segment index. ✓
    run_prefix_sum(device, queue, pipes, &segment_id_buf, n);

    // Read back number of unique kmers from the last element of segment_id + boundary
    // (easier: just use n_unique = last(segment_id) + last(boundary))
    // We'll read it back after we have the result count.

    // We need to know num_unique_kmers to size seg_count. Conservative upper bound = n.
    // In practice we could readback but that stalls pipeline; use n as upper bound.
    let num_unique = n; // upper bound; seg_count[n_unique..n] stays zero, harmless

    // seg_count[num_unique] (atomic)
    let seg_count_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("seg_count"),
        size:               (num_unique as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&seg_count_buf, 0, &vec![0u8; (num_unique as usize) * 4]);

    // ---- Dispatch 5: count_segments ----
    {
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  None,
            layout: &pipes.count_seg.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: segment_id_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: seg_count_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: boundary_buf.as_entire_binding() },
            ],
        });
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipes.count_seg);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(num_wg, 1, 1);
        }
        queue.submit([enc.finish()]);
    }

    // min_count uniform
    let min_count_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label:    Some("min_count"),
        contents: bytemuck::bytes_of(&min_count),
        usage:    wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // output_mask[n]
    let output_mask_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("output_mask"),
        size:               (n as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // ---- Dispatch 6: filter_mark ----
    {
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  None,
            layout: &pipes.filter_mark.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: boundary_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: segment_id_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: seg_count_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: output_mask_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: min_count_buf.as_entire_binding() },
            ],
        });
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipes.filter_mark);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(num_wg, 1, 1);
        }
        queue.submit([enc.finish()]);
    }

    // output_index = prefix_sum(output_mask) → exclusive, so output_index[i] = write position
    let output_index_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("output_index"),
        size:               (n as u64) * 4,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    {
        let mut enc = device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(&output_mask_buf, 0, &output_index_buf, 0, (n as u64) * 4);
        queue.submit([enc.finish()]);
    }
    run_prefix_sum(device, queue, pipes, &output_index_buf, n);

    // Read back output count: last(output_mask) + last(output_index)
    // We read them back to size the output buffer exactly.
    let readback_size = 2u64 * 4;
    let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("readback_count"),
        size:               readback_size,
        usage:              wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    {
        let mut enc = device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(&output_index_buf, ((n - 1) as u64) * 4, &readback_buf, 0, 4);
        enc.copy_buffer_to_buffer(&output_mask_buf,  ((n - 1) as u64) * 4, &readback_buf, 4, 4);
        queue.submit([enc.finish()]);
    }
    // Poll + map
    let output_count = {
        let slice = readback_buf.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| { tx.send(r).unwrap(); });
        let _ = device.poll(wgpu::PollType::wait_indefinitely());
        rx.recv().unwrap().unwrap();
        let data = slice.get_mapped_range();
        let vals: &[u32] = bytemuck::cast_slice(&data);
        vals[0] + vals[1] // exclusive_prefix[last] + mask[last]
    };
    readback_buf.unmap();

    if output_count == 0 {
        // Return empty buffer
        let empty = device.create_buffer(&wgpu::BufferDescriptor {
            label:              Some("output_empty"),
            size:               4,
            usage:              wgpu::BufferUsages::STORAGE,
            mapped_at_creation: false,
        });
        return (empty, 0);
    }

    // Output buffer
    let output_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("kmer_out"),
        size:               (output_count as u64) * std::mem::size_of::<GpuKmerOut>() as u64,
        usage:              wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // ---- Dispatch 10: scatter_output ----
    {
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  None,
            layout: &pipes.scatter_out.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: sorted.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: output_mask_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: output_index_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: segment_id_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: seg_count_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: output_buf.as_entire_binding() },
            ],
        });
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipes.scatter_out);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(num_wg, 1, 1);
        }
        queue.submit([enc.finish()]);
    }

    (output_buf, output_count)
}

/// Readback GPU output buffer to CPU Vec<GpuKmerOut>.
fn readback_output(
    device: &wgpu::Device,
    queue:  &wgpu::Queue,
    buf:    &wgpu::Buffer,
    count:  u32,
) -> Vec<GpuKmerOut> {
    if count == 0 { return Vec::new(); }
    let size = count as u64 * std::mem::size_of::<GpuKmerOut>() as u64;
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label:              Some("staging"),
        size,
        usage:              wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut enc = device.create_command_encoder(&Default::default());
    enc.copy_buffer_to_buffer(buf, 0, &staging, 0, size);
    queue.submit([enc.finish()]);

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| { tx.send(r).unwrap(); });
    let _ = device.poll(wgpu::PollType::wait_indefinitely());
    rx.recv().unwrap().unwrap();

    let data = slice.get_mapped_range();
    let result: Vec<GpuKmerOut> = bytemuck::cast_slice(&data).to_vec();
    drop(data);
    staging.unmap();
    result
}

/// Process one bucket: GPU radix sort → count+filter → readback.
pub fn process_bucket(
    device:    &wgpu::Device,
    queue:     &wgpu::Queue,
    pipes:     &GpuPipelines,
    bucket:    &[(u64, u64, u8)],
    min_count: u16,
) -> Vec<GpuKmerOut> {
    if bucket.is_empty() { return Vec::new(); }

    let entries = to_gpu_entries(bucket);
    let n = entries.len() as u32;

    // Sort
    let sorted_buf = gpu_radix_sort(device, queue, pipes, &entries);

    // Count + filter
    let (out_buf, out_count) = gpu_count_filter(
        device, queue, pipes, &sorted_buf, n, min_count as u32,
    );

    // Readback
    readback_output(device, queue, &out_buf, out_count)
}

// ──────────────────────────────────────────────────────────────────────────────
// Public API
// ──────────────────────────────────────────────────────────────────────────────

/// GPU-accelerated radix sort + count + filter over partitioned kmer buckets.
///
/// Returns `themap` accumulating all kmers with count ≥ min_count.
pub async fn gpu_sort_count_filter(
    buckets:   &[Vec<(u64, u64, u8)>; BUCKET_COUNT],
    min_count: u16,
    _do_fit:   bool,
    _out_path: &mut Option<PathBuf>,
) -> HashMap<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>> {
    // Init GPU
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
        backends: wgpu::Backends::all(),
        ..Default::default()
    });
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        })
        .await
        .expect("No GPU adapter found");

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor {
            label:    Some("sparrowhawk_gpu"),
            ..Default::default()
        })
        .await
        .expect("Failed to get GPU device");

    log::info!("GPU: {}", adapter.get_info().name);

    // Compile pipelines once
    let pipes = compile_pipelines(&device);

    let mut themap: HashMap<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>> =
        HashMap::with_hasher(BuildHasherDefault::default());

    for (bucket_idx, bucket) in buckets.iter().enumerate() {
        if bucket.is_empty() { continue; }
        log::info!("GPU processing bucket {}/{} ({} entries)", bucket_idx + 1, BUCKET_COUNT, bucket.len());

        let results = process_bucket(&device, &queue, &pipes, bucket, min_count);

        for r in results {
            let hc  = (r.canon_lo as u64) | ((r.canon_hi as u64) << 32);
            let hnc = (r.nc_lo    as u64) | ((r.nc_hi    as u64) << 32);
            let b   = r.bases as u8;
            let cnt = r.count.min(u16::MAX as u32) as u16;

            themap.entry(hc).or_insert_with(|| RefCell::new(HashInfoSimple {
                hnc,
                b,
                pre:    Vec::new(),
                post:   Vec::new(),
                counts: cnt,
            }));
        }
    }

    themap
}
