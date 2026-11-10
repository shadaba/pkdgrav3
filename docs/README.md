# Documentation map

This fork adds compressed parallel snapshot output and related tools to pkdgrav3 3.5.1
(https://github.com/shadaba/pkdgrav3, branch `main`). The Markdown files in `docs/` describe that work.
The `.rst` files are upstream pkdgrav3's Sphinx manual; Sphinx does not read the `.md` files.

## Where to start

| You want to | Read |
|---|---|
| Know what was added and run a simulation with it | [compressed_output_guide.md](compressed_output_guide.md) |
| Set up and run on pawna (TIFR) | [running_on_pawna.md](running_on_pawna.md) |
| Set up and run on pegasus (IUCAA) | [running_on_pegasus.md](running_on_pegasus.md) |
| Read `.pkc` files, or write a reader | [snapcomp_format.md](snapcomp_format.md), then [../tools/pkcreader/README.md](../tools/pkcreader/README.md) |
| Find how something was measured or decided | [compressed_io_plan.md](compressed_io_plan.md) |

## The Markdown files

| File | Contents | Use it for |
|---|---|---|
| [compressed_output_guide.md](compressed_output_guide.md) | **User guide.** Feature list, measured results, every parameter, example parameter files (defaults, coarser steps, potential, potential only at low z, validation), laptop and pawna runs, staged checkpoints, stopping and restarting, a 2048^3 production example, reading and converting snapshots, recovery, tests. | Running simulations. Start here. |
| [snapcomp_format.md](snapcomp_format.md) | **File format specification**, versions 1 and 2: file naming, header fields with units, particle quantities with units and accuracy, the info JSON (simulation, box, ICs, cosmology, numerics, steps) with units, block and stream layout, decoding, subsample hash. | Writing or checking a reader; looking up a unit or header key. |
| [running_on_pawna.md](running_on_pawna.md) | **pawna operations:** one-command setup from GitHub, the build and its pitfalls, I/O parameters, `run_staged.pbs` options (local-SSD staging and the drain to NFS, checkpoints, restart, production binary), the drain tool, recovering data left on a node. | Day-to-day work on pawna. |
| [running_on_pegasus.md](running_on_pegasus.md) | **pegasus operations:** cluster facts, conda-based setup, launching on one or several nodes, test results (smoke, 2 nodes, 512^3 I/O on 4 nodes: write straight to `/data`, no staging needed). | Day-to-day work on pegasus. |
| [compressed_io_plan.md](compressed_io_plan.md) | **Project log and plan** (sections 1-20): goal, design, phases and exit criteria, results of every benchmark with how it was measured, the potential, the Ewald bug, checkpoint and restart tests, run-time and memory estimates (2048^3 to 3250^3). | Why a choice was made; the evidence behind a number. Sections are dated. |
| [memory_persistent_ids.md](memory_persistent_ids.md) | Persistent particle IDs without the 8-byte ID field (ordered mode, `iOrder`), with memory measurements. | Choosing memory settings for large runs. |
| [rockstar_pkc_session_prompt.md](rockstar_pkc_session_prompt.md) | Hand-off prompt for a separate session that adds native `.pkc` reading to Rockstar in sahyadri-sandbox. | Starting that work. |
| [upstream/ewald_bug_report.md](upstream/ewald_bug_report.md) and [upstream/ewald_simd_remainder_sign.patch](upstream/ewald_simd_remainder_sign.patch) | Report and one-line patch for the pkdgrav3 Ewald potential bug, to send to the pkdgrav3 developers. | Reporting the bug upstream. |

## Elsewhere in the repository

| Path | Contents |
|---|---|
| `io/snapcomp/` | Codec library (C++, no pkdgrav3 headers) and tools: `pkc_check` (IDs, subsamples, persistence), `pkc_to_gadget` (Gadget-2 for Rockstar), `pkc_convert` (HDF5 to `.pkc`), `snapcomp_test` (unit tests). |
| `core/writecompressed.*`, `msr_writecompressed.cxx` | The writer inside pkdgrav3 (per-thread compression, header and info JSON). |
| `parameters.toml` | All parameters with help text. The compressed-output ones are `bCompressedOutput`, `dComp*`, `lstComp*`, `iComp*`, `bCompAlsoStandard` and `achSim*`. |
| `tools/pkcreader/` | Python reader (`pkcreader.py`, numpy + zstandard) with its README and test against HDF5. |
| `tools/staging/` | `run_staged.pbs` (job wrapper), `pkd_drain.py` (copies snapshots and checkpoints off local disk), `test_drain.py`. |
| `tools/pawna/`, `tools/pegasus/` | Setup, build and environment scripts per cluster. |
| `tools/snapcomp_proto/` | Phase-1 Python prototype of the codec (kept for reference; see its README). |
| `bench/phase0` to `bench/phase5` | PBS jobs behind the numbers in compressed_io_plan.md, by project phase. |
| `bench/pegasus/` | pegasus tests: smoke, 2-node MPI, 512^3 I/O, large-N memory and timing. |

## Conventions
- Every measurement in these files says how it was made (machine, job, size, method).
- Section numbers in compressed_io_plan.md are stable; other files refer to them.
- Paths such as `/storage/shadab/...` and `/data/project/...` are the authors' work directories. Replace
  them with your own (`PKD_IO_ROOT` in the scripts).
