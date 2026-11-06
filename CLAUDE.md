# pkdgrav3 compressed parallel snapshot output
Read docs/compressed_io_plan.md first and follow its phases and exit criteria. Stop at each exit criterion and summarize.

## Desktop (Mac mini M4 Pro, 24 GB RAM, about 74 GB free disk)
- Keep repo + build + test data under 15 GB. Tests at 256^3 or smaller; 512^3 only as a single run, deleted afterwards.
- Test data goes in /Users/shadab/Documents/Projects/pkdgrav3_dev/tmp. Delete outputs when done.
- macOS has no /dev/shm and no O_DIRECT. No CUDA.
- If pkdgrav3 does not build on macOS after a reasonable effort, say so and develop the codec standalone (no pkdgrav3 headers), then integrate and build on pawna.

## pawna (ssh alias: pawna)
- Work dir: /storage/shadab/pkdgrav3_io/ (source, builds, data). Home: do not use until Phase 0 compares it.
- Login node: read-only probes and short builds only. Benchmarks and runs go through PBS (qsub). Show me the job script and wait for approval before submitting.
- Check node-local scratch from inside a short PBS job, not the login node.
- Production pkdgrav3 installs and simulation directories are read-only. Never modify or delete anything outside the work dir.
- Sync: commit on the Mac, push to origin (my fork), git pull on pawna.

## Code rules
- Branch feature/compressed-output. Small commits. Never push to upstream.
- New code in new files. Keep edits to existing files minimal. The codec must not include pkdgrav3 headers.
- State how every measurement was made.
