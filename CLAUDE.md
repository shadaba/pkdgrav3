# pkdgrav3 compressed parallel snapshot output
Read docs/compressed_io_plan.md first and follow its phases and exit criteria. Stop at each exit criterion and summarize.

## Desktop (Mac mini M4 Pro, 24 GB RAM, about 74 GB free disk)
- Keep repo + build + test data under 15 GB. Tests at 256^3 or smaller; 512^3 only as a single run, deleted afterwards.
- Repo root: /Users/shadab/Documents/Projects/pkdgrav3_dev/pkdgrav3. Build in build.nosync/, test data in tmp.nosync/ (both inside the repo, untracked). Delete outputs when done.
- macOS has no /dev/shm and no O_DIRECT. No CUDA.
- Mac build (pkdgrav3 3.5.1): `cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DMETAL=OFF -DHDF5_DIR=/opt/homebrew/opt/hdf5/lib/cmake/hdf5 -DBoost_DIR=/opt/homebrew/opt/boost/lib/cmake/Boost-1.92.0 -DPython3_EXECUTABLE=/Users/shadab/Documents/Projects/pkdgrav3_dev/venv/bin/python` (METAL=OFF: only Command Line Tools, no metal compiler; HDF5_DIR/Boost_DIR stop CMake mixing anaconda headers with brew libs; the venv has numpy+Cython).
- If pkdgrav3 does not build on macOS after a reasonable effort, say so and develop the codec standalone (no pkdgrav3 headers), then integrate and build on pawna.

## pawna (ssh alias: pawna)
- Work dir: /storage/shadab/pkdgrav3_io/ (NFS). Clone: /storage/shadab/pkdgrav3_io/pkdgrav3 (branch feature/compressed-output); builds in /storage/shadab/pkdgrav3_io/build; benchmark data/logs in /storage/shadab/pkdgrav3_io/bench. Home (/user/shadab, local ext4 on the login node): do not use until Phase 0 compares it.
- Login node: read-only probes and short builds only. Benchmarks and runs go through PBS (qsub). Show me the job script and wait for approval before submitting.
- Check node-local scratch from inside a short PBS job, not the login node.
- Never run a simulation larger than 512^3 (pawna or Mac) without my explicit approval for that run.
- Production pkdgrav3 installs and simulation directories are read-only. Never modify or delete anything outside the work dir.
- Sync: commit on the Mac, push to origin (my fork), git pull on pawna.

## Code rules
- Branch feature/compressed-output. Small commits. Never push to upstream.
- New code in new files. Keep edits to existing files minimal. The codec must not include pkdgrav3 headers.
- State how every measurement was made.
