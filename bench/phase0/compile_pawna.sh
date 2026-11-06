#!/bin/bash -l
# Build the fork (pkdgrav3 3.5.1 + our branch) on the pawna login node (short build).
# Deps: openmpi-5.0.1 module (same as production), system HDF5 1.10.7 (threadsafe), GSL, FFTW,
# Boost 1.74, and our own venv (numpy + Cython >= 3) in the work dir.
# Production links conda's HDF5 2.1.0 instead; the HDF5 version is a known difference.
set -eu
W=/storage/shadab/pkdgrav3_io
module purge
module load openmpi-5.0.1
cd $W/pkdgrav3
git pull --ff-only && git submodule update --init
cmake -S . -B $W/build -DCMAKE_BUILD_TYPE=Release \
      -DPython3_EXECUTABLE=$W/venv/bin/python > $W/build.cmake.log 2>&1
cmake --build $W/build -j 16 > $W/build.log 2>&1
ls -la $W/build/pkdgrav3
