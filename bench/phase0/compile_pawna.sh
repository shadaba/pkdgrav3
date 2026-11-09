#!/bin/bash -l
# Build pkdgrav3 (3.5.1 + branch feature/compressed-output), the snapcomp tools and the HDF5 gzip shim on
# the pawna login node (short builds). PKD_IO_ROOT (default /storage/shadab/pkdgrav3_io) holds pkdgrav3/
# (the clone) and venv/; builds go to build/ and build_snapcomp/. CLEAN=1 rebuilds pkdgrav3 from scratch.
# Deps: openmpi-5.0.1 module (as production), system HDF5 1.10.7, GSL, FFTW, Boost 1.74, the venv.
set -eu
source "$(dirname "$(readlink -f "$0")")/env_pawna.sh"
W=$PKD_IO_ROOT
cd $W/pkdgrav3
git submodule update --init
cmake -S . -B $W/build -DCMAKE_BUILD_TYPE=Release -DPython3_EXECUTABLE=$W/venv/bin/python > $W/build.cmake.log 2>&1
if [ "${CLEAN:-0}" = 1 ]; then cmake --build $W/build --target clean; fi
cmake --build $W/build -j 16 > $W/build.log 2>&1
# standalone snapcomp tools: tests, pkc_check, pkc_convert (zstd is fetched; pawna has no zstd headers)
cmake -S $W/pkdgrav3/io/snapcomp -B $W/build_snapcomp -DCMAKE_BUILD_TYPE=Release > $W/build_snapcomp.cmake.log 2>&1
cmake --build $W/build_snapcomp -j 16 > $W/build_snapcomp.log 2>&1
make -s -C $W/pkdgrav3/tools/h5shim
ls -la $W/build/pkdgrav3 $W/build_snapcomp/snapcomp_test $W/build_snapcomp/pkc_check $W/build_snapcomp/pkc_convert $W/pkdgrav3/tools/h5shim/libh5_deflate_shim.so
