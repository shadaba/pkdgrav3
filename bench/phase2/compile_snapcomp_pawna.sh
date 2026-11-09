#!/bin/bash -l
# Build the standalone snapcomp library, tests and pkc_convert on the pawna login node (short build).
# zstd headers are not installed on pawna, so CMake fetches zstd 1.5.6; HDF5 is the system 1.10 (serial).
set -eu
W=/storage/shadab/pkdgrav3_io
source $W/pkdgrav3/bench/phase0/env_pawna.sh
cmake -S $W/pkdgrav3/io/snapcomp -B $W/build_snapcomp -DCMAKE_BUILD_TYPE=Release > $W/build_snapcomp.cmake.log 2>&1
cmake --build $W/build_snapcomp -j 16 > $W/build_snapcomp.log 2>&1
ls -la $W/build_snapcomp/snapcomp_test $W/build_snapcomp/pkc_convert
