#!/bin/bash -l
# Build pkdgrav3 (branch feature/compressed-output) and the snapcomp tools on the pegasus login node.
# PKD_IO_ROOT holds pkdgrav3/ (clone) and env/ (conda, made by setup_pegasus.sh); builds go to build/ and
# build_snapcomp/. CLEAN=1 rebuilds pkdgrav3 from scratch.
# ARCH: CPU target (-march). Default skylake-avx512 runs on every compute node type (Skylake hpc001-060,
# Ice Lake hpc061-080, Sapphire Rapids hpc081-100 incl. the debug queue). ARCH=auto = -march=native of the
# login node (Ice Lake on ln3/ln4: would not run on the Skylake nodes).
set -eu
source "$(dirname "$(readlink -f "$0")")/env_pegasus.sh"
W=$PKD_IO_ROOT
ARCH=${ARCH:-skylake-avx512}
NJ=${NJ:-8}
cd $W/pkdgrav3
git submodule update --init
cmake -S . -B $W/build -DCMAKE_BUILD_TYPE=Release -DTARGET_ARCHITECTURE=$ARCH \
      -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python -DCMAKE_PREFIX_PATH=$CONDA_PREFIX > $W/build.cmake.log 2>&1
if [ "${CLEAN:-0}" = 1 ]; then cmake --build $W/build --target clean; fi
cmake --build $W/build -j $NJ > $W/build.log 2>&1
cmake -S $W/pkdgrav3/io/snapcomp -B $W/build_snapcomp -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
      > $W/build_snapcomp.cmake.log 2>&1
cmake --build $W/build_snapcomp -j $NJ > $W/build_snapcomp.log 2>&1
ls -la $W/build/pkdgrav3 $W/build_snapcomp/snapcomp_test $W/build_snapcomp/pkc_check $W/build_snapcomp/pkc_to_gadget
