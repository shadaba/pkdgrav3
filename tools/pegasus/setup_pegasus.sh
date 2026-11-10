#!/bin/bash -l
# One-time setup of pkdgrav3 (compressed-output branch) on pegasus from GitHub, in a new directory.
#   bash setup_pegasus.sh /data/project/hpc2502016/shadab/pkdgrav3_io
# Creates <root>/env (conda environment, ~3 GB; package cache in <root>/conda_pkgs, not in $HOME which has a
# small quota), <root>/pkdgrav3 (clone), <root>/build, <root>/build_snapcomp. Needs conda ($CONDA_BASE,
# default ~/miniconda3). Run on a login node; then submit bench/pegasus/smoke.pbs to the debug queue.
set -eu
ROOT=${1:?usage: setup_pegasus.sh <root directory>}
REPO=${REPO:-https://github.com/shadaba/pkdgrav3.git}
BRANCH=${BRANCH:-feature/compressed-output}
CONDA_BASE=${CONDA_BASE:-$HOME/miniconda3}
mkdir -p "$ROOT" && cd "$ROOT"
source "$CONDA_BASE/etc/profile.d/conda.sh"
export CONDA_PKGS_DIRS=$ROOT/conda_pkgs
if [ ! -d env ]; then
  conda create -y -q -p "$ROOT/env" --override-channels -c conda-forge \
    python=3.11 "numpy<2" "cython>=3" dill tomli h5py zstandard \
    c-compiler cxx-compiler make cmake pkg-config git \
    "openmpi=5" "fftw=*=mpi_openmpi*" "hdf5=*=nompi*" gsl "libboost-devel" libhwloc zstd
fi
conda activate "$ROOT/env"
python -c "import numpy, Cython, dill, tomli, h5py, zstandard; print('env ok: numpy', numpy.__version__, 'Cython', Cython.__version__)"
[ -e pkdgrav3 ] || git clone -b "$BRANCH" "$REPO" pkdgrav3
(cd pkdgrav3 && git submodule update --init)
mkdir -p bench/logs
PKD_IO_ROOT="$ROOT" CONDA_BASE="$CONDA_BASE" bash -l pkdgrav3/tools/pegasus/compile_pegasus.sh
echo "setup done in $ROOT. Next: qsub -v PKD_IO_ROOT=$ROOT $ROOT/pkdgrav3/bench/pegasus/smoke.pbs"
