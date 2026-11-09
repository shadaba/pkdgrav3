#!/bin/bash -l
# One-time setup of pkdgrav3 (compressed-output branch) on pawna from GitHub, in a new directory.
#   bash setup_pawna.sh /storage/$USER/pkdgrav3_io        (directory must not contain pkdgrav3/ yet)
# Creates <root>/pkdgrav3 (clone), <root>/venv (Python packages), <root>/build, <root>/build_snapcomp.
# Run on the login node; takes about 20-40 min on NFS. Then submit the smoke test (see docs).
set -eu
ROOT=${1:?usage: setup_pawna.sh <root directory>}
REPO=${REPO:-https://github.com/shadaba/pkdgrav3.git}
BRANCH=${BRANCH:-feature/compressed-output}
mkdir -p "$ROOT" && cd "$ROOT"
[ -e pkdgrav3 ] && { echo "$ROOT/pkdgrav3 exists; use another directory or update it with git pull"; exit 1; }
git clone -b "$BRANCH" "$REPO" pkdgrav3
(cd pkdgrav3 && git submodule update --init)
# venv on the system Python 3.10, with the system NumPy 1.21 (its headers are what the compiler finds)
python3 -m venv --system-site-packages venv
venv/bin/python -m pip install -q "Cython>=3.0" "dill>=0.3.6" tomli "h5py<3.12" zstandard
venv/bin/python -c "import numpy, Cython, dill, tomli, h5py, zstandard; print('venv ok: numpy', numpy.__version__, 'Cython', Cython.__version__)"
mkdir -p bench/logs
PKD_IO_ROOT="$ROOT" bash -l pkdgrav3/bench/phase0/compile_pawna.sh
echo "setup done in $ROOT. Next: qsub -v PKD_IO_ROOT=$ROOT $ROOT/pkdgrav3/bench/phase0/smoke32.pbs"
