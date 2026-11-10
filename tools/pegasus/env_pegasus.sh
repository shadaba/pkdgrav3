# Runtime and build environment for pkdgrav3 on pegasus (IUCAA). Source from job scripts (bash -l).
# Everything comes from one conda environment in the work directory ($PKD_IO_ROOT/env): compilers,
# Open MPI, HDF5, FFTW (MPI + threads), GSL, Boost, hwloc, CMake, Python with numpy/Cython/dill.
# Pegasus has no system MPI, CMake, HDF5 or FFTW; the pawna scripts (bench/phase0/env_pawna.sh) are separate.
# PKD_IO_ROOT: directory holding pkdgrav3/ (clone), env/, build/, build_snapcomp/.
export PKD_IO_ROOT=${PKD_IO_ROOT:-/data/project/hpc2502016/shadab/pkdgrav3_io}
CONDA_BASE=${CONDA_BASE:-$HOME/miniconda3}
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$PKD_IO_ROOT/env"
# keep the cluster's Intel oneAPI paths (if any profile adds them) out of compile and run time
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
# pkdgrav3 embeds the env's Python; make sure it finds the env's packages and libraries
export PYTHONHOME=$CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
# conda's UCX also probes CUDA; on CPU nodes that prints "cuDeviceGetCount ... initialization error". Skip it.
export UCX_TLS=${UCX_TLS:-^cuda}
# mpirun on several nodes (ssh launch): pass these to the remote processes with mpirun_pegasus_args
mpirun_pegasus_args="--prefix $CONDA_PREFIX -x PATH -x LD_LIBRARY_PATH -x PYTHONHOME -x UCX_TLS"
