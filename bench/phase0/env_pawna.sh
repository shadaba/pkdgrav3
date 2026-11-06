# Runtime environment for our pawna pkdgrav3 build. Source from job scripts (bash -l).
module purge
module load openmpi-5.0.1
# login profile adds Intel oneAPI MPI paths; keep them out of compile and run time
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
# pkdgrav3 embeds system python3.10; give it the venv's packages (dill, tomli, Cython)
export PYTHONPATH=/storage/shadab/pkdgrav3_io/venv/lib/python3.10/site-packages${PYTHONPATH:+:$PYTHONPATH}
