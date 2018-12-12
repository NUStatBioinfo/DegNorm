# Release Notes

## ``v0.1.2`` (December 2018)
- ``degnorm_mpi`` entrypoint available for distributed pipeline runs.
- change `-c` (cores flag) to `-p` (processes per node flag) to reflect that degnorm_mpi
runs processes within a node.
- Update `install.sh` so that if `mpiexec` is in `PATH`, automatically install `mpi4py`, otherwise
don't.
- Update installation documentation to reflect MPI requirements to run `degnorm_mpi`.


## ``v0.1.1`` (November 2018)
- First release of DegNorm as a package.
- ``degnorm`` entrypoint is the only supported CLI tool
- ``get_coverage_plots`` and ``get_coverage_data`` are the only supported posthoc analysis tools