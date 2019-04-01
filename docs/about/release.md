# Release Notes

## ``v0.1.4`` (April 2019)
- DegNorm accepted to Genome Biology!
- Accounting bug fixes: correctly use pysam-loaded .bam files as 0-indexed.
- Separate read counting for genes with mutual overlap vs. isolated genes.
- In `coverage_read_counts` save parsed counts, coverage to disk instead of returning complex file tree dict.
- Don't drop genes just because they may have a sample with 0 coverage.
- Add notes about 0-indexing (for .bam files) and 1-indexing (for .gtf files) to docs.

## ``v0.1.3`` (January 2019)
- Improve upon memory efficiency in `coverage.gene_coverage` - now sparse
chromosome-wide coverage matrices are diced into approximately equally-sized chunks so that
 the dense versions of the coverage chunks are approximately ~0.5Gb memory, instead of
 calling <sparse>.todense() on a whole chromosome's coverage matrix.
- Enlist `HTSeq` package in .gtf file parsing to speed up identification of genes with
intersecting exons. GeneAnnotationProcessor now single-process, much faster.
- GH Readme points to DegNorm homepage.

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