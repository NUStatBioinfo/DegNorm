===========================================================
DegNorm TODO
===========================================================

## As of ``v0.1.4`` release
- Add "supply your own reads" support for users who have already computed their own read counts (e.g. with htseq)
- Add ability to drop genes that do not have a minimum number or fraction of samples with nonzero coverage (e.g. drop gene if > 2 / 5 samples have no reads for said gene)
- Determine why reads that do not fall between any gene_start and gene_end still persist after dot product trick with transcript vector in reads.py (requiring sqlite3 hack)