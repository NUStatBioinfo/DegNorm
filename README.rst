===========================================================
DegNorm: Normalizing RNA degradation in RNA-Seq experiments
===========================================================

RNA-Seq transcriptome sequencing experiments often undergo gene- *and* sample-specific transcript degradation, thus
potentially biasing results in standard RNA-Seq analysis (e.g. differential expression analysis). The DegNorm pipeline
uses non-negative matrix factorization with an over-approximation constraint to normalized degraded gene transcript
coverage curves.

Users supply `.sam` files from paired-read RNA-seq experiments and a genome annotation file, and the DegNorm
pipeline will generate everything from coverage matrices, normalized coverage matrices, plots, and a report
(see **`degnorm pipeline output`**).

**Pipeline steps**


**`degnorm` usage**


**`degnorm` pipeline output**

============
Installation
============
