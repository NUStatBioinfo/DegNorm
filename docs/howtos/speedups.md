# DegNorm speed enhancements

## Downsampling coverage matrices
The backbone of the DegNorm normalization pipeline is the non-negative matrix factorization with over-approximation (NMFOA) algorithm.
 NMFOA is a computationally expensive algorithm, requiring `--nmf-iter`-many singular value decompositions of the coverage matrix.
 The SVD is O(n^3) in terms of flop count (according to Trefethen and Bau), so genes with smaller coverage matrices will 
 naturally run through NMFOA faster than larger ones.
 
In an effort to trim down the coverage matrices, the `--downsample-rate` argument in `degnorm` allows you to use a [systematic sample](https://en.wikipedia.org/wiki/Systematic_sampling) 
for each gene's coverage matrix in lieu of the larger original ones. The `--downsample-rate` is just the integer-valued
"take-every" step size, so that the original coverage matrix is reduced by a factor of 1 / (downsample rate). At a modest
 downsampling rate, the overall shape of the coverage curves should still be recognizable.
 
![systematic_sample](img/systematic_sample.png)

Obviously, the higher the downsampling rate, the faster the DegNorm iterations will complete, because DegNorm
executes SVD approximations on smaller and smaller gene coverage matrices. The figure below illustrates the time savings
at various downsampling rates from the 9-sample gliablastoma (GBM) cell line study data set:

![gbm_dsample_speed_times](img/gbm_dsample_speed_times.png)

We recreated the DE analysis using the GBM data as well as
the Sequencing Quality Control (SEQC) AA vs. AB data at various downsampling rates.
In comparison to the non-downsampled results (as in, 1/1 downsampling), we see that moderate downsampling - e.g. at rates of 1/10, 1/20 - 
does not substantially impact DE analysis and still allow for DegNorm retains its edge over existing RNA-Seq normalization methods, while
offering 50%+ DegNorm iteration time savings.





## Skipping baseline selection
The most expensive part of the DegNorm normalization algorithm is in finding a "baseline" region of non-degraded coverage
across all samples, a process referred to as **baseline selection**. In the event that that region is small, the approximation of a single gene's coverage curves may require
over a thousand SVDs. Skipping baseline selection induces the largest speedup, more than any level of downsampling,
but likely with a more heavy impact on DegNorm's performance. This impact is experiment-specific, and is likely
due to the nature and extent of RNA degradation present in the samples. As is shown in the figure below,
the DegNorm-without-baseline-selection ROC curve in the SEQC AA vs. AB is negatively impacted over the original, while the GBM p-value ECDF
 maintains its advantage over the existing methods' curves.
 
 



#### v.0.1.1 MPI release
The next `degnorm` release will feature an MPI implementation for cross-node .bam file processing.
Further, individual DegNorm iterations will be distributed, as the per-gene NMF-OA approximation
algorithm is embarrassingly parallel. Stay tuned!
