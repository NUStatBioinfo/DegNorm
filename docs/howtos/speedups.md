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
 
 


#### v.0.1.1 MPI release
The next `degnorm` release will feature an MPI implementation for cross-node .bam file processing.
Further, individual DegNorm iterations will be distributed, as the per-gene NMF-OA approximation
algorithm is embarrassingly parallel. Stay tuned!
