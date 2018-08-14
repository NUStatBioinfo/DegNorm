from scipy.sparse.linalg import svds
from degnorm.utils import *
import warnings
import tqdm

class GeneNMFOA():

    def __init__(self, iter=5, grid_points=None,
                 tol=1e-3, nmf_iter=100, bins=20, n_jobs=max_cpu()):
        """
        Initialize an NMF-over-approximator object.

        :param iter: int maximum number of NMF-OA iterations to run.
        :param grid_points: int maximum number of base pairs considered per gene coverage matrix; use to
        down-sample coverage matrices at a rate of grid_points/Li. Can also be used to effectively limit
        the maximum size of a coverage matrix considered during NMF computations.
        :param tol: relative difference in Frobenius norms in coverage matrix approximations between
        successive NMF-OA loops before approximation stops for a particular gene.
        :param nmf_iter: int number of iterations to run per NMF-OA approximation per gene's coverage matrix.
        :param bins: int number of bins to use during baseline selection step of NMF-OA loop.
        :param n_jobs: int number of cores used for distributing NMF computations over gene coverage matrices.
        """
        self.iter = np.abs(int(iter))
        self.nmf_iter = np.abs(int(nmf_iter))
        self.tol = np.abs(tol)
        self.n_jobs = np.abs(int(n_jobs))
        self.bins = np.abs(int(bins))
        self.min_bins = np.ceil(self.bins * 0.3)
        self.grid_points = np.abs(int(grid_points)) if grid_points else None
        self.x = None
        self.p = None
        self.gene_cov_dat = None
        self.degnorm_idx = None
        self.degnorm_genes = None
        self.degnorm_coverage = None
        self.exclude_idx = None
        self.exclude_genes = None
        self.exclude_coverage = None
        self.n_genes = None
        self.degnorm_coverage_sums = None
        self.scale_factors = None
        self.rho = None

    def rank_one_approx(self, x):
        """
        Decompose a matrix X via truncated SVD into (K)(E^t) = U_{1} \cdot \sigma_{1}V_{1}

        :param x: numpy 2-d array
        :return: 2-tuple (K, E) matrix factorization
        """
        u, s, v = svds(x, k=1)
        return u[::-1], s*v[::-1]
        # return u, s*v

    def nmf(self, x, factors=False):
        """
        Run NMF-OA approximation. See "Normalization of generalized transcript degradation
        improves accuracy in RNA-seq analysis" supplement section 1.2.

        :param x: numpy 2-d array
        :param factors: boolean return 2-tuple of K, E matrix factorization? If False,
        return K.dot(E)
        :return: depending on factors, return (K, E) matrices or K.dot(E) over-approximation to x
        """
        K, E = self.rank_one_approx(x)
        est = K.dot(E)
        lmbda = np.zeros(shape=x.shape)
        c = np.sqrt(self.nmf_iter)

        for _ in range(self.nmf_iter):
            res = est - x
            lmbda -= res * (1 / c)
            lmbda[lmbda < 0.] = 0.
            K, E = np.abs(self.rank_one_approx(est + lmbda))
            est = K.dot(E)

        if factors:
            return K, E

        return est

    # def nmf_reldiff(self, x, factors=False):
    #     K, E = self.rank_one_approx(x)
    #     est0 = K.dot(E)
    #     lmbda = np.zeros(shape=x.shape)
    #     c = np.sqrt(self.nmf_iter)
    #
    #     for _ in range(self.nmf_iter):
    #         res = x - est0
    #         lmbda -= res * (1 / c)
    #         lmbda[lmbda < 0.] = 0.
    #         K, E = np.abs(self.rank_one_approx(est0 + lmbda))
    #         est1 = K.dot(E)
    #         delta = self.delta_norm(est0, est1)
    #
    #         if delta < self.tol:
    #             break
    #
    #         est0 = est1
    #
    #     if factors:
    #         return {'K': K
    #                 , 'E': E}
    #
    #     return est1

    def run_nmf_serial(self, x, factors):
        return list(map(lambda z: self.nmf(z, factors), x))

    def par_apply_nmf(self, dat):
        """
        Apply NMF-OA to a list of coverage matrices.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p)
        :return: list of dict output from self.nmf, one per numpy.array in dat
        """

        # split list of coverage matrices into chunks for distribution over n_jobs cores.
        dat_split = split_into_chunks(dat, self.n_jobs)

        p = mp.Pool(processes=self.n_jobs)
        nmf_ests = [p.apply_async(self.run_nmf_serial
                                  , args=(x, False)) for x in dat_split]
        p.close()

        # Consolidate worker results.
        nmf_ests = [x.get() for x in nmf_ests]
        return [elt for lst1d in nmf_ests for elt in lst1d]

    def compute_scale_factors(self, estimates):
        """
        Update read count scaling factors by computing n x p matrix of degradation index scores from
        a set of NMF-OA estimates.

        Compute s_{j} = \frac{sum_{i=1}^{n}\tilde{X}_{ij}}{median_{j}(\sum_{i=1}^{n}\tilde{X}_{ij})}
        vector (p x 1) of read count adjustment factors.

        :param estimates: list of rank-one approximations of coverage curves
        """
        # construct block rho matrix: len(self.degnorm_genes) x p matrix of DI scores,
        # followed by len(self.exclude_genes) x p matrix of sample-average DI score from top block.
        # From Supplement, section 2: "For [excluded genes], we adjust the read count based on the average DI score
        # of the corresponding sample."
        est_sums = list(map(lambda x: x.sum(axis=1), estimates))
        rho = 1 - self.degnorm_coverage_sums / np.vstack(est_sums)
        self.rho = np.vstack([rho, np.repeat(rho.mean(axis=0).reshape([1, -1]), self.n_genes - len(estimates), axis=0)])

        # scale read count matrix by DI scores.
        scaled_counts_mat = self.x / (1 - self.rho)
        scaled_counts_sums = scaled_counts_mat.sum(axis=0)
        self.scale_factors = scaled_counts_sums / np.median(scaled_counts_sums)

    def adjust_coverage_curves(self):
        """
        Adjust coverage curves by dividing the coverage values by corresponding adjustment factor, s_{j}

        :return: list of adjusted (scaled) numpy reads coverage arrays, shapes are (Li x p)
        """
        return [(F.T / self.scale_factors).T for F in self.degnorm_coverage]

    def matrix_sst(self, x, left_idx, right_idx):
        """
        Compute sum of squared elements for a slice of a 2-d numpy array.

        :param x: numpy 2-dimensional array
        :param left_idx: left index, beginning of slice
        :param right_idx: right index, end of slice (not inclusive)
        :return: float sum of squared values in matrix slice
        """
        return np.linalg.norm(x[:, left_idx:right_idx])

    def baseline_selection(self, F, KE):
        """
        Find "baseline" region for a gene's coverage curves - a region where it is
        suspected that degradation is minimal, so that the coverage envelope function
        may be estimated more accurately.

        :param F: gene coverage curves, a numpy 2-dimensional array
        :param KE: gene coverage curve estimate, numpy 2-dimensional array
        :return:
        """
        Li = F.shape[1]

        # Decide where to bin up regions of the gene.
        bin_vec = np.unique(np.linspace(0, Li, num=(self.bins + 1), endpoint=True, dtype=int))
        bin_vec.sort()
        bin_segs = [[bin_vec[i], bin_vec[i + 1]] for i in range(len(bin_vec) - 1)]

        KE_bin = np.copy(KE)
        F_bin = np.copy(F)

        # compute degradation index scores: closer to 1 -> more degradation.
        rho_vec = 1 - F_bin.sum(axis=1) / KE_bin.sum(axis=1)

        # Run baseline selection algorithm while the degradation index score is still high
        # or while there are still non-baseline regions left to trim.
        # If max DI score is <= 0.1, then the baseline selection algorithm has converged; if >= min_bins
        # have been dropped, still exit baseline selection, although this is not convergence in a true sense.
        while (len(bin_segs) >= self.min_bins) and (np.nanmax(rho_vec) > 0.1):

            # compute normalized residuals for each bin.
            ss_r = np.array(list(map(lambda z: self.matrix_sst(KE_bin - F_bin, z[0], z[1]), bin_segs)))
            ss_f = np.array(list(map(lambda z: self.matrix_sst(F_bin, z[0], z[1]), bin_segs)))

            # safely divide binned residual norms by original binned norms.
            with np.errstate(divide='ignore', invalid='ignore'):
                res_norms = np.true_divide(ss_r, ss_f)
                res_norms[~np.isfinite(res_norms)] = 0

            # if perfect approximation, exit loop.
            if np.nanmax(res_norms) == 0:
                break

            # drop the bin corresponding to the bin with the maximum normalized residual.
            del bin_segs[np.nanargmax(res_norms)]

            # collapse the remaining bins into an array of indices to keep, then subset F matrix.
            keep_idx = [np.arange(bin[0], bin[1], dtype=int) for bin in bin_segs]
            keep_idx = np.unique(flatten_2d(keep_idx))
            F_bin = np.copy(F[:, keep_idx])

            # estimate coverage curves from remaining regions.
            KE_bin = self.nmf(F_bin)

            # compute DI scores from remaining regions.
            rho_vec = 1 - F_bin.sum(axis=1) / KE_bin.sum(axis=1)

            # TODO: determine how to best proceed if there are divide by 0 errors.
            if any(np.isnan(rho_vec)):
                break

        # Run NMF on baseline regions.
        K, E = self.nmf(F_bin, factors=True)

        # Use refined estimate of K (p x 1 vector) to refine the envelope estimate,
        # return refined estimate of KE^{T}
        E = (F.T / K.ravel()).max(axis=1).reshape(-1, 1)

        return K.dot(E.T)

    def run_baseline_selection_serial(self, x, y):
        return list(map(self.baseline_selection, x, y))

    def par_apply_baseline_selection(self, dat):
        """
        Apply baseline selection algorithm to all gene coverage curve matrices in parallel.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p). Ideally
        coverage arrays will be adjusted for sequencing-depth.
        :return: list of numpy reads coverage arrays estimates
        """
        # Obtain initial coverage estimates. Split estimates
        # and coverage matrices into chunked sublists.
        estimates = split_into_chunks(self.par_apply_nmf(dat), self.n_jobs)
        dat = split_into_chunks(dat, self.n_jobs)

        # Dispatch serial jobs over parallel workers.
        p = mp.Pool(processes=self.n_jobs)
        baseline_ests = [p.apply_async(self.run_baseline_selection_serial
                                  , args=(dat[i], estimates[i])) for i in range(len(dat))]
        p.close()

        # Consolidate worker results.
        baseline_ests = [x.get() for x in baseline_ests]
        return [elt for lst1d in baseline_ests for elt in lst1d]

    def downsample_2d(self, x, by_row=True):
        """
        Downsample a coverage matrix at evenly spaced base position indices.

        :param x: 2-d numpy array; a coverage matrix
        :param by_row: bool indicator - downsample by rows or columns?
        :return: 2-d numpy array with self.grid_points columns
        """
        Li = x.shape[0 if by_row else 1]

        # if downsample rate high enough to cover entire matrix,
        # just return matrix.
        if self.grid_points >= Li:
            return x

        # otherwise, evenly sample column indices and return downsampled matrix.
        downsample_idx = np.sort(np.unique(np.linspace(0, Li, num=self.grid_points, endpoint=False, dtype=int)))

        if by_row:
            return x[downsample_idx, :]

        return x[:, downsample_idx]

    def delta_norm(self, mat_t0, mat_t1):
        """
        Compute the relative difference between two matrices (one at time = 0,
         the other at time = 1) by using a Frobenius norm:
        ||X0 - X1||_{F}^{2} / ||X0||_{F}^{2}.

        If X0 is a matrix of 0's, return ||X1||_{F}^{2}.

        :param mat_t0: numpy nd array, represents matrix at iteration time = 0
        :param mat_t1: numpy nd array, represents matrix at iteration time = 1
        :return: float relative change in Frobenius norms between time 0 and time 1
        """
        t0_norm = np.linalg.norm(mat_t0)
        diff_norm = np.linalg.norm(mat_t0 - mat_t1)

        if t0_norm == 0:
            return diff_norm

        return diff_norm / t0_norm

    def apply_delta_norm(self, mats_t0, mats_t1):
        """
        Apply self.delta_norm to a list of matrices at time = 0 and a list of matrices at
        time = 1.

        :param mats_t0: list of numpy nd arrays at time 0
        :param mats_t1: list of numpy nd arrays at time 0
        :return: 1-d numpy array of relative changes in Frobenius norms between corresponding
        matrices in mats_t0 and mats_t1
        """
        if len(mats_t0) != len(mats_t1):
            raise ValueError('len(matsT0) != len(matsT1); cannot find '
                             'relative difference in matrix norms.')
        return np.array(list(map(self.delta_norm, mats_t0, mats_t1)))

    def fit(self, gene_cov_dat, reads_dat):
        """
        Initialize estimates for the DegNorm iteration loop.

        :param gene_cov_dat: Dict of {gene: coverage matrix} pairs;
        coverage matrix shapes are (p x Li) (wide matrices). For each coverage matrix,
        row index is sample number, column index is gene's relative base position on chromosome.
        :param reads_dat: n (genes) x p (samples) numpy array of gene read counts
        """
        self.gene_cov_dat = gene_cov_dat
        self.n_genes = len(gene_cov_dat)

        # split apart genes to be run in DegNorm from those that should not.
        self.degnorm_idx = np.where([self.gene_cov_dat[gene]['run_degnorm'] for gene in self.gene_cov_dat])[0]
        self.degnorm_genes = np.array(list(self.gene_cov_dat.keys()))[self.degnorm_idx]
        self.degnorm_coverage = [self.gene_cov_dat[gene]['filtered_coverage'] for gene in self.degnorm_genes]
        self.p = self.degnorm_coverage[0].shape[0]

        if len(self.degnorm_genes) < self.n_genes:
            self.exclude_idx = np.array(list(set(range(self.n_genes)) - set(self.degnorm_idx)))
            self.exclude_genes = np.array(list(self.gene_cov_dat.keys()))[self.exclude_idx]
            self.exclude_coverage = [self.gene_cov_dat[gene]['raw_coverage'] for gene in self.exclude_genes]

        # order read counts into two blocks: top block for genes to be run in DegNorm.
        self.x = np.copy(reads_dat[np.concatenate([self.degnorm_idx, self.exclude_idx]), :])

        # ---------------------------------------------------------------------------- #
        # Run data checks:
        # 1. Check that read count matrix has same number of rows as number of coverage arrays.
        # 2. Check that all coverage arrays are 2-d matrices.
        # 3. Check that all coverage arrays are wider than they are tall.
        # 4. Check for genes to be run in DegNorm that have experiments with entirely-zero coverage.
        # ---------------------------------------------------------------------------- #
        if self.x.shape[0] != self.n_genes:
            raise ValueError('Number of genes in read count matrix not equal to number of coverage matrices!')

        if not all(map(lambda z: z.ndim == 2, self.degnorm_coverage + self.exclude_coverage)):
            raise ValueError('Not all coverage matrices are 2-d arrays!')

        if np.sum(np.array(list(map(lambda x: x.shape[1], self.degnorm_coverage))) / self.p < 1) > 0:
            warnings.warn('At least one coverage matrix slotted for DegNorm is longer than it is wide.'
                          'Ensure that coverage matrices are (p x Li).')

        zero_coverage_idx = np.where(np.array(list(map(lambda x: any(x.sum(axis=1) == 0)
                                                       , self.degnorm_coverage))))[0].astype(int)
        if len(zero_coverage_idx) > 0:
            raise ValueError('The following genes to run in DegNorm have experiments with zero coverage: {0}'
                             .format('\n\t'.join([self.degnorm_genes[idx] for idx in zero_coverage_idx])))

        # downsample coverage matrices if desired.
        if self.grid_points:
            self.degnorm_coverage = list(map(lambda z: self.downsample_2d(z, by_row=False), self.degnorm_coverage))

        # assemble len(self.degnorm_genes) x p matrix of sums over coverage arrays
        # (sum coverage over positions, one sum per sample included in DegNorm iterations).
        self.degnorm_coverage_sums = np.array(list(map(lambda z: z.sum(axis=1), self.degnorm_coverage)))

    def transform(self):
        # initialize NMF-OA estimates of coverage curves
        estimates = self.par_apply_nmf(self.degnorm_coverage)

        # Instantiate progress bar.
        pbar = tqdm.tqdm(total = self.iter
                         , leave=False
                         , desc='NMF-OA iteration progress')

        # Run DegNorm iterations.
        i = 0
        while i < self.iter:

            # update vector of read count adjustment factors
            self.compute_scale_factors(estimates)

            # scale coverage curves by 1 / adjustment factors
            coverage_adjusted = self.adjust_coverage_curves()

            # run NMF-OA + baseline selection; obtain refined estimate of K*Et
            estimates = self.par_apply_baseline_selection(coverage_adjusted)

            # # determine which genes coverage curve estimates have converged.
            # delta_vec = self.apply_delta_norm(estimates_t0
            #                                   , mats_t1=estimates_t1)
            # converged_idx = np.where(delta_vec <= self.tol)[0]
            # non_converged_idx = np.where(delta_vec > self.tol)[0]

            # # store the converged genes' data and gene index.
            # if len(converged_idx) > 0:
            #     for idx in converged_idx:
            #         gene_dict[self.genes[idx]] = dict()
            #         gene_dict[self.genes[idx]]['estimate'] = estimates_t1[idx]
            #         gene_dict[self.genes[idx]]['converged'] = True

            # if all genes have converged, yay exit and continue.
            # if len(non_converged_idx) == 0:
            #     pbar.update(self.iter - i)
            #     break

            # drop the converged genes' coverage matrices, coverage sums, and estimates for next iteration.
            # else:
            #     self.coverage_dat = [self.coverage_dat[idx] for idx in non_converged_idx]
            #     self.genes = [self.genes[idx] for idx in non_converged_idx]
            #     estimates_t0 = [estimates_t1[idx] for idx in non_converged_idx]
            #     self.coverage_sums = self.coverage_sums[non_converged_idx]
            #     self.x = np.delete(self.x
            #                        , obj=converged_idx
            #                        , axis=0)

            i += 1
            pbar.update()

        # # additionally, store all of the non-converged genes.
        # if len(non_converged_idx) > 0:
        #     for idx in non_converged_idx:
        #         gene_dict[self.genes[idx]] = dict()
        #         gene_dict[self.genes[idx]]['estimate'] = estimates_t1[idx]
        #         gene_dict[self.genes[idx]]['converged'] = False

        pbar.close()

        for i in range(len(self.degnorm_genes)):
            self.gene_cov_dat[self.degnorm_genes[i]]['estimate'] = estimates[i]

    def fit_transform(self, coverage_dat, reads_dat):
        self.fit(coverage_dat, reads_dat)
        return self.transform()


if __name__ == '__main__':
    import pandas as pd
    import pickle as pkl

    data_path = os.path.join(os.getenv('HOME'), 'nu/jiping_research/exploration')
    X_dat = np.load(os.path.join(data_path, 'read_counts.npz'))
    X = X_dat['X']
    genes_df = pd.read_csv(os.path.join(data_path, 'genes_metadata.csv'))
    with open(os.path.join(data_path, 'gene_cov_dict.pkl'), 'rb') as f:
        gene_cov_dict = pkl.load(f)

    print('read counts matrix shape -- {0}'.format(X.shape))
    print('genes_df shape -- {0}'.format(genes_df.shape))
    print('number of coverage matrices -- {0}'.format(len(gene_cov_dict.values())))

    nmfoa = GeneNMFOA(nmf_iter=30, grid_points=5000, n_jobs=1)
    nmfoa.fit({k: v for (k, v) in list(gene_cov_dict.items())[0:100]}
              , reads_dat=X[0:100, :])
    print('Successfully fitted GeneNMFOA object.')

    print('Running GeneNMFOA.transform()')
    nmfoa.transform()

    import pickle as pkl
    with open(os.path.join(data_path, 'test_output.pkl'), 'wb') as f:
        pkl.dump(nmfoa.gene_cov_dat, f)

    estimates = [nmfoa.gene_cov_dat[gene]['estimate'] for gene in nmfoa.degnorm_genes]

    # adjust read counts based on DI scores computed from final coverage estimates.
    nmfoa.compute_scale_factors(estimates)
    rho = nmfoa.rho
    X_adj = nmfoa.x / (1 - rho)
    print(X_adj)