from scipy.sparse.linalg import svds
from degnorm.utils import *
from sklearn.exceptions import NotFittedError
import warnings
import tqdm
import pickle as pkl

class GeneNMFOA():

    def __init__(self, iter=5, grid_points=None, min_high_coverage=50,
                 nmf_iter=100, bins=20, n_jobs=max_cpu()):
        """
        Initialize an NMF-over-approximator object.

        :param iter: int maximum number of NMF-OA iterations to run.
        :param grid_points: int maximum number of base pairs considered per gene coverage matrix; use to
        down-sample coverage matrices at a rate of grid_points/Li. Can also be used to effectively limit
        the maximum size of a coverage matrix considered during NMF computations.
        :param min_high_coverage: int minimum number of "high coverage" base positions for a gene
        to be considered for baseline selection algorithm. Refer to Supplement, section 2, for definition
        of high-enough coverage.
        :param nmf_iter: int number of iterations to run per NMF-OA approximation per gene's coverage matrix.
        :param bins: int number of bins to use during baseline selection step of NMF-OA loop.
        :param n_jobs: int number of cores used for distributing NMF computations over gene coverage matrices.
        """
        self.iter = np.abs(int(iter))
        self.nmf_iter = np.abs(int(nmf_iter))
        self.n_jobs = np.abs(int(n_jobs))
        self.bins = np.abs(int(bins))
        self.min_high_coverage = np.abs(int(min_high_coverage))
        self.min_bins = np.ceil(self.bins * 0.3)
        self.grid_points = np.abs(int(grid_points)) if grid_points else None
        self.downsample_idx = list()
        self.x = None
        self.x_adj = None
        self.p = None
        self.n_genes = None
        self.scale_factors = None
        self.rho = None
        self.Lis = None
        self.estimates = list()
        self.cov_mats = list()
        self.cov_mats_adj = list()
        self.fitted = False
        self.transformed = False

        # Try to prevent situation where downsampling rate is so high
        # that we can't pick high-coverage regions for baseline selection.
        if self.grid_points:
            if self.grid_points < self.min_high_coverage + 50:
                self.grid_points = self.min_high_coverage + 50

    def rank_one_approx(self, x):
        """
        Decompose a matrix X via truncated SVD into (K)(E^t) = U_{1} \cdot \sigma_{1}V_{1}

        :param x: numpy 2-d array
        :return: 2-tuple (K, E) matrix factorization
        """
        u, s, v = svds(x, k=1)
        # return u[::-1], s*v[::-1]
        return u, s*v

    def get_high_coverage_idx(self, x):
        """
        Find positions of high coverage in a gene's coverage matrix, defined
        as positions where the sample-wise maximum coverage is at least 10%
        the maximum of the entire coverage array.

        :param x: numpy 2-d array: coverage matrix (p x Li)
        :return: numpy 1-d array of integer base positions
        """
        return np.where(x.max(axis=0) > 0.1 * x.max())[0]

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

        # quality control - ensure an over-approximation.
        est[est < x] = x[est < x]

        return est

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
        # construct block rho matrix: len(self.degnorm_genes) x p matrix of DI scores
        est_sums = list(map(lambda x: x.sum(axis=1), estimates))
        self.rho = 1 - self.cov_sums / np.vstack(est_sums)

        # quality control.
        self.rho[self.rho < 0] = 0.
        self.rho[self.rho >= 1] = 1. - 1e-5

        # scale read count matrix by DI scores.
        scaled_counts_mat = self.x / (1 - self.rho)
        scaled_counts_sums = scaled_counts_mat.sum(axis=0)
        self.scale_factors = scaled_counts_sums / np.median(scaled_counts_sums)

    def adjust_coverage_curves(self):
        """
        Adjust coverage curves by dividing the coverage values by corresponding adjustment factor, s_{j}
        """
        self.cov_mats_adj = [(F.T / self.scale_factors).T for F in self.cov_mats]

    def adjust_read_counts(self):
        """
        Adjust read counts by scale factors, while being mindful of
        which genes were run through baseline algorithm and which ones were excluded.
        """

        # determine which genes will be adjusted by self.scale_factors, and which will
        # be updated by 1 - average(DI score) per sample.
        # From Supplement section 2: "For [genes with insufficient coverage], we adjust
        # the read count based on the average DI score of the corresponding sample.
        adj_standard = list(map(lambda x: len(self.get_high_coverage_idx(x)) >= self.min_high_coverage
                                , self.cov_mats_adj))
        adj_standard = np.array(adj_standard)

        # update the genes that had sufficient coverage by the estimated scale factors.
        adj_standard_idx = np.where(adj_standard)[0]
        if len(adj_standard) > 0:
            self.x_adj[adj_standard_idx, :] = self.x[adj_standard_idx, :] / (1 - self.rho[adj_standard_idx, :])

        # update the genes that had insufficient coverage by 1 - sample-average DI score
        adj_low_idx = np.where(~adj_standard)[0]
        if len(adj_low_idx) > 0:
            avg_di_score = self.rho.mean(axis=0)
            self.x_adj[adj_low_idx, :] = self.x[adj_low_idx, :] / (1 - avg_di_score)

    def matrix_sst(self, x, left_idx, right_idx):
        """
        Compute sum of squared elements for a slice of a 2-d numpy array.

        :param x: numpy 2-dimensional array
        :param left_idx: left index, beginning of slice
        :param right_idx: right index, end of slice (not inclusive)
        :return: float sum of squared values in matrix slice
        """
        return np.linalg.norm(x[:, left_idx:right_idx])

    def baseline_selection(self, F):
        """
        Find "baseline" region for a gene's coverage curves - a region where it is
        suspected that degradation is minimal, so that the coverage envelope function
        may be estimated more accurately.

        This is typically applied once coverage curves have been scaled by 1 / s_{j}.

        :param F: numpy 2-d array, gene coverage curves, a numpy 2-dimensional array
        :return: numpy 2-d array, estimate of F post baseline-selection algorithm
        """

        # determine if gene has sufficient coverage to execute baseline selection.
        # if not, do not try to estimate coverage matrix.
        hi_cov_idx = self.get_high_coverage_idx(F)
        if len(hi_cov_idx) < self.min_high_coverage:
            return F

        # select high-coverage positions from coverage matrix.
        F_bin = F[:, hi_cov_idx]
        Li = F_bin.shape[1]

        # obtain initial coverage curve estimate.
        KE_bin = self.nmf(F_bin)

        # decide where to bin up regions of the gene.
        bin_vec = np.unique(np.linspace(0, Li, num=(self.bins + 1), endpoint=True, dtype=int))
        bin_vec.sort()
        bin_segs = [[bin_vec[i], bin_vec[i + 1]] for i in range(len(bin_vec) - 1)]

        # compute DI scores.
        rho_vec = 1 - F_bin.sum(axis=1) / KE_bin.sum(axis=1)
        rho_vec[rho_vec < 0] = 0.
        rho_vec[rho_vec >= 1] = 1. - 1e-5

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

            # safely compute degradation index scores: closer to 1 -> more degradation,
            # although this should. not. be a problem given transcript filtering above.
            with np.errstate(divide='ignore', invalid='ignore'):
                rho_vec = 1 - np.true_divide(F_bin.sum(axis=1), KE_bin.sum(axis=1))
                rho_vec[~np.isfinite(rho_vec)] = 1. - 1e-5

            # quality control.
            rho_vec[rho_vec < 0] = 0.
            rho_vec[rho_vec >= 1] = 1. - 1e-5

        # Run NMF on baseline regions.
        K, E = self.nmf(F_bin, factors=True)

        # Use refined estimate of K (p x 1 vector) to refine the envelope estimate,
        # return refined estimate of KE^{T}
        E = (F.T / K.ravel()).max(axis=1).reshape(-1, 1)

        # return estimate of F.
        return K.dot(E.T)

    def run_baseline_selection_serial(self, x):
        return list(map(self.baseline_selection, x))

    def par_apply_baseline_selection(self, dat):
        """
        Apply baseline selection algorithm to all gene coverage curve matrices in parallel.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p). Ideally
        coverage arrays will be adjusted for sequencing-depth.
        :return: list of numpy reads coverage arrays estimates
        """
        # Obtain initial coverage estimates. Split estimates
        # and coverage matrices into chunked sublists.
        dat = split_into_chunks(dat, self.n_jobs)

        # Dispatch serial jobs over parallel workers.
        p = mp.Pool(processes=self.n_jobs)
        baseline_ests = [p.apply_async(self.run_baseline_selection_serial
                                  , args=(dat[i],)) for i in range(len(dat))]
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
        # return input matrix and consider everything sampled.
        if self.grid_points >= Li:
            return x, np.arange(0, Li)

        # otherwise, evenly sample column indices and return downsampled matrix.
        downsample_idx = np.sort(np.unique(np.linspace(0, Li, num=self.grid_points, endpoint=False, dtype=int)))

        if by_row:
            return x[downsample_idx, :], downsample_idx

        return x[:, downsample_idx], downsample_idx

    @staticmethod
    def restore_from_downsample(x, Li, xp, by_row=True):
        """
        Restore a coverage matrix from a downsampled version to its
        original transcript size.

        :param x: 2-d numpy.array, downsampled coverage matrix.
        :param Li: int target length of original coverage matrix prior to downsample.
        :param xp: list of int or 1-d numpy array of int, points from original which downsample was obtained.
        :return: 2-d numpy array of expanded coverage matrix.
        """
        ax = 1 if not by_row else 0

        # quality control.
        if len(xp) != x.shape[ax]:
            raise ValueError('x.shape == {0} is incompatible with downsampled'
                             'indices of length {1}'.format(x.shape, len(xp)))

        # if x has not been downsampled (shape of x is <= target length, Li), exit.
        if Li <= x.shape[ax]:
            return x

        # restore (expand-interpolate) each column or row of x at points xp with values
        # fp coming from the column or row of x.
        z = np.arange(0, Li)
        x_restored = np.apply_along_axis(lambda fp: np.interp(z
                                                              , xp=xp
                                                              , fp=fp)
                                         , axis=ax
                                         , arr=x)
        return x_restored


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

    def fit(self, cov_dat, reads_dat):
        """
        Initialize estimates for the DegNorm iteration loop.

        :param cov_dat: Dict of {gene: coverage matrix} pairs;
        coverage matrix shapes are (p x Li) (wide matrices). For each coverage matrix,
        row index is sample number, column index is gene's relative base position on chromosome.
        :param reads_dat: n (genes) x p (samples) numpy array of gene read counts
        """
        self.n_genes = len(cov_dat)
        self.genes = list(cov_dat.keys())
        self.p = next(iter(cov_dat.values())).shape[0]
        self.cov_mats = list(cov_dat.values())
        self.x = np.copy(reads_dat)
        self.x_adj = np.copy(self.x)

        # ---------------------------------------------------------------------------- #
        # Run data checks:
        # 1. Check that read count matrix has same number of rows as number of coverage arrays.
        # 2. Check that all coverage arrays are 2-d matrices.
        # 3. Check that all coverage arrays are wider than they are tall.
        # ---------------------------------------------------------------------------- #
        if self.x.shape[0] != self.n_genes:
            raise ValueError('Number of genes in read count matrix not equal to number of coverage matrices!')

        if not all(map(lambda z: z.ndim == 2, self.cov_mats)):
            raise ValueError('Not all coverage matrices are 2-d arrays!')

        if np.sum(np.array(list(map(lambda x: x.shape[1], self.cov_mats))) / self.p < 1) > 0:
            warnings.warn('At least one coverage matrix slotted for DegNorm is longer than it is wide.'
                          'Ensure that coverage matrices are (p x Li).')

        # Store array of gene transcript lengths.
        self.Lis = np.array(list(map(lambda x: x.shape[1], self.cov_mats)))

        # downsample coverage matrices (if desired), saving the indices
        # at which each matrix was sampled, so they can be restored to original size later.
        if self.grid_points:
            for i in range(self.n_genes):
                cov_mat = self.cov_mats[i]
                cov_mat, down_idx = self.downsample_2d(cov_mat, by_row=False)
                self.cov_mats[i] = cov_mat
                self.downsample_idx.append(down_idx)

        # sum coverage per sample, for later.
        self.cov_sums = list(map(lambda x: x.sum(axis=1), self.cov_mats))

        self.fitted = True

    def transform(self):
        if not self.fitted:
            raise NotFittedError('self.fit not yet executed.')

        # initialize NMF-OA estimates of coverage curves.
        # Only use on genes that meet baseline selection criteria.
        self.estimates = self.par_apply_nmf(self.cov_mats)

        # update vector of read count adjustment factors; update rho (DI) matrix.
        self.compute_scale_factors(self.estimates)

        # scale baseline_cov coverage curves by 1 / (updated adjustment factors).
        self.adjust_coverage_curves()

        # impose 1 / (1 - DI scores) scaling post self.compute_scale_factors.
        self.adjust_read_counts()

        # Instantiate progress bar.
        pbar = tqdm.tqdm(total = self.iter
                         , leave=False
                         , desc='NMF-OA iteration progress')

        # Run DegNorm iterations.
        i = 0
        while i < self.iter:

            # run NMF-OA + baseline selection; obtain refined estimates of F.
            self.estimates = self.par_apply_baseline_selection(self.cov_mats_adj)

            # update scale factors, adjust coverage curves and read counts.
            self.compute_scale_factors(self.estimates)
            self.adjust_coverage_curves()
            self.adjust_read_counts()

            i += 1
            pbar.update()

        pbar.close()

        # restore grid-sampled gene coverage curve estimates to original gene transcript sizes.
        if self.downsample_idx:
            self.estimates = list(map(lambda i: self.restore_from_downsample(self.estimates[i]
                                                                             , Li=self.Lis[i]
                                                                             , xp=self.downsample_idx[i]
                                                                             , by_row=False), range(self.n_genes)))

        self.transformed = True

    def fit_transform(self, coverage_dat, reads_dat):
        self.fit(coverage_dat, reads_dat)
        self.transform()

    def save_results(self, gene_manifest_df, output_dir='.',
                     sample_ids=None):
        """
        Save GeneNMFOA output to disk:
            - self.estimates: for each gene's estimated coverage matrix, find the chromosome to which the
            gene belongs. Per chromosome, assemble dictionary of {gene: estimated coverage matrix} pairs,
            and save that dictionary in a directory labeled by the chromosome in the output_dir directory,
            in a file "estimated_coverage_matrices_<chromosome name>.pkl"
            - self.rho: save degradation index scores to "degradation_index_scores.csv" using sample_ids
            as the header.
            - self.x_adj: save adjusted read counts to "adjusted_read_counts.csv" using sample_ids as the header.

        :param gene_manifest_df: pandas.DataFrame establishing chromosome-gene map. Must contain,
        at a minimum, the columns `chr` (str, chromosome) and `gene` (str, gene name).
        :param output_dir: str output directory to save GeneNMFOA output.
        :param sample_ids: (optional) list of str names of RNA_seq experiments, to be used
         as headers for adjusted read counts matrix and degradation index score matrix.
        """
        # quality control.
        if not self.transformed:
            raise NotFittedError('Model not fitted, not transformed -- NMF-OA has not been run.')

        if not os.path.isdir(output_dir):
            raise IOError('Directory {0} not found.'.format(output_dir))

        # ensure that gene_manifest_df has `chr` and `gene` columns to establish chromosome-gene map.
        if not all([col in gene_manifest_df.columns.tolist() for col in ['chr', 'gene']]):
            raise ValueError('gene_manifest_df must have columns `chr` and `gene`.')

        # if RNA-SEQ sample IDs were not provided, name them.
        if sample_ids:
            if len(sample_ids) != self.p:
                raise ValueError('Number of supplied sample IDs does not match number'
                                 'of samples used to fit GeneNMFOA object.')
        else:
            sample_ids = ['sample_{0}'.format(i + 1) for i in range(self.p)]

        # subset gene manifest to genes run through DegNorm.
        gene_df = gene_manifest_df[gene_manifest_df.gene.isin(self.genes)].copy()
        if gene_df.empty:
            raise ValueError('No genes used in DegNorm were found in gene manifest dataframe!')

        manifest_chroms = gene_df.chr.unique().tolist()

        # assemble chromosome-gene:estimates partitions.
        chrom_gene_dict = {chrom: dict() for chrom in manifest_chroms}
        for gene in self.genes:
            chrom = gene_df[gene_df.gene == gene].chr.iloc[0]
            if chrom not in chrom_gene_dict:
                chrom_gene_dict[chrom] = dict()

            gene_idx = np.where(np.array(self.genes) == gene)[0][0]
            chrom_gene_dict[chrom][gene] = self.estimates[gene_idx]

        # Instantiate results-save progress bar.
        pbar = tqdm.tqdm(total=(self.n_genes + 2)
                         , leave=False
                         , desc='GeneNMFOA results save progress')

        for chrom in manifest_chroms:
            chrom_dir = os.path.join(output_dir, chrom)
            if not os.path.isdir(chrom_dir):
                os.makedirs(chrom_dir)

            with open(os.path.join(chrom_dir, 'estimated_coverage_matrices_{0}.pkl'.format(chrom)), 'wb') as f:
                pkl.dump(chrom_gene_dict[chrom], f)

            pbar.update()

        np.savetxt(os.path.join(output_dir, 'degradation_index_scores.csv')
                   , X=self.rho
                   , header=','.join(sample_ids)
                   , comments=''
                   , delimiter=',')
        pbar.update()

        np.savetxt(os.path.join(output_dir, 'adjusted_read_counts.csv')
                   , X=self.x_adj
                   , header=','.join(sample_ids)
                   , comments=''
                   , delimiter=',')
        pbar.update()
        pbar.close()


# if __name__ == '__main__':
#     import pandas as pd
#     import pickle as pkl
#
#     data_path = os.path.join(os.getenv('HOME'), 'nu/jiping_research/exploration')
#     X_dat = np.load(os.path.join(data_path, 'read_counts.npz'))
#     X = X_dat['X']
#     genes_df = pd.read_csv(os.path.join(data_path, 'genes_metadata.csv'))
#     with open(os.path.join(data_path, 'gene_cov_dict.pkl'), 'rb') as f:
#         gene_cov_dict = pkl.load(f)
#
#     print('read counts matrix shape -- {0}'.format(X.shape))
#     print('genes_df shape -- {0}'.format(genes_df.shape))
#     print('number of coverage matrices -- {0}'.format(len(gene_cov_dict.values())))
#
#     print('Executing GeneNMFOA.fit_transform')
#     nmfoa = GeneNMFOA(nmf_iter=3, grid_points=1000, n_jobs=1)
#     nmfoa.fit_transform({k: v for (k, v) in list(gene_cov_dict.items())[0:100]}
#               , reads_dat=X[0:100, :])
#
#     with open(os.path.join(data_path, 'test_output.pkl'), 'wb') as f:
#         pkl.dump(nmfoa, f)
#
#     print('Saving NMF-OA results.')
#     nmfoa.save_results(gene_manifest_df=genes_df
#                        , output_dir=data_path)
#
#     print(nmfoa.rho)
#     print(nmfoa.x_adj)