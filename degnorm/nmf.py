from scipy.sparse.linalg import svds
from pandas import DataFrame, concat
from degnorm.utils import *
import warnings
import tqdm
import pickle as pkl
from joblib import Parallel, delayed


class GeneNMFOA():

    def __init__(self, degnorm_iter=5, downsample_rate=1, min_high_coverage=50,
                 nmf_iter=100, bins=20, n_jobs=max_cpu(), skip_baseline_selection=False, random_state=123):
        """
        Initialize an NMF-over-approximator object.

        :param degnorm_iter: int maximum number of NMF-OA iterations to run.
        :param downsample_rate: int reciprocal of downsample rate; a "take every" systematic sampling rate. Use to
        down-sample coverage matrices by taking every r-th nucleotide. When downsample_rate == 1 this is equivalent
        to running DegNorm sans any downsampling.
        :param min_high_coverage: int minimum number of "high coverage" base positions for a gene
        to be considered for baseline selection algorithm. Refer to Supplement, section 2, for definition
        of high-enough coverage.
        :param nmf_iter: int number of iterations to run per NMF-OA approximation per gene's coverage matrix.
        :param bins: int number of bins to use during baseline selection step of NMF-OA loop.
        :param n_jobs: int number of cores used for distributing NMF computations over gene coverage matrices.
        :param skip_baseline_selection: Boolean should DegNorm skip baseline selection process?
        :param random_state: int seed for random number generator, useful if downsampling coverage matrices.
        """
        self.degnorm_iter = np.abs(int(degnorm_iter))
        self.nmf_iter = np.abs(int(nmf_iter))
        self.n_jobs = np.abs(int(n_jobs))
        self.bins = np.abs(int(bins))
        self.min_high_coverage = max(2, np.abs(int(min_high_coverage)))
        self.min_bins = np.ceil(self.bins * 0.2)
        self.downsample_rate = np.abs(int(downsample_rate))
        self.mem_splits = None
        self.x = None
        self.x_weighted = None
        self.x_adj = None
        self.p = None
        self.n_genes = None
        self.norm_factors = None
        self.scale_factors = None
        self.rho = None
        self.fitted = False
        self.ran_baseline_selection = None
        self.skip_baseline_selection = skip_baseline_selection
        self.random_state = random_state

        # all coverage matrices must have >= 2 high-coverage indices if downsampling (svds function limitation).
        if self.downsample_rate > 1:
            self.min_high_coverage = 2

    @staticmethod
    def rank_one_approx(x):
        """
        Decompose a matrix X via truncated SVD into (K)(E^t) = U_{1} \cdot \sigma_{1}V_{1}

        :param x: numpy 2-d array
        :return: 2-tuple (K, E) matrix factorization
        """
        u, s, v = svds(x, k=1)
        return u*s, v

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
        c = 1. / np.sqrt(self.nmf_iter)

        for _ in range(self.nmf_iter):
            res = est - x
            lmbda -= c * res
            lmbda[lmbda < 0.] = 0.
            K, E = self.rank_one_approx(x + lmbda)
            est = K.dot(E)

        if factors:
            # return np.abs(K), np.abs(E)
            return K, E

        # quality control - ensure an over-approximation.
        # est[est < x] = x[est < x]

        return est

    def ratio_svd(self, x):
        """
        One-iteration SVD over-approximation, but not the NMFOA algorithm.
        See https://bit.ly/2zR4XEn.

        :param x: 2-d numpy array
        :return: 2-d numpy array estimate, elements are at least as large as those in x.
        """
        K, E = self.rank_one_approx(x)
        est = K.dot(E)
        est[est < x] = x[est < x]

        return est

    def run_ratio_svd_serial(self, x):
        return list(map(self.ratio_svd, x))

    def par_apply(self, fun, dat):
        """
        Apply a function to a list in parallel with threading.

        :param fun: function to be applied to elements of list.
        :param dat: list of elements to be chunked up, distributed to threads, and hit with fun.
        :return: list of output from fun
        """
        # split up list of coverage matrices to feed to workers.
        dat = split_into_chunks(dat, self.mem_splits)
        par_output = Parallel(n_jobs=self.n_jobs
                              , verbose=0
                              , backend='threading')(map(delayed(fun), dat))

        return [est for est1d in par_output for est in est1d]

    def adjust_coverage_curves(self, dat):
        """
        Adjust coverage curves by dividing the coverage values by corresponding scale factor, s_{j}
        """
        return [(F.T / self.scale_factors).T for F in dat]

    def correct_di_scores(self, save=False):
        """
        Correct degradation index (DI) scores, i.e. the rho matrix, by assigning the sample average DI
        score to genes that were not run through baseline selection.
        """

        # determine indices of genes that were not sent through baseline selection; if they exist, do the swap.
        non_baseline_gene = self.rho.max(axis=1) == 0
        if np.sum(non_baseline_gene) > 0:
            sample_avg_di_scores = 1 - (self.x_weighted.sum(axis=0) / self.x_adj.sum(axis=0))
            self.rho[non_baseline_gene, :] = sample_avg_di_scores

    @staticmethod
    def shift_bins(bins, dropped_bin):
        """
        Shift bins of consecutive numbers after a bin has been dropped.
        Example:

        shift_bins([[0, 1], [2, 3], [4, 5], [6, 7]], dropped_bin=1)

        [[0, 1], [2, 3], [4, 5]]

        :param bins: list of list of int
        :param dropped_bin: int, the bin that has been deleted from bins.
        :return: list of consecutive integer bins, preserving to the fullest extent
        possible the original bins.
        """
        if (dropped_bin == len(bins)) or (len(bins) == 1):
            return bins
        else:
            if dropped_bin == 0:
                delta = bins[0][0]

            else:
                delta = bins[dropped_bin][0] - bins[(dropped_bin - 1)][-1] - 1

            for bin_idx in range(dropped_bin, len(bins)):
                bins[bin_idx] = [idx - delta for idx in bins[bin_idx]]

        return bins

    def baseline_selection(self, F):
        """
        Find "baseline" region for a gene's coverage curves - a region where it is
        suspected that degradation is minimal, so that the coverage envelope function
        may be estimated more accurately, in addition to finding the gene's DI scores
        across the samples.

        This is typically applied once coverage curves have been scaled by 1 / s_{j}.

        :param F: numpy 2-d array, gene coverage curve matrix, a numpy 2-dimensional array, perferrably
        the coverage curves have been scaled by the degradation normalization scale factor.
        :return: 3-tuple --
        (numpy 1-d array (DI scores)
        , numpy 2-d array (estimate of F post baseline-selection algorithm)
        , Boolean indicator of whether gene was sent through baseline selection)
        """

        # ------------------------------------------------------------------------- #
        # Overhead:
        # 1. Downsample if desired.
        # 2. Find high-coverage regions of the gene's transcript.
        # 3. If not downsampling, determine whether or not there are sufficient
        #   high-coverage nucleotide positions found in 2.
        # 4. Subset coverage curve to (possibly downsampled) high-coverage regions.
        # 5. Bin-up consecutive regions of high-coverage regions.
        # ------------------------------------------------------------------------- #
        output = dict()
        output['rho'] = np.zeros(self.p)  # default degradation scores are all 0.
        output['estimate'] = F
        output['ran_baseline_selection'] = False

        hi_cov_idx = self.get_high_coverage_idx(F)

        # Run systematic downsample if desired, prior to filtering out low-coverage regions.
        if self.downsample_rate > 1:
            _, downsample_idx = self.downsample_2d(F, by_row=False)

            # intersect sampled indices with the high-coverage indices.
            hi_cov_idx = np.intersect1d(downsample_idx, hi_cov_idx)

        n_hi_cov = len(hi_cov_idx)

        # if there are not sufficiently-many high-coverage base pairs, return null degradation.
        if n_hi_cov < self.min_high_coverage:
            return output['rho'], output['estimate'], output['ran_baseline_selection']

        # select high-coverage positions from (possibly downsampled) coverage matrix.
        hi_cov_idx.sort()
        F_start = np.copy(F[:, hi_cov_idx])
        F_bin = np.copy(F_start)

        # if any sample sans coverage after filtering, return defaults.
        if np.sum(F_bin.sum(axis=1) > 0) < self.p:
            return output['rho'], output['estimate'], output['ran_baseline_selection']

        # run NMF on filtered coverage, obtain initial coverage curve estimate.
        K, E = self.nmf(F_bin
                        , factors=True)
        KE_bin = K.dot(E)

        # keep original NMFOA-estimated coverage in case we do not run baseline selection.
        K_start, E_start = np.copy(K), np.copy(E)
        estimate = np.copy(KE_bin)

        # compute DI scores. High score ->> high degradation.
        rho_vec = 1 - F_bin.sum(axis=1) / (KE_bin.sum(axis=1) + 1)

        # exclude extreme cases where NMF result doesn't converge.
        if np.nanmedian(1 - rho_vec) > 1:
            return output['rho'], output['estimate'], output['ran_baseline_selection']

        # decide whether to search for a gene's baseline regions.
        # If any of these criteria are satisfied, do not run baseline selection and return what we have.
        if (n_hi_cov >= 200) and (np.nanmin(rho_vec) <= 0.2) and (not self.skip_baseline_selection):

            # split up consecutive regions of the high-coverage gene.
            # even in downsampling regime, drop batches of sample points on each baseline selection iteration.
            bin_segs = split_into_chunks(list(range(F_bin.shape[1]))
                                         , n=self.bins)
            n_bins = len(bin_segs)

            while np.nanmax(rho_vec) > 0.1:

                # identify gene as having gone through baseline selection.
                output['ran_baseline_selection'] = True

                # compute relative residuals from NMF output,
                # then compute average weighted squared relative residual.
                res_vec = np.apply_along_axis(lambda z: np.nanmax(z ** 2)
                                              , axis=0
                                              , arr=(KE_bin - F_bin) / (F_bin + 1))
                ss_r = np.array(list(map(lambda idx: np.nanmean(res_vec[bin_segs[idx]]), range(n_bins))))

                # if perfect approximation, exit loop.
                if np.nanmax(ss_r) == 0:
                    break

                # drop the bin corresponding to the bin with the maximum normalized residual,
                # along with corresponding indices of F.
                drop_idx = np.nanargmax(ss_r)
                F_bin = np.delete(F_bin
                                  , obj=bin_segs[drop_idx]
                                  , axis=1)
                del bin_segs[drop_idx]

                n_hi_cov = F_bin.shape[1]

                # update binned segments so that all bins are referencing indices in newly shrunken coverage matrices.
                bin_segs = self.shift_bins(bin_segs
                                           , dropped_bin=drop_idx)
                n_bins = len(bin_segs)

                # shrink F matrix to the indices not dropped, cast back to wide matrix.
                # estimate coverage curves from said indices with NMFOA.
                K, E = self.nmf(F_bin
                                , factors=True)
                KE_bin = K.dot(E)

                # stop is fitted values are all zero for any sample (extreme cases).
                if np.min(KE_bin.sum(axis=1)) == 0:
                    break

                KE_bin[KE_bin < F_bin] = F_bin[KE_bin < F_bin]

                # recompute DI scores: closer to 1 ->> more degradation. Then run quality control.
                rho_vec = 1 - F_bin.sum(axis=1) / (KE_bin.sum(axis=1) + 1)  # + 1 as per Bin's code.

                if (n_bins <= self.min_bins) or (n_hi_cov < 200):
                    break

            # determine whether a baseline region has been identified.
            if np.nanmax(rho_vec) < 0.2:  # baseline region successfully converged upon
                # quality control: ensure we never divide F by 0.
                K = np.abs(K)
                K[K < 1.e-5] = np.min(K[K >= 1.e-5])

                # Use refined estimate of K (p x 1 vector) to refine the envelope estimate.
                E = np.true_divide(F_start.T, K.ravel()).max(axis=1).reshape(-1, 1).T
                estimate = K.dot(E)

                # update DI scores based on whole filtered gene transcript.
                rho_vec = 1 - F_start.sum(axis=1) / (estimate.sum(axis=1) + 1)

                # eliminate extreme cases with high degradation but where it's unlikely
                # that this degradation is more than just noise. Typically happens for long
                # genes with low sequencing depth.
                if np.nanmax(rho_vec) > 0.9:
                    K, E = K_start, E_start
                    estimate = K.dot(E)
                    estimate[estimate < F_start] = F_start[estimate < F_start]
                    rho_vec = 1 - F_start.sum(axis=1) / (estimate.sum(axis=1) + 1)

            # otherwise when baseline has not been identified, use the original NMF results.
            else:
                K, E = K_start, E_start
                estimate = K.dot(E)
                estimate[estimate < F_start] = F_start[estimate < F_start]
                rho_vec = 1 - F_start.sum(axis=1) / (estimate.sum(axis=1) + 1)

        # problem: estimate most likely not the same length as original gene transcript if baseline selection ran.
        # current solution: use best K estimate to back out E = F / K; use K^{T}E as estimate for visual purposes.
        # TODO: determine best way get final estimate of entire transcript after baseline selection.
        if estimate.shape[1] < F.shape[1]:

            # quality control: ensure we never divide F by 0.
            K = np.abs(K)
            K[K < 1.e-5] = np.min(K[K >= 1.e-5])
            E = np.true_divide(F.T, K.ravel()).max(axis=1).reshape(-1, 1).T
            estimate = K.dot(E)
            estimate[estimate < F] = F[estimate < F]

        # assemble final output data.
        output['rho'] = rho_vec
        output['estimate'] = estimate

        # return DI vector, estimate of F, True/False whether gene sent through baseline selection.
        return output['rho'], output['estimate'], output['ran_baseline_selection']

    def run_baseline_selection_serial(self, x):
        return list(map(self.baseline_selection, x))

    def par_apply_baseline_selection(self, dat, degnorm_iter):
        """
        Apply baseline selection algorithm to all gene coverage matrices in parallel,
        building the DI score matrix (rho) in the process and
        returning the estimated coverage matrices for visualization purposes.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p). Ideally
        coverage arrays will have been already adjusted for sequencing-depth by scale factors.
        :return: list of 2-d numpy arrays used to visualize DegNorm estimated coverage matrices.
        """
        # split up coverage matrices so that no worker gets much more than 50Mb.
        dat = split_into_chunks(dat, self.mem_splits)
        baseline_dat = Parallel(n_jobs=self.n_jobs
                                , verbose=0
                                , backend='threading')(map(delayed(self.run_baseline_selection_serial), dat))

        # flatten results.
        baseline_dat = [est for est1d in baseline_dat for est in est1d]

        # Assemble DI score matrix, apply QA thresholding on output.
        self.rho = np.vstack([x[0] for x in baseline_dat])
        self.rho[self.rho > 0.9] = 0.9
        self.rho[self.rho < 0.] = 0.

        # Update current column of baseline selection tracker matrix for which genes were
        # sent through baseline selection algorithm.
        self.ran_baseline_selection[:, degnorm_iter] = np.array([x[2] for x in baseline_dat])

        # return coverage curve estimates; for visualization purposes (not currently for DI score calculations).
        return [x[1] for x in baseline_dat]

    @staticmethod
    def _systematic_sample(n, take_every):
        """
        Take an honest to goodness systematic sample of n objects using a take_every gap.
        Randomly initialize a starting point < take every.

        :param n: int size of population to sample from
        :param take_every: int, preferably < n
        :return: 1-d numpy array of sampled integers, length ceiling(n / take_every)
        """
        # if downsample rate larger than n, randomly sample a single point.
        if take_every >= n:
            return int(np.random.choice(n))

        start = np.random.choice(take_every)
        sample_idx = np.arange(start, n
                               , step=take_every
                               , dtype=int)
        return sample_idx

    def downsample_2d(self, x, by_row=True):
        """
        Downsample a coverage matrix at evenly spaced base position indices via systematic sampling.

        :param x: 2-d numpy array; a coverage matrix
        :param by_row: bool indicator - downsample by rows or columns?
        :return: 2-d numpy array with self.grid_points columns
        """
        Li = x.shape[0 if by_row else 1]

        # if downsample rate high enough to cover entire matrix,
        # return input matrix and consider everything sampled.
        if self.downsample_rate == 1:
            return x, np.arange(0, Li)

        if self.downsample_rate >= Li:
            raise ValueError('Cannot downsample at a rate < 1 / length(gene)')

        # obtain systematic sample of transcript positions.
        downsample_idx = self._systematic_sample(Li
                                                 , take_every=self.downsample_rate)

        if by_row:
            return x[downsample_idx, :], downsample_idx

        return x[:, downsample_idx], downsample_idx

    def check_input(self, cov_mats):
        """
        Run data checks:
        1. Check that read count matrix has same number of rows as number of coverage arrays.
        2. Check that all coverage arrays are 2-d matrices.
        3. Check that all coverage arrays are wider than they are tall.
        4. Check that downsample rate < length(gene) for all genes, if downsampling.

        :param cov_mats: list of gene coverage matrices, i.e. p x Li 2-d numpy arrays
        :return: None
        """
        # Store array of gene transcript lengths.
        li_vec = np.array(list(map(lambda x: x.shape[1], cov_mats)))

        if self.x.shape[0] != self.n_genes:
            raise ValueError('Number of genes in read count matrix not equal to number of coverage matrices!')

        if not all(map(lambda z: z.ndim == 2, cov_mats)):
            raise ValueError('Not all coverage matrices are 2-d arrays!')

        if np.sum(li_vec / self.p < 1) > 0:
            warnings.warn('At least one coverage matrix is longer than it is wide.'
                          'Ensure that coverage matrices are (p x Li).')

        if self.downsample_rate > 1:
            if not np.min(li_vec) >= self.downsample_rate:
                raise ValueError('downsample_rate is too large; take-every size > at least one gene.')

    def run(self, cov_dat, reads_dat):
        """
        Run DegNorm degradation normalization pipeline: adjust read counts, compute degradation index scores,
        and compute normalized coverage curve estimates.

        :param cov_dat: OrderedDict of {gene: coverage matrix} pairs;
        coverage matrix shapes are (p x Li) (wide matrices). For each coverage matrix,
        row index is sample number, column index is gene's relative base position on chromosome.
        :param reads_dat: n (genes) x p (samples) numpy array of gene read counts

        :return: list of 2-d numpy arrays, estimated coverage matrices. In same order as self.genes.
        """
        self.n_genes = len(cov_dat)
        self.genes = list(cov_dat.keys())
        self.x = np.copy(reads_dat)

        cov_mats = list(cov_dat.values())
        self.p = cov_mats[0].shape[0]

        # initialize baseline selection tracker entirely False (no genes have gone through baseline selection yet).
        self.ran_baseline_selection = np.zeros(shape=[self.n_genes, self.degnorm_iter]).astype(bool)

        # check validity of input data.
        out = self.check_input(cov_mats)

        # determine (integer) number of data splits for parallel workers (50Mb per worker).
        mem_splits = int(np.ceil(np.sum(list(map(lambda x: x.nbytes, cov_mats))) / 5e7))
        self.mem_splits = mem_splits if mem_splits > self.n_jobs else self.n_jobs

        # ---------------------------------------------------------------------------- #
        # DegNorm initialization:
        # 1. Estimate initial scale factors from single-iteration rank-1 SVD estimates,
        #    obtain initial DI scores.
        # 2. Use initial DI scores to compute initial normalization factors.
        # 3. Normalize read counts by initial normalization factors.
        # 4. Initialize coverage scale factors with initial normalization factors.
        # ---------------------------------------------------------------------------- #

        # use self.ratio_svd to obtain first coverage matrix estimates, use to compute initial DI scores.
        estimates = self.par_apply(fun=self.run_ratio_svd_serial
                                   , dat=cov_mats)
        est_sums = np.vstack(list(map(lambda x: x.sum(axis=1), estimates)))
        cov_sums = np.vstack(list(map(lambda x: x.sum(axis=1), cov_mats)))
        self.rho = 1 - (cov_sums / (est_sums + 1))

        # estimate normalization factors from initial DI scores.
        low_di_gene = self.rho.max(axis=1) < 0.1
        count_sums = self.x[low_di_gene, :].sum(axis=0) if np.any(low_di_gene) else self.x.sum(axis=0)
        self.norm_factors = count_sums / np.median(count_sums)

        # adjust read counts by initial normalization factors, set initial scale factors.
        self.x_weighted = self.x / self.norm_factors
        self.scale_factors = np.copy(self.norm_factors)

        logging.info('Initial sequencing depth scale factors -- \n\t{0}'
                     .format(', '.join([str(x) for x in self.scale_factors])))

        # ---------------------------------------------------------------------------- #
        # DegNorm iterations:
        # 1. Scale coverage curves by scale factors from iteration t-1.
        # 2. Run baseline selection, obtain degradation-normalized over-approximated coverage curves.
        # 3. Compute DI scores from NMF-OA estimates, treating genes that went through
        #    baseline selection differently than those who did not.
        # 4. Re-normalize read counts based on DI scores at time t.
        # 5. Update sequencing depth scale factors based on re-normalized read counts.
        # ---------------------------------------------------------------------------- #

        # instantiate progress bar.
        pbar = tqdm.tqdm(total=self.degnorm_iter
                         , leave=False
                         , desc='NMF-OA iteration progress')

        # set random number generator seed.
        np.random.seed(self.random_state)

        # Run DegNorm iterations.
        i = 0
        while i < self.degnorm_iter:

            # scale baseline_cov coverage curves by 1 / (updated adjustment factors).
            cov_mats_adj = self.adjust_coverage_curves(cov_mats)

            # run NMF-OA + baseline selection; obtain refined estimates of F; update rho matrix.
            estimates = self.par_apply_baseline_selection(cov_mats_adj
                                                          , degnorm_iter=i)

            # declare number of genes sent through baseline selection on this iteration.
            if not self.skip_baseline_selection:
                logging.info('DegNorm iteration {0} -- {1} genes sent through baseline selection'
                             .format(i + 1, np.sum(self.ran_baseline_selection[:, i])))

            # adjust (weighted) read counts.
            self.x_adj = self.x_weighted / (1 - self.rho)

            # update scale factors, adjust read counts.
            self.correct_di_scores()

            # adjust norm-factor-weighted read counts by DI scores.
            self.x_adj = self.x_weighted / (1 - self.rho)

            # get new norm factors.
            self.norm_factors = self.x_adj.sum(axis=0) / np.median(self.x_adj.sum(axis=0))

            # update read counts by adjusting degradation effect into sequencing depth
            self.x_weighted = self.x_weighted / self.norm_factors

            # increment new scale_factors.
            self.scale_factors = self.scale_factors * self.norm_factors

            logging.info('DegNorm iteration {0} -- sequencing depth scale factors: \n\t{1}'
                         .format(i + 1, ', '.join([str(x) for x in self.scale_factors])))

            i += 1
            pbar.update()

        pbar.close()
        self.fitted = True

        return estimates

    def save_results(self, estimates, gene_manifest_df,
                     output_dir='.', sample_ids=None):
        """
        Save GeneNMFOA output to disk:
            - estimates: for each gene's estimated coverage matrix, find the chromosome to which the
            gene belongs. Per chromosome, assemble dictionary of {gene: estimated coverage matrix} pairs,
            and save that dictionary in a directory labeled by the chromosome in the output_dir directory,
            in a file "estimated_coverage_matrices_<chromosome name>.pkl"
            - self.rho: save degradation index scores to "degradation_index_scores.csv" using sample_ids
            as the header.
            - self.x_adj: save adjusted read counts to "adjusted_read_counts.csv" using sample_ids as the header.

        :param estimates: list of 2-d numpy arrays, estimated coverage matrices. In same order as self.genes.
        :param gene_manifest_df: pandas.DataFrame establishing chromosome-gene map. Must contain,
        at a minimum, the columns `chr` (str, chromosome) and `gene` (str, gene name).
        :param output_dir: str output directory to save GeneNMFOA output.
        :param sample_ids: (optional) list of str names of RNA_seq experiments, to be used
         as headers for adjusted read counts matrix and degradation index score matrix.
        """
        # quality control.
        if not self.fitted:
            raise ValueError('Model not yet fit. NMF-OA has not been run.')

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

        gene_intersect = np.intersect1d(gene_manifest_df.gene.unique(), self.genes)
        if len(gene_intersect) < len(self.genes):
            warnings.warn('Gene manifest data does not encompass set of genes sent through DegNorm.')

        # subset gene manifest to genes run through DegNorm.
        if len(gene_intersect) == 0:
            raise ValueError('No genes used in DegNorm were found in gene manifest dataframe!')

        gene_df = gene_manifest_df[gene_manifest_df.gene.isin(gene_intersect)]
        manifest_chroms = gene_df.chr.unique().tolist()

        # nest per-gene coverage curve estimates within chromosomes.
        chrom_gene_dict = {chrom: dict() for chrom in manifest_chroms}
        for gene in gene_intersect:
            chrom = gene_df[gene_df.gene == gene].chr.iloc[0]
            if chrom not in chrom_gene_dict:
                chrom_gene_dict[chrom] = dict()

            gene_idx = np.where(np.array(self.genes) == gene)[0][0]
            chrom_gene_dict[chrom][gene] = estimates[gene_idx]

        # save estimated coverage matrices to genes nested within chromosomes.
        chrom_gene_dfs = list()
        for chrom in manifest_chroms:
            chrom_dir = os.path.join(output_dir, chrom)

            if not os.path.isdir(chrom_dir):
                os.makedirs(chrom_dir)

            with open(os.path.join(chrom_dir, 'estimated_coverage_matrices_{0}.pkl'.format(chrom)), 'wb') as f:
                pkl.dump(chrom_gene_dict[chrom], f)

            # keep track of chromosome-gene mapping.
            chrom_gene_dfs.append(DataFrame({'chr': chrom
                                             , 'gene': list(chrom_gene_dict[chrom].keys())}))

        # order chromosome-gene index according to nmfoa gene order.
        chrom_gene_df = concat(chrom_gene_dfs)
        chrom_gene_df.set_index('gene'
                                , inplace=True)
        chrom_gene_df = chrom_gene_df.loc[self.genes]
        chrom_gene_df.reset_index(inplace=True)
        chrom_gene_df = chrom_gene_df[['chr', 'gene']]

        # append chromosome-gene index to DI scores and save.
        rho_df = DataFrame(self.rho
                           , columns=sample_ids)
        rho_df = concat([chrom_gene_df, rho_df]
                        , axis=1)
        rho_df = rho_df[['chr', 'gene'] + sample_ids]
        rho_df.to_csv(os.path.join(output_dir, 'degradation_index_scores.csv')
                      , index=False)

        # append chromosome-gene index to adjusted read counts and save.
        x_adj_df = DataFrame(self.x_adj
                             , columns=sample_ids)
        x_adj_df = concat([chrom_gene_df, x_adj_df]
                          , axis=1)
        x_adj_df = x_adj_df[['chr', 'gene'] + sample_ids]
        x_adj_df.to_csv(os.path.join(output_dir, 'adjusted_read_counts.csv')
                        , index=False)

        # append chromosome-gene index to ran-through-baseline-selection tracker and save.
        iter_names = ['iter_{0}'.format(i) for i in range(self.degnorm_iter)]
        ran_baseline_selection_df = DataFrame(self.ran_baseline_selection
                                              , columns=iter_names)
        ran_baseline_selection_df = concat([chrom_gene_df, ran_baseline_selection_df]
                                           , axis=1)
        ran_baseline_selection_df = ran_baseline_selection_df[['chr', 'gene'] + iter_names]
        ran_baseline_selection_df.to_csv(os.path.join(output_dir, 'ran_baseline_selection.csv')
                                         , index=False)
