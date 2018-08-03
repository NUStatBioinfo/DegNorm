import numpy as np
from scipy.sparse.linalg import svds
from degnorm.utils import *
from collections import OrderedDict

class GeneNMFOA():

    def __init__(self, nmf_iter=100, tol=1e-3,
                 loop_iter=5, bins=20, n_jobs=max_cpu()):
        self.nmf_iter = nmf_iter
        self.tol = tol
        self.loop_iter = loop_iter
        self.bins = bins
        self.min_bins = np.ceil(self.bins * 0.3)
        self.x = None
        self.p = None
        self.n_genes = None
        self.coverage_sums = None
        self.scale_factors = None

        if self.min_bins <= 0:
            raise ValueError('bins is not large enough for baseline selection algorithm. Try setting bins=20.')

    # TODO: look up whether I need to flip order of singular vectors.
    def rank_one_approx(self, x):
        """
        Decompose a matrix X via truncated SVD into (K)(E^t) = U_{1} \cdot \sigma_{1}V_{1}

        :param x: numpy 2-dimensional array
        :return: 2-tuple (K, E) matrix factorization
        """
        u, s, v = svds(x, k=1)
        return u[::-1], s*v[::-1]

    def nmf(self, x, factors=False):
        K, E = self.rank_one_approx(x)
        est = K.dot(E)
        lmbda = np.zeros(shape=x.shape)
        c = np.sqrt(self.nmf_iter)

        for _ in range(self.nmf_iter):
            res = x - est
            lmbda -= res * (1 / c)
            lmbda[lmbda < 0.] = 0.
            K, E = np.abs(self.rank_one_approx(est + lmbda))

        if factors:
            return {'K': K
                    , 'E': E}

        return K.dot(E)

    def par_apply_nmf(self, dat):
        """
        Apply NMF-OA to a list of coverage matrices.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p)
        :return: list of dict output from self.nmf, one per numpy.array in dat
        """
        p = mp.Pool(processes=self.n_jobs)
        nmf_ests = [p.apply_async(self.nmf
                                  , args=(F)) for F in dat]
        p.close()

        return [x.get() for x in nmf_ests]

    def compute_scale_factors(self, estimates):
        """
        Update rho, the n x p matrix of degradation index scores.

        Compute s_{j} = \frac{sum_{i=1}^{n}\tilde{X}_{ij}}{median_{j}(\sum_{i=1}^{n}\tilde{X}_{ij})}
        vector (p x 1) of read count adjustment factors.

        :param estimates: list of rank-one approximations of coverage curves
        """

        # find rho matrix: n x p matrix of degradation index scores
        est_sums = list(map(lambda x: x.sum(axis=0), estimates))
        rho = 1 - self.coverage_sums / np.vstack(est_sums)

        # scale read count matrix by di scores
        scaled_counts_mat = self.x / (1 - rho)
        scaled_counts_sums= scaled_counts_mat.sum(axis=0)
        self.scale_factors = scaled_counts_sums / np.median(scaled_counts_sums)

    def adjust_coverage_curves(self, dat):
        """
        Adjust coverage curves by dividing the coverage values by corresponding adjustment factor, s_{j}

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p)
        :return: list of adjusted (scaled) numpy reads coverage arrays, shapes are (Li x p)
        """
        return [(F.T / self.scale_factors).T for F in dat]

    def matrix_sst(self, x, left_idx, right_idx):
        """
        Compute sum of squared elements for a slice of a 2-d numpy array.

        :param x: numpy 2-dimensional array
        :param left_idx: left index, beginning of slice
        :param right_idx: right index, end of slice (not inclusive)
        :return: float sum of squared values in matrix slice
        """
        return (x[:, left_idx:right_idx]**2).sum()

    def baseline_selection(self, F, KE):
        """

        :param F: gene coverage curves, a numpy 2-dimensional array
        :param KE: gene coverage curve estimate, numpy 2-dimensional array
        :return:
        """
        Li = F.shape[1]

        # Decide where to bin up regions of the gene.
        bin_vec = np.linspace(0, Li, num=(self.bins + 1), endpoint=True, dtype=int)
        bin_segs = [[bin_vec[i], bin_vec[i + 1]] for i in range(self.bins)]

        KE_bin = np.copy(KE)
        F_bin = np.copy(F)

        rho_vec = 1 - F_bin.sum(axis=1) / KE_bin.sum(axis=1)

        # Run baseline selection algorithm while the degradation index score is still high
        # or while there are still non-baseline regions left to trim.
        while (len(bin_segs) >= self.min_bins) and (np.nanmax(rho_vec) > 0.1):

            # compute normalized residuals for each bin.
            ss_r = np.array(list(map(lambda z: self.matrix_sst(KE_bin - F_bin, z[0], z[1]), bin_segs)))
            ss_f = np.array(list(map(lambda z: self.matrix_sst(F_bin, z[0], z[1]), bin_segs)))
            res_norms = ss_r / ss_f

            # drop the bin corresponding to the bin with the maximum normalized residual.
            bin_segs = np.delete(bin_segs, res_norms.argmax())

            # collapse the remaining bins into an array of indices to keep.
            keep_idx = [np.arange(bin[0], bin[1]) for bin in bin_segs]
            keep_idx = np.unique(flatten_2d(keep_idx))
            F_bin = F_bin[:, keep_idx]

            # estimate coverage curves from remaining regions.
            KE_bin = self.nmf(F_bin)

            # obtain DI scores from remaining regions.
            rho_vec = 1 - F_bin.sum(axis=1) / KE_bin.sum(axis=1)

        # Run NMF on baseline regions.
        K, E = self.nmf(F_bin, factors=True)

        # Use refined estimate of K (p x 1 vector) to get new KE product.
        return {'estimate': E.dot((F_bin.T / K).max(axis=1))
                , 'rho': rho_vec}

    def par_apply_baseline_selection(self, dat):
        """
        Apply baseline selection algorithm to all gene coverage curve matrices in parallel.

        :param dat: list of numpy reads coverage arrays, shapes are (Li x p). Ideally
        coverage arrays will be adjusted for sequencing-depth.
        :return: list of numpy reads coverage arrays estimates
        """

        # Obtain initial coverage estimates.
        estimates = self.par_apply_nmf(dat)

        p = mp.Pool(processes=self.n_jobs)
        ests = [p.apply_async(self.baseline_selection
                                  , args=(dat[i], estimates[i])) for i in range(len(dat))]
        p.close()

        return [x.get() for x in ests]

    def fit(self, coverage_dat, reads_dat):
        """
        Initialize estimates for the DegNorm iteration loop.

        :param coverage_dat: list of numpy reads coverage arrays, shapes are (Li x p)
        :param reads_dat: n (genes) x p (samples) numpy array of gene read counts
        """
        self.x = reads_dat
        self.n_genes = len(coverage_dat)
        self.p = coverage_dat[0].shape[0]

        if not all(np.array(list(map(lambda z: z.shape[0], coverage_dat))) == self.p):
            raise ValueError('Not all coverage matrices contain the same number of samples!')

        if not self.x.shape[0] == self.n_genes:
            raise ValueError('Number of genes in read count matrix not equal to number of coverage matrices!')

        # assemble n x p matrix of sums over coverage arrays (sum coverage over positions, one sum per sample)
        self.coverage_sums = np.array(list(map(lambda z: z.sum(axis=1), coverage_dat)))

    def fit_transform(self, coverage_dat, reads_dat):
        self.fit(coverage_dat, reads_dat)

        # initialize NMF-OA estimates of coverage curves
        estimates = self.par_apply_nmf(coverage_dat)
        norms_t0 = np.array(list(map(np.linalg.norm, estimates)))

        # initialize output storage.
        gene_dict = OrderedDict

        # Run DegNorm iterations.
        i = 0
        while i < self.loop_iter:

            # obtain vector of read count adjustment factors
            self.compute_scale_factors(estimates)

            # scale coverage curves by 1 / adjustment factors
            coverage_adjusted = self.adjust_coverage_curves(coverage_dat)

            # run NMF-OA + baseline selection; obtain refined estimate of K*Et
            baseline_dat = self.par_apply_baseline_selection(coverage_adjusted)
            estimates = [baseline_dat[i]['estimate'] for i in len(baseline_dat)]

            # determine which genes coverage curve estimates have converged.
            norms_t1 = np.array(list(map(np.linalg.norm, estimates)))
            delta_vec = np.abs(norms_t1 - norms_t0) / norms_t0
            converged_idx = np.where(delta_vec <= self.tol)[0]
            non_converged_idx = np.where(delta_vec > self.tol)[0]

            # store the converged genes' data and gene index.
            if len(converged_idx) > 0:
                for idx in converged_idx:
                    gene_dict[idx] = baseline_dat[idx]
                    gene_dict[idx]['converged'] = True

            # if all genes have converged, yay exit and continue.
            if len(non_converged_idx) == 0:
                break

            # drop the converged genes and update the coverage sums to reflect dropped genes.
            else:
                coverage_dat = [coverage_dat[idx] for idx in non_converged_idx]
                self.coverage_sums = self.coverage_sums[non_converged_idx]

            i += 1

        # additionally, store all of the non-converged genes.
        if len(non_converged_idx) > 0:
            for idx in non_converged_idx:
                gene_dict[idx] = baseline_dat[idx]
                gene_dict[idx]['converged'] = False

        return OrderedDict(sorted(gene_dict.items(), key=lambda z: z[0]))