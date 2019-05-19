from scipy.sparse.linalg import svds
from pandas import DataFrame, concat
from degnorm.utils import *
import warnings
from collections import OrderedDict
import pickle as pkl
from joblib import Parallel, delayed


def rank_one_approx(x):
    """
    Decompose a matrix X via truncated SVD into (K)(E^t) = U_{1} \cdot \sigma_{1}V_{1}

    :param x: numpy 2-d array
    :return: 2-tuple (K, E) matrix factorization
    """
    u, s, v = svds(x, k=1)
    return u * s, v


def get_high_coverage_idx(x):
    """
    Find positions of high coverage in a gene's coverage matrix, defined
    as positions where the sample-wise maximum coverage is at least 10%
    the maximum of the entire coverage array.

    :param x: numpy 2-d array: coverage matrix (p x Li)
    :return: numpy 1-d array of integer base positions
    """
    return np.where(x.max(axis=0) > 0.1 * x.max())[0]


def nmf(x, factors=False, nmf_iter=100):
    """
    Run NMF-OA approximation. See "Normalization of generalized transcript degradation
    improves accuracy in RNA-seq analysis" supplement section 1.2.

    :param x: numpy 2-d array
    :param factors: boolean return 2-tuple of K, E matrix factorization? If False,
    return K.dot(E)
    :param nmf_iter: int number of SVD iterations per NMF approximation.
    :return: depending on factors, return (K, E) matrices or K.dot(E) over-approximation to x
    """
    K, E = rank_one_approx(x)
    est = K.dot(E)
    lmbda = np.zeros(shape=x.shape)
    c = 1. / np.sqrt(nmf_iter)

    for _ in range(nmf_iter):
        res = est - x
        lmbda -= c * res
        lmbda[lmbda < 0.] = 0.
        K, E = rank_one_approx(x + lmbda)
        est = K.dot(E)

    if factors:
        return K, E

    return est


def ratio_svd(x):
    """
    One-iteration SVD over-approximation, but not the NMFOA algorithm.
    See https://bit.ly/2zR4XEn.

    :param x: 2-d numpy array
    :return: 2-d numpy array estimate, elements are at least as large as those in x.
    """
    K, E = rank_one_approx(x)
    est = K.dot(E)
    est[est < x] = x[est < x]

    return est


def run_ratio_svd_serial(x):
    return list(map(ratio_svd, x))


def par_apply(fun, dat, n_jobs, mem_splits):
    """
    Apply a function to a list in parallel with threading.

    :param fun: function to be applied to elements of list.
    :param dat: list of elements to be chunked up, distributed to threads, and hit with fun.
    :param n_jobs: number of threads for within-node parallelization of work.
    :param mem_splits: number of splits to break data into. A split is sent to one thread for processing.
    :return: list of output from fun
    """
    # split up list of coverage matrices to feed to workers.
    dat = split_into_chunks(dat, n=mem_splits)
    par_output = Parallel(n_jobs=min(n_jobs, len(dat))
                          , verbose=0
                          , backend='threading')(map(delayed(fun), dat))

    return flatten_2d(par_output
                      , arr=False)


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


def adjust_coverage_curves(dat, scale_factors):
    """
    Adjust coverage curves by dividing the coverage values by corresponding scale factor, s_{j}
    """
    return [(F.T / scale_factors).T for F in dat]


def correct_di_scores(rho, x_weighted, x_adj):
    """
    Correct degradation index (DI) scores, i.e. the rho matrix, by assigning the sample average DI
    score to genes that were not run through baseline selection.
    """

    # determine indices of genes that were not sent through baseline selection; if they exist, do the swap.
    non_baseline_gene = rho.max(axis=1) == 0

    if np.sum(non_baseline_gene) > 0:
        sample_avg_di_scores = 1 - (x_weighted.sum(axis=0) / x_adj.sum(axis=0))
        rho[non_baseline_gene, :] = sample_avg_di_scores

    return rho


def systematic_sample(n, take_every=1):
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


def run_baseline_selection_serial(x, **kwargs):
    return list(map(lambda z: baseline_selection(z, **kwargs), x))


def baseline_selection(F, nmf_iter=100, downsample_rate=1, min_high_coverage=20,
                       bins=20, bin_frac=0.2, skip_baseline_selection=False):
    """
    Find "baseline" region for a gene's coverage curves - a region where it is
    suspected that degradation is minimal, so that the coverage envelope function
    may be estimated more accurately, in addition to finding the gene's DI scores
    across the samples.

    This is typically applied once coverage curves have been scaled by 1 / s_{j}.

    :param F: numpy 2-d array, gene coverage curve matrix, a numpy 2-dimensional array: coverage curve matrix,
     scaled by the degradation normalization scale factor. Wide-format (shape p x Li).
    :param nmf_iter: int number of SVD iterations per NMF approximation.
    :param downsample_rate: int systematic sampling take every rate. See run_nmfoa* functions.
    :param min_high_coverage: int minimum number of high-coverage base positions for F
    to be sent through baseline selection.
    :param bins: int number of bins to use for breaking up transcript into distinct coverage regions.
    :param bin_frac: float in (0, 1], fraction of bins required to be kept before exiting baseline selection procedure.
    :param skip_baseline_selection: Bool skip baseline selection?

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
    p = F.shape[0]

    output = dict()
    output['rho'] = np.zeros(p)  # default degradation scores are all 0.
    output['estimate'] = F
    output['ran_baseline_selection'] = False

    hi_cov_idx = get_high_coverage_idx(F)

    # Run systematic downsample if desired, prior to filtering out low-coverage regions.
    if downsample_rate > 1:
        _, downsample_idx = downsample_2d(F
                                          , downsample_rate=downsample_rate
                                          , by_row=False)

        # intersect sampled indices with the high-coverage indices.
        hi_cov_idx = np.intersect1d(downsample_idx, hi_cov_idx)

    n_hi_cov = len(hi_cov_idx)

    # if there are not sufficiently-many high-coverage base pairs, return null degradation.
    if n_hi_cov < min_high_coverage:
        return output['estimate'], output['rho'], output['ran_baseline_selection']

    # select high-coverage positions from (possibly downsampled) coverage matrix.
    hi_cov_idx.sort()
    F_start = np.copy(F[:, hi_cov_idx])
    F_bin = np.copy(F_start)

    # if any sample sans coverage after filtering, return defaults.
    if np.sum(F_bin.sum(axis=1) > 0) < p:
        return output['estimate'], output['rho'], output['ran_baseline_selection']

    # run NMF on filtered coverage, obtain initial coverage curve estimate.
    K, E = nmf(F_bin
               , factors=True
               , nmf_iter=nmf_iter)
    KE_bin = K.dot(E)

    # keep original NMFOA-estimated coverage in case we do not run baseline selection.
    K_start, E_start = np.copy(K), np.copy(E)
    estimate = np.copy(KE_bin)

    # compute DI scores. High score ->> high degradation.
    rho_vec = 1 - F_bin.sum(axis=1) / (KE_bin.sum(axis=1) + 1)

    # exclude extreme cases where NMF result doesn't converge.
    if np.nanmedian(1 - rho_vec) > 1:
        return output['estimate'], output['rho'], output['ran_baseline_selection']

    # establish minimum length of gene to work with. Bin relic.
    min_gene_len = max(2, np.ceil(200.0 * (1 / downsample_rate)))  # scale 200 default by downsample rate.

    # minimum number of binned regions, below which exit baseline selection procedure.
    min_bins = np.ceil(bins * bin_frac)

    # decide whether to search for a gene's baseline regions.
    # If any of these criteria are satisfied, do not run baseline selection and return what we have.
    if (n_hi_cov >= min_gene_len) and (np.nanmin(rho_vec) <= 0.2) and (not skip_baseline_selection):

        # split up consecutive regions of the high-coverage gene.
        # even in downsampling regime, drop batches of sample points on each baseline selection iteration.
        bin_segs = split_into_chunks(list(range(F_bin.shape[1]))
                                     , n=bins)
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
            bin_segs = shift_bins(bin_segs
                                  , dropped_bin=drop_idx)
            n_bins = len(bin_segs)

            # shrink F matrix to the indices not dropped, cast back to wide matrix.
            # estimate coverage curves from said indices with NMFOA. Try except in case
            # baseline selection has left us with 1-column coverage matrix.
            try:
                K, E = nmf(F_bin
                           , factors=True
                           , nmf_iter=nmf_iter)
            except ValueError:
                break

            KE_bin = K.dot(E)

            # stop is fitted values are all zero for any sample (extreme cases).
            if np.min(KE_bin.sum(axis=1)) == 0:
                break

            KE_bin[KE_bin < F_bin] = F_bin[KE_bin < F_bin]

            # recompute DI scores: closer to 1 ->> more degradation. Then run quality control.
            rho_vec = 1 - F_bin.sum(axis=1) / (KE_bin.sum(axis=1) + 1)  # + 1 as per Bin's code.

            if (n_bins <= min_bins) or (n_hi_cov < min_gene_len):
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
    return output['estimate'], output['rho'], output['ran_baseline_selection']


def par_apply_baseline_selection(dat, n_jobs, mem_splits, **kwargs):
    """
    Apply baseline selection algorithm to all gene coverage matrices in parallel,
    building the DI score matrix (rho) in the process and
    returning the estimated coverage matrices for visualization purposes.

    :param dat: list of numpy reads coverage arrays, shapes are (Li x p). Ideally
    coverage arrays will have been already adjusted for sequencing-depth by scale factors.
    :param n_jobs: int number of threads to use for running baseline selection in parallel on a single node.
    :param mem_splits: int number of splits to apply to dat. Threads work on individual splits, so
    more splits means less work per worker, but more times each worker needs to get a new split to work on.
    :return: 3-tuple --
    (list of 2-d numpy arrays used to visualize DegNorm estimated coverage matrices
    , rho (DI score matrix)
    , baseline selection Boolean indicator vector, one T/F per gene)
    """
    # split up coverage matrices so that no worker gets much more than 50Mb.
    dat = split_into_chunks(dat
                            , n=mem_splits)
    baseline_dat = Parallel(n_jobs=min(n_jobs, len(dat))
                            , verbose=0
                            , backend='threading')(delayed(run_baseline_selection_serial)(
        x=x
        , **kwargs) for x in dat)

    # flatten results.
    baseline_dat = [est for est1d in baseline_dat for est in est1d]

    # Assemble DI score matrix, apply QA thresholding on output.
    rho = np.vstack([x[1] for x in baseline_dat])
    rho[rho > 0.9] = 0.9
    rho[rho < 0.] = 0.

    # return coverage curve estimates (for visualization purposes, not for DI score calculations),
    # boolean array indicating whether or not each gene was sent through baseline selection
    return [x[0] for x in baseline_dat], rho, np.array([x[2] for x in baseline_dat])


def downsample_2d(x, downsample_rate=1, by_row=True):
    """
    Downsample a coverage matrix at evenly spaced base position indices via systematic sampling.

    :param x: 2-d numpy array; a coverage matrix
    :param downsample_rate: int systematic sampling rate. See systematic_sample.
    :param by_row: bool indicator - downsample by rows or columns?
    :return: 2-d numpy array
    """
    Li = x.shape[0 if by_row else 1]

    # if downsample rate high enough to cover entire matrix,
    # return input matrix and consider everything sampled.
    if downsample_rate == 1:
        return x, np.arange(0, Li)

    if downsample_rate >= Li:
        raise ValueError('Cannot downsample at a rate < 1 / length(gene)')

    # obtain systematic sample of transcript positions.
    downsample_idx = systematic_sample(Li
                                       , take_every=downsample_rate)

    if by_row:
        return x[downsample_idx, :], downsample_idx

    return x[:, downsample_idx], downsample_idx


def save_results(gene_manifest_df,
                 estimates, rho, x_adj, ran_baseline_selection,
                 sample_ids, output_dir='.'):
    """
    DegNorm output data to disk:
        - estimates: for each gene's estimated coverage matrix, find the chromosome to which the
        gene belongs. Per chromosome, assemble dictionary of {gene: estimated coverage matrix} pairs,
        and save that dictionary in a directory labeled by the chromosome in the output_dir directory,
        in a file "estimated_coverage_matrices_<chromosome name>.pkl"
        - rho: save degradation index scores to "degradation_index_scores.csv" using sample_ids
        as the header.
        - x_adj: save adjusted read counts to "adjusted_read_counts.csv" using sample_ids as the header.
        - ran_baseline_selection: save records of which genes were run through baseline selection on which
        DegNorm iteration to "ran_baseline_selection.csv". Columns are iterations of the DegNorm algorithm.

    :param gene_manifest_df: pandas.DataFrame establishing chromosome-gene map. Must contain,
    at a minimum, the columns `chr` (str, chromosome) and `gene` (str, gene name).
    :param estimates: dictionary of the form {gene: estimated coverage matrix}
    :param rho: 2-d numpy array of DI scores. Columns correspond to sample_ids.
    :param x_adj: 2-d numpy array of adjusted read counts. Columns correspond to sample_ids.
    :param ran_baseline_selection: 2-d numpy array Boolean indicators for which genes were run
    through baseline selection on which DegNorm iteration.
    :param output_dir: str output directory to save GeneNMFOA output.
    :param sample_ids: (optional) list of str names of RNA_seq experiments, to be used
     as headers for adjusted read counts matrix and degradation index score matrix.
    """
    genes = list(estimates.keys())

    if not os.path.isdir(output_dir):
        raise IOError('Directory {0} not found.'.format(output_dir))

    # ensure that gene_manifest_df has `chr` and `gene` columns to establish chromosome-gene map.
    if not all([col in gene_manifest_df.columns.tolist() for col in ['chr', 'gene']]):
        raise ValueError('gene_manifest_df must have columns `chr` and `gene`.')

    gene_intersect = np.intersect1d(gene_manifest_df.gene.unique(), genes)
    if len(gene_intersect) < len(genes):
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

        chrom_gene_dict[chrom][gene] = estimates.get(gene)

    # save estimated coverage matrices to genes nested within chromosomes.
    chrom_gene_dfs = list()
    for chrom in manifest_chroms:
        chrom_dir = os.path.join(output_dir, str(chrom))

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
    chrom_gene_df = chrom_gene_df.loc[genes]
    chrom_gene_df.reset_index(inplace=True)
    chrom_gene_df = chrom_gene_df[['chr', 'gene']]

    # append chromosome-gene index to DI scores and save.
    rho_df = DataFrame(rho
                       , columns=sample_ids)
    rho_df = concat([chrom_gene_df, rho_df]
                    , axis=1)
    rho_df = rho_df[['chr', 'gene'] + sample_ids]
    rho_df.to_csv(os.path.join(output_dir, 'degradation_index_scores.csv')
                  , index=False)

    # append chromosome-gene index to adjusted read counts and save.
    x_adj_df = DataFrame(x_adj
                         , columns=sample_ids)
    x_adj_df = concat([chrom_gene_df, x_adj_df]
                      , axis=1)
    x_adj_df = x_adj_df[['chr', 'gene'] + sample_ids]
    x_adj_df.to_csv(os.path.join(output_dir, 'adjusted_read_counts.csv')
                    , index=False)

    # append chromosome-gene index to ran-through-baseline-selection tracker and save.
    iter_names = ['iter_{0}'.format(i) for i in range(ran_baseline_selection.shape[1])]
    ran_baseline_selection_df = DataFrame(ran_baseline_selection
                                          , columns=iter_names)
    ran_baseline_selection_df = concat([chrom_gene_df, ran_baseline_selection_df]
                                       , axis=1)
    ran_baseline_selection_df = ran_baseline_selection_df[['chr', 'gene'] + iter_names]
    ran_baseline_selection_df.to_csv(os.path.join(output_dir, 'ran_baseline_selection.csv')
                                     , index=False)


def run_gene_nmfoa_mpi(comm, cov_dat, reads_dat, degnorm_iter=5, downsample_rate=1, min_high_coverage=50,
                       nmf_iter=100, bins=20, n_jobs=1, skip_baseline_selection=False, random_state=123):
    """
    Run DegNorm degradation normalization pipeline: adjust read counts, compute degradation index scores,
    and compute normalized coverage curve estimates.

    :param comm: mpi4py communicator, group of intracommunicating processes.
    :param cov_dat: OrderedDict of {gene: coverage matrix} pairs;
    coverage matrix shapes are (p x Li) (wide matrices). For each coverage matrix,
    row index is sample number, column index is gene's relative base position on chromosome.
    :param reads_dat: n (genes) x p (samples) numpy array of gene read counts
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

    :return: list of 2-d numpy arrays, estimated coverage matrices. In same order as the keys (genes) of cov_dat.
    """
    size = comm.size
    rank = comm.rank

    # everyone needs DegNorm algorithm params. quickly qc them.
    degnorm_iter = np.abs(int(degnorm_iter))
    nmf_iter = np.abs(int(nmf_iter))
    bins = np.abs(int(bins))
    min_high_coverage = max(2, np.abs(int(min_high_coverage)))
    downsample_rate = np.abs(int(downsample_rate))
    x_adj = None
    skip_baseline_selection = skip_baseline_selection
    random_state = random_state

    # all coverage matrices must have >= 2 high-coverage indices if downsampling (svds function limitation).
    if downsample_rate > 1:
        min_high_coverage = 2

    n_genes = len(cov_dat)
    x = np.copy(reads_dat)

    # master: to distribute gene names, coverage matrices out to workers,
    # set up ran_baseline_selection tracker array.
    if rank == 0:
        all_genes = list(cov_dat.keys())
        genes_partition = split_into_chunks(all_genes
                                            , n=size)
        my_genes = genes_partition[0]

        # distribute {gene: coverage matrix} pairs.
        for worker_id in range(size):
            my_cov_dat = OrderedDict()

            for gene in genes_partition[worker_id]:
                my_cov_dat[gene] = cov_dat.get(gene)

            # display how many genes for which each worker is responsible.
            msg = '{node} will be responsible for {n} genes.'\
                .format(node='host' if worker_id == 0 else 'worker node ' + str(worker_id), n=len(my_cov_dat))
            logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

            # master keeps its own gene coverage matrix data.
            if worker_id == 0:
                master_cov_dat = my_cov_dat

            # master sends worker nodes their data.
            else:
                comm.send(my_cov_dat
                          , dest=worker_id
                          , tag=333 + worker_id)

        # initialize baseline selection tracker entirely False (no genes have gone through baseline selection yet).
        ran_baseline_selection = np.zeros(shape=[n_genes, degnorm_iter]).astype(bool)

        # master housekeeping.
        my_cov_dat = master_cov_dat
        p = my_cov_dat.get(my_genes[0]).shape[0]

        # ---------------------------------------------------------------------------- #
        # Master runs data quality checks:
        # 1. Check that read count matrix has same number of rows as number of coverage arrays.
        # 2. Check that all coverage arrays are 2-d matrices.
        # 3. Check that all coverage arrays are wider than they are tall.
        # 4. Check that downsample rate < length(gene) for all genes, if downsampling.
        # ---------------------------------------------------------------------------- #
        if x.shape[0] != n_genes:
            raise ValueError('Number of genes in read count matrix not equal to number of coverage matrices!')

        if not all(map(lambda z: z.ndim == 2, list(cov_dat.values()))):
            raise ValueError('Not all coverage matrices are 2-d arrays!')

        li_vec = np.array(list(map(lambda z: z.shape[1], list(cov_dat.values()))))
        if np.sum(li_vec / p < 1) > 0:
            logging.warning('At least one coverage matrix is taller than it is wide.'
                            'Ensure that coverage matrices are shaped (p x L_i).')

        if downsample_rate > 1:
            if not np.min(li_vec) >= downsample_rate:
                raise ValueError('downsample_rate is too large; take-every size > at least one gene.')

    # workers pick up their gene coverage matrix data.
    else:
        my_cov_dat = comm.recv(source=0
                               , tag=333 + rank)
        my_genes = list(my_cov_dat.keys())

    # determine (integer) number of data splits for threaded workers (50Mb per worker).
    mem_splits = int(np.ceil(np.sum(list(map(lambda z: z.nbytes, list(my_cov_dat.values())))) / 5e7))
    mem_splits = max(mem_splits, n_jobs)

    # ---------------------------------------------------------------------------- #
    # DegNorm initialization:
    # 1. Estimate initial scale factors from single-iteration rank-1 SVD estimates,
    #    obtain initial DI scores.
    # 2. Use initial DI scores to compute initial normalization factors.
    # 3. Normalize read counts by initial normalization factors.
    # 4. Initialize coverage scale factors with initial normalization factors.
    # ---------------------------------------------------------------------------- #

    # an individual worker has their own gene matrices by now,
    # so use ratio_svd to obtain first coverage matrix estimates, use to compute initial DI scores.
    estimates = par_apply(fun=run_ratio_svd_serial
                          , dat=list(my_cov_dat.values())
                          , n_jobs=n_jobs
                          , mem_splits=mem_splits)
    estimates = {my_genes[i]: estimates[i] for i in range(len(my_genes))}

    # everyone sends their estimates to master. Getting memory/pickle issues with .gather,
    # so have individual workers send their estimates to master.
    if rank > 0:
        comm.send(estimates
                  , dest=0
                  , tag=444 + rank)

    # master constructs rho and normalizes read counts.
    if rank == 0:

        # master receives ratio_svd estimate dictionaries from workers.
        for worker_id in range(1, size):
            estimates.update(comm.recv(source=worker_id
                                       , tag=444 + worker_id))

        est_sums = np.vstack(list(map(lambda x: x.sum(axis=1), [estimates.get(gene) for gene in all_genes])))
        cov_sums = np.vstack(list(map(lambda x: x.sum(axis=1), list(cov_dat.values()))))
        rho = 1 - (cov_sums / (est_sums + 1))

        # estimate normalization factors from initial DI scores.
        low_di_gene = rho.max(axis=1) < 0.1
        count_sums = x[low_di_gene, :].sum(axis=0) if np.any(low_di_gene) else x.sum(axis=0)
        norm_factors = count_sums / np.median(count_sums)

        # adjust read counts by initial normalization factors, set initial scale factors.
        x_weighted = x / norm_factors
        scale_factors = np.copy(norm_factors)

        # display initialized scale factors.
        msg = 'Initial sequencing depth scale factors -- \n\t{0}' \
            .format(', '.join([str(x) for x in scale_factors]))
        logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

    # ---------------------------------------------------------------------------- #
    # DegNorm iterations:
    # 1. Scale coverage curves by scale factors from iteration t-1.
    # 2. Run baseline selection, obtain degradation-normalized over-approximated coverage curves.
    # 3. Compute DI scores from NMF-OA estimates, treating genes that went through
    #    baseline selection differently than those who did not.
    # 4. Re-normalize read counts based on DI scores at time t.
    # 5. Update sequencing depth scale factors based on re-normalized read counts.
    # ---------------------------------------------------------------------------- #

    # set random number generator seed and wait up.
    np.random.seed(random_state)
    comm.Barrier()

    # run DegNorm iterations.
    i = 0
    while True:

        # master scales baseline_cov coverage curves by 1 / (updated adjustment factors),
        # then distributes adjusted coverage matrices among workers.
        if rank == 0:
            msg = 'DegNorm iteration {0} -- master adjusting coverage matrices for scale factors.' \
                .format(i + 1)
            logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

            adj_cov_dat = adjust_coverage_curves(list(cov_dat.values()), scale_factors)
            adj_cov_dat = {all_genes[i]: adj_cov_dat[i] for i in range(n_genes)}

            # master distributes dicts of {gene: adjusted coverage matrix} pairs.
            for worker_id in range(1, size):
                my_adj_cov_dat = OrderedDict()

                for gene in genes_partition[worker_id]:
                    my_adj_cov_dat[gene] = adj_cov_dat.get(gene)

                # getting MPI timeouts upon sending adjusted coverage matrix dicts.
                # Just write data to disk for particular worker, have them read in later.
                if worker_id > 0:
                    comm.send(my_adj_cov_dat
                              , dest=worker_id
                              , tag=555 + worker_id)

            # get host's adjusted coverage matrices.
            my_adj_cov_dat = OrderedDict()
            for gene in genes_partition[0]:
                my_adj_cov_dat[gene] = adj_cov_dat.get(gene)

        # workers receive their adjusted coverage matrices.
        if rank > 0:
            my_adj_cov_dat = comm.recv(source=0
                                       , tag=555 + rank)

        # everyone runs NMF-OA + baseline selection on their genes.
        msg = 'DegNorm iteration {0} -- running baseline selection on {1} genes.' \
            .format(i + 1, len(my_adj_cov_dat))
        logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

        baseline_output = par_apply_baseline_selection(list(my_adj_cov_dat.values())
                                                       , n_jobs=n_jobs
                                                       , mem_splits=mem_splits
                                                       , nmf_iter=nmf_iter
                                                       , downsample_rate=downsample_rate
                                                       , min_high_coverage=min_high_coverage
                                                       , bins=bins
                                                       , bin_frac=0.2
                                                       , skip_baseline_selection=skip_baseline_selection)

        # declare number of genes sent through baseline selection on this iteration.
        if not skip_baseline_selection:
            msg = 'DegNorm iteration {0} -- {1} genes were sent through baseline selection.' \
                .format(i + 1, np.sum(baseline_output[2]))
            logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

        # everyone sends their estimates, DI scores, and baseline selection trackers to master. Do this
        # with a specific tag though, so we can assemble everything in the correct gene order, as we
        # know which worker has been assigned which genes.
        if rank > 0:
            comm.send(baseline_output
                      , dest=0
                      , tag=666 + rank)

        # begin host's work of organizing everyone's baseline selection data and re-scaling read counts.
        else:
            # obtain host's baseline selection output.
            estimates, rho, bs_bool = [baseline_output[i] for i in range(3)]

            # get each worker's baseline selection output data, in order.
            # construct correctly ordered, comprehensive set of coverage matrix estimates, DI scores,
            # and baseline selection tracking.
            for worker_id in range(1, size):
                baseline_output = comm.recv(source=worker_id
                                            , tag=666 + worker_id)

                estimates += baseline_output[0]
                rho = np.vstack([rho, baseline_output[1]])
                bs_bool = np.append(bs_bool, baseline_output[2])

            # update indicators whether gene was sent thru baseline selection on DegNorm iteration i
            ran_baseline_selection[:, i] = bs_bool

            # adjust (weighted) read counts.
            x_adj = x_weighted / (1 - rho)

            # update DI scores based on who went through baseline selection.
            rho = correct_di_scores(rho
                                    , x_weighted=x_weighted
                                    , x_adj=x_adj)

            # adjust norm-factor-weighted read counts by DI scores.
            x_adj = x_weighted / (1 - rho)

            # get new norm factors.
            norm_factors = x_adj.sum(axis=0) / np.median(x_adj.sum(axis=0))

            # update read counts by adjusting degradation effect into sequencing depth
            x_weighted = x_weighted / norm_factors

            # increment new scale_factors.
            scale_factors = scale_factors * norm_factors

            # display scale_factors.
            msg = 'DegNorm iteration {0} -- sequencing depth scale factors: \n\t{1}' \
                .format(i + 1, ', '.join([str(x) for x in scale_factors]))
            logging.info('({rank}/{size}) -- {msg}'.format(rank=rank + 1, size=size, msg=msg))

        # increment the current degnorm iteration for everyone.
        i += 1

        # everyone catches up before the next iteration.
        comm.Barrier()

        # if sufficient iterations have passed, have master return results.
        if i == degnorm_iter:

            if rank == 0:
                output = {'estimates': dict(zip(all_genes, estimates))
                          , 'rho': rho
                          , 'x_adj': x_adj
                          , 'ran_baseline_selection': ran_baseline_selection}

                return output

            else:
                return None
