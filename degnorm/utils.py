import multiprocessing as mp


def subset_to_chrom(df, chrom):
    return df[df['chr'] == chrom]


def max_cpu():
    return mp.cpu_count() - 1