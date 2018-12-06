# Running DegNorm

DegNorm is a command line interface (CLI) tool. You can access it through the `degnorm` command and the MPI-enabled
distributed version is `degnorm_mpi`.
 
## Inputs
You only need two types of files to supply `degnorm`: .bam (and their corresponding .bai files - [bam index files](https://www.biostars.org/p/15847/)), and a .gtf file.
If you have `samtools` in your `$PATH`, .bai files can be created for you automatically, if they are not supplied.

#### 1. (Paired or single end) aligned and sorted reads data
At least two alignment files (.bam) must be specified, as inter-sample degradation normalization cannot happen on a standalone 
RNA-Seq sample. Use `p` to refer to the total number of samples.

It is assumed your .bam files are **sorted** by leftmost start coordinate (i.e. with `samtools sort`), contain a header, and abide by the [conventions](http://samtools.sourceforge.net/SAM1.pdf). If .bai files are not submitted,
`degnorm` will look for .bai files named after the .bam files only with the ".bai" extension. If no such file is found, `degnorm` will attempt to build one with the `samtools index` command. This will only work if the .bam files are sorted.
Instead of specifying individual .bam and .bai files, you can just specify `--bam-dir`, a path to a directory holding the relevant .bam and .bai files.
  With `--bam-dir`, it is assumed that the .bai files are named according to the .bam files, only with the ".bai" extension.

**You can supply paired reads files or single end reads files**.

Argument    | Required? |    Meaning
----------- | --------- | ------------
`--bam-files` | If neither `--warm-start-dir` nor `--bam-dir` are specified | Set of individual .bam files
`--bai-files` | If `samtools` is not in the `$PATH` and neither `--warm-start-dir` nor `--bam-dir` are specified | Set of individual .bai files. If specified, they must be in the order corresponding to `--bam-files`.
`--bam-dir`   | If neither `--warm-start-dir` nor `--bam-files` are specified | Directory containing .bam and .bai files for a pipeline run. It is assumed the .bai files have the same name as the .bam files, just with a different extension.
`-u`, `--unique-alignments` | optional flag | If specified, tells `degnorm` to remove reads aligned to more than one location on the genome. Suggested for use with single end reads data.

#### 2. Genome annotation file
DegNorm needs a .gtf file in order to construct the total transcript and compute all per-gene coverage score curves. It is assumed your .gtf file abides by the standard [conventions](https://useast.ensembl.org/info/website/upload/gff.html).

Argument    | Required? |    Meaning
----------- | --------- | ------------
`-g`, `--genome-annotation` | Yes | .gtf file for relevant genome.


## Using a warm start directory
Loading multiple .bam files, parsing a .gtf file, and computing per-gene cross-sample coverage matrices and read counts, can take some time, but this is a one-time
 preprocessing cost given a specific set of RNA-Seq samples. If you run the `degnorm` pipeline once, you can leverage
the stored coverage, read counts, and parsed transcript data from an existing DegNorm output directory to start a a new DegNorm run where all of the
preprocessing is completed ahead of time. When using `--warm-start-dir` you do *not* need to supply a .gtf file.

Argument    | Required? |    Meaning
----------- | --------- | ------------
`--warm-start-dir` | No | A directory holding the contents of a successful previous `degnorm` pipeline run. When specified, there is no need to supply a .gtf file.

## Runtime and algorithmic parameters

Argument    | Required? |    Meaning
----------- | --------- | ------------
`-o`, `--output-dir` | No | Defaults to the current working directory. Use to specify location where pipeline output directory will be written.
`--plot-genes` | No | Names of genes for which to render coverage plots. Sequence of explictly stated gene names or a .txt file containing one gene name per line.
`-d`, `--downsample-rate` | No | EXPERIMENTAL. Integer downsampling rate. Systematic samples of a coverage matrix are used to speed up NMF iterations.
`--nmf-iter` | No | Number of iterations per NMF-OA approximation. The higher the more accurate the approximation, but the more costly in terms of time.
`--iter` | No | Number of whole DegNorm iterations. Default is 5.
`--minimax-coverage` | No | Minimum cross-sample maximum coverage for a gene before it is included in the DegNorm pipeline. Can be used to exclude relatively low-coverage genes.
 `-s`, `--skip-baseline-selection` | No | EXPERIMENTAL. Flag to skip baseline selection, will greatly speed up DegNorm iterations.
 `-p`, `--proc-per-node` | No | Integer number of processes to spawn per compute node. The more the better.


## Example usage

Run `degnorm` "the hard way" by enumerating each .bam file individually. Suppose we have two alignment files, `degnorm_data/GBM/S1.bam` and `degnorm_data/GBM/S2.bam`.
 If we do not specify any .bai files, `degnorm` will search for "degnorm_data/GBM/S1.bai" and  "degnorm_data/GBM/S2.bai." If they cannot be found, they will be created within `degnorm_data/GBM` if `samtools` is in your `$PATH`.

    degnorm --bam-files degnorm_data/GBM/S1.bam degnorm_data/GBM/S2.bam \
        -g human.gtf \
        -p 20 \
        --nmf-iter 50 \
        -o degnorm_output

It's probably just easier to supply `degnorm` with an input directory (containing multiple .bam files) with 20 threads, use 5 DegNorm iterations and 50 NMF iterations per gene per NMF iteration.
Route output to a directory besides the current working directory.

    degnorm --bam-dir degnorm_data/GBM \
        -g human.gtf \
        -p 20 \
        --nmf-iter 50 \
        -o degnorm_output

After one pipeline run, we could start `degnorm` from the output directory resulting from the first run - this directory is known as the "warm start" directory. This time, let's exclude genes with a maximum coverage
(across all samples) less than 20, and run DegNorm with only 3 iterations. In this example, `DegNorm_102118_081045` was generated automatically from one of the prior DegNorm runs.

    degnorm --warm-start-dir degnorm_output/DegNorm_102118_081045 \
        -p 20 \
        -o degnorm_output \
        --minimax-coverage 20 \
        --iter 3
        
 If the `mpi4py` python package has been installed (meaning that MPICH is available in your compute environment), try distributing the DegNorm
 workload over multiple servers (also often referred to as "nodes") with `degnorm_mpi`. You can still specify the number of threads to use within each node with `-p`.
 
    mpiexec -n 4 degnorm_mpi \
        --warm-start-dir degnorm_output/DegNorm_GBM_102018 \
        -p 20 \
        -o degnorm_output \
        --minimax-coverage 20
        --nmf-iter 100
        
 The functionality of `degnorm_mpi` is the same as the single-node `degnorm`, so you can continue using warm start directories just like you could with `degnorm`.
        
## Output

Raw and estimated coverage data are stored in separate `.pkl` files, one file per chromosome. See the [posthoc analysis](../howtos/posthoc_analysis.md)
 documentation for helper functions to access coverage data.

In addition to per-gene coverage matrices, `degnorm` will produce raw and degradation-adjusted read counts, a matrix of degradation index scores,
a matrix describing which genes were sent through the baseline selection procedure and when, along with various summary graphics and a pipeline summary report. 
The report will render as an .html file, but if you have `pandoc` installed
it will be converted to a .pdf.


    ├── degnorm.log
    ├── degradation_index_scores.csv
    ├── ran_baseline_selection.csv
    ├── adjusted_read_counts.csv
    ├── read_counts.csv
    ├── adjusted_read_counts.csv
    ├── SAMPLE_1
    │   ├── sample_SAMPLE_1_chr1.npz
    │   ├── sample_SAMPLE_1_chr2.npz
    │   ├── <more sample chromosome coverage .npz files>
    ├── SAMPLE_2
    │   ├── sample_SAMPLE_2_chr1.npz
    │   ├── sample_SAMPLE_2_chr2.npz
    │   ├── <more per-sample chromosome coverage .npz files>
    ├── <more sample directories with chromosome coverage>
    ├── chr1
    │   ├── GAPDH_coverage.png
    │   ├── <more coverage plots>
    │   ├── coverage_matrices_chr1.pkl
    │   └── estimated_coverage_matrices_chr1.pkl
    ├── chr2
    │   ├── coverage_matrices_chr2.pkl
    │   └── estimated_coverage_matrices_chr2.pkl
    ├── chr3
    │   ├── coverage_matrices_chr3.pkl
    │   └── estimated_coverage_matrices_chr3.pkl
    ├── chr4
    │   ├── GAPDH_coverage.png
    │   ├── coverage_matrices_chr4.pkl
    │   └── estimated_coverage_matrices_chr4.pkl
    ├── <more chromosome directories>
    │   ├── <more raw coverage matrix .pkl files>
    │   └── <more estimated coverage matrix .pkl files>
    └── report
        ├── degnorm_summary.pdf # (or .html if pandoc not available)
        ├── di_boxplots.png
        ├── di_correlation.png
        └── di_heatmap.png