# Running DegNorm

DegNorm is a CLI tool. You can access it through the `degnorm` command.
 
## Inputs
You only need two types of files to supply `degnorm`

#### 1. (Paired or single eng) aligned reads data
At least two aligned reads files must be supplied, as inter-sample degradation normalization can't happen on a standalone 
RNA-Seq expermient. Use `p` to refer to the total number of experiments. If you have `samtools` installed, feel free
to provide .bam files - DegNorm will convert them to their text equivalent via `samtools view -h -o {.sam} {.bam}`.

It is assumed your .sam files abide by the [conventions](https://en.wikipedia.org/wiki/SAM_(file_format)).
**You can supply paired reads files or single end reads files**.

Argument    | Required? |    Meaning
----------- | --------- | ------------
`-i`, `--input-files` |   Only if `--input-dir` not supplied | Individual .sam or .bam files. Any .bam files will be converted to .sam if `samtools` is available.
`--input-dir` | Only if input files not supplied | A directory containing a set of .sam files. In this case, .bam files will be ignored.
`-u`, `--unique-alignments` | optional flag | If specified, tells DegNorm to remove reads aligned to more than one area of the genome. Suggested for use with single end reads data.

#### 2. Genome annotation file
DegNorm needs a .gtf file to determine the transcript for computing the per-gene coverage curves, which span a 
concatenation of the coding regions only (introns are ignored). It is assumed your .gtf file abides by the [conventions](https://useast.ensembl.org/info/website/upload/gff.html)

Argument    | Required? |    Meaning
----------- | --------- | ------------
`-g`, `--genome-annotation` | Yes | .gtf file for relevant genome.


## Using a "warm start" directory
Loading multiple .sam files, parsing a .gtf files, and computing per-gene cross-sample coverage matrices and read counts, can take some time, but there is really no need
to run this process more than once for a specific set of RNA-Seq experiments. If you run the `degnorm` pipeline once, you can leverage
the stored coverage, read counts, and parsed transcript data from an existing DegNorm output directory to start a a new DegNorm run where all of the
data preprocessing is already completed. When using `--warm-start-dir` you do *not* need to supply a .gtf file.

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
 `-c`, `--cpu` | No | Integer number of threads. The more the better.


## Example usage

Run DegNorm with an input directory (containing multiple .sam files) with 20 threads, use 5 DegNorm iterations and 50 NMF iterations per gene per NMF iteration.
Route output to a directory besides the current working directory.

    degnorm --input-dir degnorm_data/GBM \
        -g human.gtf \
        -c 20 \
        --nmf-iter 50 \
        -o degnorm_output

Run DegNorm on from a warm start directory (from a previous run). This time, do not include genes with a maximum coverage
(across all samples) less than 20, and only run 3 DegNorm iterations.

    degnorm --warm-start-dir degnorm_output/DegNorm_GBM_102018 \
        -c 20 \
        -o degnorm_output \
        -minimax-coverage \
        --iter 3
        
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