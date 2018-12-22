# Running DegNorm

DegNorm is a command line interface (CLI) tool. You can access it through the `degnorm` command and the MPI-enabled
distributed version is `degnorm_mpi`.

-----

 
## Inputs
You only need two types of files to supply `degnorm` - .bam files and a .gtf file
if you have `samtools` in your `$PATH`. If you do not have `samtools`, you will also need [bam index files]([bam index files](https://www.biostars.org/p/15847/)) to accompany your .bam files.

#### 1. Align and sort reads data
Before running DegNorm, you will need to sort (and optionally, index) them with `samtools sort`. This command just re-orders reads
by chromosome and starting index, so that they can be more useful to DegNorm later on.

To [sort](https://davetang.org/wiki/tiki-index.php?page=SAMTools#Sorting_a_BAM_file) an alignment file, suppose it's called `S1.bam`, you can sort it with the following command:

    samtools sort S1.bam -o S1_sorted.bam
        
DegNorm requires sorted .bam files because it really needs [bam index files](https://www.biostars.org/p/15847/). If you don't create them
first, DegNorm will run `samtools index` to create an index file for each input alignment file:

    samtools index S1_sorted.bam S1_sorted.bai


#### 2. (Paired or single end) aligned and sorted reads data
At least two alignment files (.bam) must be specified, as inter-sample degradation normalization cannot happen on a standalone 
RNA-Seq sample. Use `p` to refer to the total number of samples.

If .bai files are not submitted, `degnorm` will look for .bai files named after the .bam files only with the ".bai" extension.
If no such file is found, `degnorm` will attempt to build one with the `samtools index` command. This will only work if the .bam files are sorted.
Instead of specifying individual .bam and .bai files, you can just specify `--bam-dir`, a path to a directory holding the relevant .bam and .bai files.
  With `--bam-dir`, it is assumed that the .bai files are named according to the .bam files, only with the ".bai" extension.

**You can supply paired reads files or single end reads files**.

Argument    | Required? |    Meaning
----------- | --------- | ------------
`--bam-files` | If neither `--warm-start-dir` nor `--bam-dir` are specified | Set of individual .bam files
`--bai-files` | If `samtools` is not in the `$PATH` and neither `--warm-start-dir` nor `--bam-dir` are specified | Set of individual .bai files. If specified, they must be in the order corresponding to `--bam-files`.
`--bam-dir`   | If neither `--warm-start-dir` nor `--bam-files` are specified | Directory containing .bam and .bai files for a pipeline run. It is assumed the .bai files have the same name as the .bam files, just with a different extension.
`-u`, `--unique-alignments` | optional flag | If specified, tells `degnorm` to remove reads aligned to more than one location on the genome. Suggested for use with single end reads data.

#### 3. Genome annotation file
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


### Sorting, indexing alignment files

Remember - we have to sort and index the alignment files first, and then we can run `degnorm`. Suppose we have set of alignment
files in our current working directory that we would like to send through DegNorm:

    # sort and index alignment files in current working directory.
    for FILE in ./*.bam
    do
        echo 'sorting '$FILE
        samtools sort $FILE -o ${FILE/.bam/}'_sorted.bam'
        echo 'indexing '${FILE/.bam/}'_sorted.bam'
        samtools index ${FILE/.bam/}'_sorted.bam' ${FILE/.bam/.bai}
    done

### Launching a full `degnorm` run "the verbose way"

Next, after the alignment files have been sorted (and optionally indexed), let's run `degnorm` "the verbose way" by enumerating 
each .bam and .bai file individually. Suppose we have two sorted alignment files, `S1_sorted.bam` and `S2_sorted.bam`.
    
    # run degnorm on two samples on 20 cores, using 100 NMF iterations per gene
    degnorm --bam-files S1_sorted.bam S2_sorted.bam \
        --bai-files S1_sorted.bai S2_sorted.bai \
        -g human.gtf \
        -p 20 \
        -o degnorm_output

If we did not specify any `--bai-files`, then DegNorm would automatically search for files "S1_sorted.bai" and "S2_sorted.bai" because they're named after
the alignment files, "S1_sorted.bam" and "S2_sorted.bam." If those index files cannot be found, they will be created in the same directory
  containing the input .bam files, only if `samtools` is in your `$PATH`.

### Launching a full `degnorm` run "the easy way"

It's probably just easier to supply `degnorm` with an input directory containing all of our .bam files, instead of each .bam file individually.
Suppose that data directory is located at `degnorm_data/GBM`.
 Let's also use 5 DegNorm iterations, 50 NMF iterations per gene per NMF iteration, and route output to a directory besides the current working directory.

    degnorm --bam-dir degnorm_data/GBM \
        -g human.gtf \
        -p 20 \
        --nmf-iter 50 \
        -o degnorm_output

### Using a warm start directory

After one pipeline run, we could start `degnorm` from the output directory resulting from the first run - this directory is known as the "warm start" directory. This time, let's exclude genes with a maximum coverage
(across all samples) less than 20, and run DegNorm with only 3 iterations. In this example, `DegNorm_102118_081045` was generated automatically from one of the prior DegNorm runs.

    degnorm --warm-start-dir degnorm_output/DegNorm_102118_081045 \
        -p 20 \
        -o degnorm_output \
        --minimax-coverage 20 \
        --iter 3
  
### Launching distributed runs with MPI

Let's speed this up by distributing DegNorm's workload over multiple servers with `degnorm_mpi`!
 If the `mpi4py` python package has been installed (meaning that MPICH is available in your compute environment).
 You can still specify the number of cores to use on each node with `-p`.
 
    # run degnorm_mpi using 4 nodes ("-n 4")
    mpiexec -n 4 degnorm_mpi \
        --warm-start-dir degnorm_output/DegNorm_GBM_102018 \
        -p 20 \
        -o degnorm_output \
        --minimax-coverage 20
        --nmf-iter 100
        
 The functionality of `degnorm_mpi` is the same as the single-node `degnorm`, so you can continue using warm start directories just like you could with `degnorm`.
 If you are unfamiliar with MPI, a great 
        
## Output

Raw and estimated coverage data are stored in separate `.pkl` files, one file per chromosome. See the [posthoc analysis](../howtos/posthoc_analysis.md)
 documentation for helper functions to access coverage data.

In addition to per-gene coverage matrices, `degnorm` will produce raw and degradation-adjusted read counts, a matrix of degradation index scores,
a matrix describing which genes were sent through the baseline selection procedure and when, along with various summary graphics and a pipeline summary report. 
The report will render as an .html file, but if you have `pandoc` installed
it will be converted to a .pdf.


    ├── degnorm.log
    ├── degradation_index_scores.csv
    ├── ran_baseline_selection.csv  # matrix of Booleans: which genes go thru baseline selection
    ├── gene_exon_metadata.csv  # chromosome, gene, exon relationship data
    ├── read_counts.csv  # raw read counts
    ├── adjusted_read_counts.csv  # degradation normalized read counts
    ├── <SAMPLE_1>
    │   ├── sample_SAMPLE_1_chr1.npz  # full chrom coverage for sample 1 in compressed matrix.
    │   ├── sample_SAMPLE_1_chr2.npz
    │   ├── <more sample chromosome coverage .npz files>
    ├── <SAMPLE_2>
    │   ├── sample_SAMPLE_2_chr1.npz
    │   ├── sample_SAMPLE_2_chr2.npz
    │   ├── <more per-sample chromosome coverage .npz files>
    ├── <more sample directories with chromosome coverage>
    ├── chr1
    │   ├── GAPDH_coverage.png
    │   ├── <more coverage plots>
    │   ├── coverage_matrices_chr1.pkl  # (gene, raw coverage matrix) pairs in serialized python object.
    │   └── estimated_coverage_matrices_chr1.pkl # (gene, estimated coverage matrix) pairs in serialized python object.
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