# Running DegNorm

DegNorm is a command line interface (CLI) tool. You can access it through the `degnorm` command and the MPI-enabled
distributed version is `degnorm_mpi`.

-----

 
## Inputs
You only need two types of files to supply `degnorm` - .bam files (0-indexed) and a .gtf file (1-indexed) if you have `samtools` in your `$PATH`. If you do not have `samtools`, you will also need [bam index files]([bam index files](https://www.biostars.org/p/15847/)) to accompany your .bam files.

### 1. (Optional) index files (.bai)
Before running DegNorm, you will need to sort (and optionally, index) them with `samtools sort`. This command just re-orders reads
by chromosome and starting index, so that they can be more useful to DegNorm later on.

To [sort](https://davetang.org/wiki/tiki-index.php?page=SAMTools#Sorting_a_BAM_file) an alignment file, suppose it's called `S1.bam`, you can sort it with the following command:

    # sort alignment files (necessary for indexing them)
    samtools sort S1.bam -o S1_sorted.bam
        
DegNorm requires sorted .bam files because it really needs [bam index files](https://www.biostars.org/p/15847/). If you don't create them
first, DegNorm will run `samtools index` to create an index file for each input alignment file:

    # create alignment index files
    samtools index S1_sorted.bam S1_sorted.bai

Within the .bam files **reads must be 0-indexed**.

### 2. Paired or single-end aligned and sorted reads (.bam)
At least two alignment files (.bam) must be specified, as inter-sample degradation normalization cannot happen on a standalone 
RNA-Seq sample. We refer to the total number of samples as `p`.

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

### 3. Genome annotation file (.gtf)
DegNorm needs a .gtf file in order to construct the total transcript and compute all per-gene coverage score curves. It is assumed that the positions in the .gtf file are 1-indexed and that the file abides by the standard [.gtf conventions](https://useast.ensembl.org/info/website/upload/gff.html).

Argument    | Required? |    Meaning
----------- | --------- | ------------
`-g`, `--genome-annotation` | Yes | .gtf file for relevant genome.

**Start and end positions with the .gtf file must be 1-indexed**. Additionally, the .gtf file must have the 9 standard .gtf fields. Here is an example of a .gtf file we use to test `degnorm`:

seqname | source | feature | start | end | score | strand | frame | attribute
------- | ------ | ------- | ----- | --- | ----- | ------ | ----- | ----------
chr1 |	unknown	| exon |	323892 | 324060	| . |	+	| . |	gene_id "LOC100133331"; gene_name "LOC100133331"; transcript_id "NR_028327_1"; tss_id "TSS14025";
chr1 | unknown	| exon | 324288 | 324345 | . | + | . | gene_id "LOC100133331"; gene_name "LOC100133331"; transcript_id "NR_028327_1"; tss_id "TSS14025";
chr1 |unknown | exon | 324439	| 326938 |	. |	+ |	. |	gene_id "LOC100133331"; gene_name "LOC100133331"; transcript_id "NR_028327_1"; tss_id "TSS14025";


Specifically, the `attribute` field must, at a minimum, contain either a `gene_name` or `gene_id` attribute so that we can pair exons to uniquely-identified genes.
Also, .gtf records will be filtered down to those whose `feature` = "exon".


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
`-o`, `--output-dir` | No | Output directory. Default is current working directory. Directory will be created if it doesn't currently. If not supplied or if directory already exists, DegNorm creates a new output directory `degnorm_mmddYYYY_HHMMSS`.
`--plot-genes` | No | Names of genes for which to render coverage plots. Sequence of explictly stated gene names or a .txt file containing one gene name per line.
`-d`, `--downsample-rate` | No | Integer downsampling rate. Systematic samples of a coverage matrix are used to speed up NMF iterations.
`--nmf-iter` | No | Number of iterations per NMF-OA approximation. The higher the more accurate the approximation, but the more costly in terms of time.
`--iter` | No | Number of whole DegNorm iterations. Default is 5.
`--minimax-coverage` | No | Minimum cross-sample maximum coverage for a gene before it is included in the DegNorm pipeline. Can be used to exclude relatively low-coverage genes.
 `-s`, `--skip-baseline-selection` | No | Flag to skip baseline selection, will greatly speed up DegNorm iterations.
 `--non-unique-alignments` | No | Flag, allow non-uniquely mapped reads. Otherwise, DegNorm only keeps reads with `NH` (number of hits) == 1 (default behavior).
 `-p`, `--proc-per-node` | No | Integer number of processes to spawn per compute node. The more the better. Defaults to a lonely 1.


## Example usage


### Required: sorting, indexing alignment files

Remember - prior to running `degnorm`, you must `sort` your aligned reads and create index files (.bai files).
Suppose we have set of .bam files in our current working directory that we would like to send through DegNorm, the following
commands will sort the .bam files and create .bai files, all with `samtools`:

    # sort and index alignment files in current working directory.
    for FILE in ./*.bam
    do
        echo 'sorting '$FILE
        samtools sort $FILE -o ${FILE/.bam/}'_sorted.bam'
        echo 'indexing '${FILE/.bam/}'_sorted.bam'
        samtools index ${FILE/.bam/}'_sorted.bam' ${FILE/.bam/}'_sorted.bai'
    done

### Launching a full `degnorm` run "the verbose way"

Next, after the alignment files have been sorted (and optionally indexed), let's run `degnorm` "the verbose way" by enumerating 
each .bam and .bai file individually. Suppose we have two sorted alignment files, `S1_sorted.bam` and `S2_sorted.bam`, and that we want
to save the results of the pipeline to a new directory called `degnorm_output`.
    
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

It's probably just easier to supply `degnorm` with an input data directory containing all of our .bam files, instead of each .bam file individually.
Suppose that data directory is located at `degnorm_data/GBM`.
 Let's also use 5 DegNorm iterations, 50 NMF iterations per gene per NMF iteration, and route output to a directory besides the current working directory.

    degnorm --bam-dir degnorm_data/GBM \
        -g human.gtf \
        -p 20 \
        --nmf-iter 50 \
        -o degnorm_output

### Using a warm start directory

After one pipeline run, we could start `degnorm` from the output directory resulting from the first run - this directory is known as the "warm start" directory. This time, let's exclude genes with a maximum coverage
(across all samples) less than 20, and run DegNorm with only 3 iterations. In this example, the `DegNorm_102118_081045` warm start directory would have been generated from a prior DegNorm run.

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
        
 The functionality of `degnorm_mpi` is the same as the single-node `degnorm` command, so you can continue using warm start directories just like you could with `degnorm`.
 If you are unfamiliar with MPI, you may find [Wes Kendall's tutorials](http://mpitutorial.com/tutorials/) very enlightening.


## Example job submission script
Below is an example of a shell script containing the contents of a job you might submit to a job manager (e.g. MOAB).
Note that additional steps would be required to convert this script into a SLURM job submission script.

<script src="https://gist.github.com/ffineis/0a93d87519c64c1d0163a8eb1403bb2c.js"></script>
        
## Output

Raw and estimated coverage data are stored in separate `.pkl` files, one file per chromosome. See the [posthoc analysis](../howtos/posthoc_analysis.md)
 documentation for helper functions to access coverage data.

A gene that has 100% missing coverage (no coverage found in any RNA-seq sample) will not be run through DegNorm for numerical stability purposes and will therefore not have a
coverage matrix stored in any of the `.pkl` files.

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
    ├── chr1
    │   ├── GAPDH_coverage.png
    │   ├── <more coverage plots>
    │   ├── coverage_matrices_chr1.pkl  # (gene, raw coverage matrix) pairs in serialized python dictionary.
    │   └── estimated_coverage_matrices_chr1.pkl # (gene, estimated coverage matrix) pairs in serialized python dictionary.
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