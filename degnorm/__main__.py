from degnorm.reads_coverage import *
from degnorm.gene_coverage import *
from degnorm.gene_processing import *
from degnorm.utils import *
from datetime import datetime
import time


def main():

    # ---------------------------------------------------------------------------- #
    # Load CLI arguments and initialize storage, output directory
    # ---------------------------------------------------------------------------- #

    args = parse_args()
    n_jobs = args.cpu
    n_samples = len(args.input_files)
    sample_ids = list()
    chroms = list()
    cov_files = dict()

    output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    if os.path.isdir(output_dir):
        time.sleep(2)
        output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(output_dir)

    # ---------------------------------------------------------------------------- #
    # Load .sam or .bam files; parse into coverage arrays and store them.
    # ---------------------------------------------------------------------------- #

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in range(n_samples):
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))
        sam_file = args.input_files[idx]

        if args.input_type == 'bam':
            logging.info('Converting {0} into .sam file format...'
                         .format(sam_file))
            sam_file = bam_to_sam(sam_file)

        reader = ReadsCoverageProcessor(sam_file=sam_file
                                        , n_jobs=args.cpu
                                        , tmp_dir=output_dir
                                        , verbose=True)

        sample_id = reader.sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id] = reader.coverage()
        chroms += [os.path.basename(f).split('.npz')[0].split('_')[-1] for f in cov_files[sample_id]]

    chroms = list(set(chroms))
    logging.info('Obtained reads coverage for {0} chromosomes:\n'
                 '{1}'.format(len(chroms), ', '.join(chroms)))

    # ---------------------------------------------------------------------------- #
    # Load .gtf or .gff files and run processing pipeline.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #

    # load and process gene annotation file
    logging.info('Begin genome annotation file processing...')
    gap = GeneAnnotationProcessor(args.genome_annotation
                                  , n_jobs=n_jobs
                                  , verbose=True)
    exon_df = gap.run(chroms=chroms)
    
    # ---------------------------------------------------------------------------- #
    # Slice up coverage matrix for each gene according to exon positioning.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #

    cov_output = None if args.disregard_coverage else output_dir
    p = mp.Pool(processes=n_jobs)
    gene_cov_mats = [p.apply_async(gene_coverage
                                   , args=(exon_df, chrom, cov_files, cov_output)) for chrom in chroms]
    p.close()
    gene_cov_mats = [x.get() for x in gene_cov_mats]

    # convert list of tuples into 2-d dictionary: {chrom: {gene: coverage matrix}}
    chrom_gene_dict = {gene_cov_mats[i][1]: gene_cov_mats[i][0] for i in range(len(gene_cov_mats))}

    # ---------------------------------------------------------------------------- #
    # Run NMF.
    # Run in parallel over chromosomes.
    # ---------------------------------------------------------------------------- #



if __name__ == "__main__":
    main()