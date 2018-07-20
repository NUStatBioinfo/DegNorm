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
    sample_files = args.input_files
    n_samples = len(sample_files)
    sample_ids = list()
    chroms = list()
    cov_files = dict()

    output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))
    if os.path.isdir(output_dir):
        time.sleep(3)
        output_dir = os.path.join(args.output_dir, 'DegNorm_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    os.makedirs(output_dir)

    # ---------------------------------------------------------------------------- #
    # Load .sam or .bam files; parse into coverage arrays and store them.
    # ---------------------------------------------------------------------------- #

    # iterate over .sam files; compute each sample's chromosomes' coverage arrays
    # and save them to .npz files.
    for idx in range(n_samples):
        logging.info('Loading RNA-seq data file {0} / {1}'.format(idx + 1, n_samples))
        if args.input_type == 'sam':
            reader = ReadsCoverageProcessor(sam_file=sample_files[idx]
                                            , n_jobs=args.cpu
                                            , tmp_dir=output_dir
                                            , verbose=True)

        sample_id = reader.sample_id
        sample_ids.append(sample_id)
        cov_files[sample_id] = reader.coverage()
        chroms += [os.path.basename(f).split('.npz')[0].split('_')[-1] for f in cov_files[sample_id]]

    chroms = list(set(chroms))
    n_chroms= len(chroms)
    logging.info('Total chromosomes: {0}'.format(n_chroms))

    # ---------------------------------------------------------------------------- #
    # Load .gtf or .gff files and run processing pipeline.
    # ---------------------------------------------------------------------------- #

    # load and process gene annotation file
    logging.info('Begin genome annotation file processing...')
    gap = GeneAnnotationProcessor(args.genome_annotation
                                  , n_jobs=n_jobs
                                  , verbose=True)
    exon_df = gap.run()

    # ---------------------------------------------------------------------------- #
    # Slice up coverage matrix for each gene according to exon positioning.
    # ---------------------------------------------------------------------------- #

    cov_output = output_dir if args.save_coverage else None
    p = mp.Pool(processes=n_jobs)
    gene_cov_mats = [p.apply_async(gene_coverage
                                   , args=(exon_df, chroms[i], cov_files, sample_ids, cov_output)) for i in range(n_chroms)]
    p.close()
    gene_cov_mats = [x.get() for x in gene_cov_mats]