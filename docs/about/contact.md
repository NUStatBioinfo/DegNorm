# DegNorm contact

### Email
- Northwestern Statistical Bioinformatics Lab: Dr. Jiping Wang (jzwang@northwestern.edu)
- Python package issues: Frank Fineis (frankfineis2022@u.northwestern.edu)
- R package issues: Dr. Jiping Wang (jzwang@northwestern.edu), Bin Xiong (binxiong2012@u.northwestern.edu)

### Python package issues and bug reporting

If you run into problems using the `degnorm` CLI, we recommend please ensuring that your .bam and .gtf files
are 0- and 1-indexed, respectively. Also please check that your .gtf input matches the format of the example .gtf in the [usage docs](../howtos/run_the_pipeline.md)
and that your `attribute` column includes a `gene_id` or `gene_name` tag for every gene.

If you're sure you've checked your input against the usage docs, [open an issue](https://github.com/NUStatBioinfo/DegNorm/issues) on the DegNorm GitHub repository! DegNorm's applications are
 wide, so we encourage users to open git issues with minimally reproducible examples. When opening an issue on GitHub, please provide the full stack trace for debugging purposes.

Also feel free to submit pull requests for additional feature support or fixes.