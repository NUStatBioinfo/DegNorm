## DegNorm TODO

### For ``v0.1.3`` release
- Improve upon the memory hogging in `coverage.gene_coverage` - dice up sparse
chromosome-wide coverage matrix into chunks prior to .todense(),
then run .todense + gene coverage splicing on chunks, while ensuring that every gene
still gets processed even if a gene spans multiple chunks.