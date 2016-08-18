SNV Density Filtering Validation
================================

Files
=====

* `FM211187.fasta` - reference genome
* `PMEN1.aln` - Original alignment, from <ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz>
* `pmen1.name_accession` - Table of genome names and accession numbers for dataset, extracted from Table S1 in <http://science.sciencemag.org/content/331/6016/430.full>.
* `pmen1.err` - Table of run ids for accessions in `pmen1.name_accession`.  Generated using [SRAdb](https://bioconductor.org/packages/release/bioc/html/SRAdb.html).
* `fastq/` - Concatnated fastq files.  Downloaded from NCBI using SRR ids from `pmen1.srr` and concatenated/re-named from names in `pmen1.name_accession`.
