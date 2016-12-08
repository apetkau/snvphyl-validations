SNV Density Filtering Validation
================================

This describes the procedures used for SNVPhyl's SNV density evaluation against the [Gubbins](https://github.com/sanger-pathogens/gubbins) software package.

Files
=====

* `FM211187.fasta` - reference genome
* `PMEN1.aln.gz` - Original alignment, from <ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz>
* `pmen1.name_accession` - Table of genome names and accession numbers for dataset, extracted from Table S1 in <http://science.sciencemag.org/content/331/6016/430.full>.
* `pmen1.err` - Table of run ids for accessions in `pmen1.name_accession`.  Generated using [SRAdb](https://bioconductor.org/packages/release/bioc/html/SRAdb.html).
* `fastq/` - Concatnated fastq files.  Downloaded from NCBI using ERR ids from `pmen1.err` and concatenated/re-named from names in `pmen1.name_accession`.

Scripts
=======

1. [../scripts/gubbinsSnps2Table.pl](../scripts/gubbinsSnps2Table.pl): Converts Gubbins VCF SNP/SNV file to SNVPhyl SNV table format.
2. [../scripts/compare_positions.pl](../scripts/compare_positions.pl): Compares converted variant table from Gubbins and table produced by SNVPhyl.  Counts TP/FP/TN/FN variants.
3. [plot_trees.R](plot_trees.R): Constructs figure comparing phylogenetic trees.

Dependencies
============

* [Gubbins](https://github.com/sanger-pathogens/gubbins)
* Perl modules: `cpanm Bio::SeqIO Set::Scalar`
* [SNVPhyl command-line-interface](https://github.com/phac-nml/snvphyl-galaxy-cli)
* [Docker](https://www.docker.com/)
* [PhyML](http://www.atgc-montpellier.fr/phyml/)
* [Ktreedist](http://molevol.cmima.csic.es/castresana/Ktreedist.html)
* <https://github.com/phac-nml/snvphyl-tools/blob/ff57489703be4ba716eb5468b6de82c808f556ad/positions2snv_invariant_alignment.pl>
* R and R modules [APE](http://ape-package.ird.fr/) and [Phytools](https://github.com/liamrevell/phytools).

Original Gubbins Results
========================

```bash
gunzip PMEN1.aln.gz
cat FM211187.fasta PMEN1.aln > PMEN1-with-reference.aln # concatenate reference genome back onto alignment
sed -i -e 's/gi|220673408|emb|FM211187.1| Streptococcus pneumoniae ATCC 700669 complete genome/reference/' PMEN1-with-reference.aln

mkdir original_gubbins/ && cd original_gubbins/
run_gubbins.py ../PMEN1-with-reference.aln

# Remove N/- characters from SNVs
head -n 4 original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf > original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash
tail -n+5 original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf | grep -v '[-N]' >> original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash

# Generate phylip-formatted alignment from SNVs after removal of N/-
sed -e 's/N/-/g' original_gubbins/PMEN1-with-reference.filtered_polymorphic_sites.fasta | perl -MBio::AlignIO -e '$i=Bio::AlignIO->new(-fh=>\*STDIN,-format=>"fasta");$o=Bio::AlignIO->new(-file=>">original_gubbins/PMEN1-with-reference.filtered_polymorphic_sites.phylip.noNDashes",-format=>"phylip");print $i;$a=$i->next_aln->remove_gaps("-");$o->write_aln($a);'
```

SNVPhyl no SNV-density filtering
================================

## Run SNVPhyl

```bash
# Disable SNV density filtering (run with threshold much larger than window size)
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 1 --filter-density-threshold 20 --output-dir snvphyl-no-filter
```

## Analyze Results

Generate the (TP,FP,TN,FN) numbers which are compiled into a table later.  Case **No SNV filtering** comes from *all-vs-valid* row.

```bash
cd snvphyl-no-filter

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

# For numerical comparisons (TP/FP, etc)
perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl With Filter 2 in 20
===========================

Assumes [SNVPhyl Galaxy CLI](https://github.com/phac-nml/snvphyl-galaxy-cli) is installed. I ran with pipeline on Galaxy instance connected to cluster, but same results can be repeated with Docker (`--deploy-docker`).

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 20 --filter-density-threshold 2 --output-dir snvphyl-2-20
```
## Analyze Results

```bash
cd snvphyl-2-20

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

# For numerical comparisons
perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl With Filter 2 in 100
============================

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 100 --filter-density-threshold 2 --output-dir snvphyl-2-100
```

## Analyze Results

```bash
cd snvphyl-2-100

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl With Filter 2 in 500
============================

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 500 --filter-density-threshold 2 --output-dir snvphyl-2-500
```

## Analyze Results

```bash
cd snvphyl-2-500

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl With Filter 2 in 1000
=============================

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 1000 --filter-density-threshold 2 --output-dir snvphyl-2-1000
```

## Analyze Results

```bash
cd snvphyl-2-1000

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl With Filter 2 in 2000
=============================

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --min-coverage 10 --filter-density-window 2000 --filter-density-threshold 2 --output-dir snvphyl-2-2000
```

## Analyze Results

```bash
cd snvphyl-2-2000

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl then Gubbins
====================

To test using an alignment from SNVs detected with SNVPhyl from Gubbins, the SNVPhyl results from no SNV density filtering was used to generate an alignment for Gubbins valid alignment positions detected by SNVPhyl.

## Analyze Results

```bash
cd snvphyl-no-filter

# Construct alignment of valid SNV alignment columns and invariant sites
# Uses script from https://github.com/phac-nml/snvphyl-tools/blob/ff57489703be4ba716eb5468b6de82c808f556ad/positions2snv_invariant_alignment.pl
positions2snv_invariant_alignment.pl -i snvTable.tsv -o invariant-alignment -f fasta --reference-file ../FM211187.fasta
cd invariant-alignment/
mv "gi|220673408|emb|FM211187.1|.fasta" no-density-alignment.fasta

# Run Gubbins
run_gubbins.py no-density-alignment.fasta

# Remove Ns from results
head -n 4 no-density-alignment.summary_of_snp_distribution.vcf > no-density-alignment.summary_of_snp_distribution.noNs.vcf
tail -n+5 no-density-alignment.summary_of_snp_distribution.vcf | grep -v '[N]' >> no-density-alignment.summary_of_snp_distribution.noNs.vcf

# Generate alignment with no Ns
# note, there is no gap characters (grep '-' no-density-alignment.filtered_polymorphic_sites.fasta) before running this, so we are only removing N
sed -e 's/N/-/g' no-density-alignment.filtered_polymorphic_sites.fasta | perl -MBio::AlignIO -e '$i=Bio::AlignIO->new(-fh=>\*STDIN,-format=>"fasta");$o=Bio::AlignIO->new(-file=>">no-density-alignment.filtered_polymorphic_sites.noNs.phylip",-format=>"phylip");print $i;$a=$i->next_aln->remove_gaps("-");$o->write_aln($a);'

# Compare results to original alignment used by Gubbins
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table ../../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf.removeNDash | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.removeNDash.tsv
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table no-density-alignment.summary_of_snp_distribution.noNs.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > no-density-alignment.ordered.summary_of_snp_distribution.noNs.tsv

perl ../../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.removeNDash.tsv --variants-detected no-density-alignment.ordered.summary_of_snp_distribution.noNs.tsv --reference-genome ../../FM211187.fasta --false-detection-output false | column -t

# For results including Ns/dashes
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table ../../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.withNDashes.tsv

perl ../../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.withNDashes.tsv --variants-detected no-density-alignment.ordered.summary_of_snp_distribution.noNs.tsv --reference-genome ../../FM211187.fasta --false-detection-output false | column -t
```

Calculating Tree distances
==========================

All phylogenetic trees were re-generated with phyml in order to compute distances between them.

```bash
mkdir tree-distances && cd tree-distances

# Copy alignment files
for i in ../snvphyl-*; do b=`basename $i`; cp $i/snvAlignment.phy $b.snvAlignment.phy; done
cp ../original_gubbins/PMEN1-with-reference.filtered_polymorphic_sites.phylip.noNDashes original_gubbins.phy
cp ../original_gubbins/PMEN1-with-reference.filtered_polymorphic_sites.phylip original_gubbins-withNDashes.phy
cp ../snvphyl-no-filter/invariant-alignment/no-density-alignment.filtered_polymorphic_sites.noNs.phylip snvphyl-gubbins.phy

# change all 'Reference' names in alignment to lower-case
sed -i -e 's/Reference/reference/' *.phy

# Generate trees with phyml (version 20131022)
for i in *.phy; do phyml -i $i -s BEST -m GTR --r_seed 42; done
```

Compare trees with [Ktreedist](http://molevol.cmima.csic.es/castresana/Ktreedist.html) version 1.0.

```bash
for i in *_tree.txt; do Ktreedist.pl -rt original_gubbins.phy_phyml_tree.txt -ct $i -t $i.tree.tsv; done
Ktreedist.pl -rt original_gubbins-withNDashes.phy_phyml_tree.txt -ct snvphyl-gubbins.phy_phyml_tree.txt -t snvphyl-gubbins.phy_phyml_tree.txt.tree.withNDashes.tsv
cat *.tree*.tsv | sort -ur | column -ts $'\t'
```


Compiling Results
=================

All results from `all-vs-valid` comparisons to Gubbins for each case were compiled into a table.

Figure S2 (plot of all trees) constructed from:

```
R CMD BATCH plot_trees.R
```
