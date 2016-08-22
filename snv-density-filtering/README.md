SNV Density Filtering Validation
================================

Files
=====

* `FM211187.fasta` - reference genome
* `PMEN1.aln` - Original alignment, from <ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz>
* `pmen1.name_accession` - Table of genome names and accession numbers for dataset, extracted from Table S1 in <http://science.sciencemag.org/content/331/6016/430.full>.
* `pmen1.err` - Table of run ids for accessions in `pmen1.name_accession`.  Generated using [SRAdb](https://bioconductor.org/packages/release/bioc/html/SRAdb.html).
* `fastq/` - Concatnated fastq files.  Downloaded from NCBI using ERR ids from `pmen1.err` and concatenated/re-named from names in `pmen1.name_accession`.

Original Gubbins Results
========================

```bash
cat FM211187.fasta PMEN1.aln > PMEN1-with-reference.aln # concatenate reference genome back onto alignment
sed -i -e 's/gi|220673408|emb|FM211187.1| Streptococcus pneumoniae ATCC 700669 complete genome/reference/' PMEN1-with-reference.aln

mkdir original_gubbins/ && cd original_gubbins/
run_gubbins.py ../PMEN1-with-reference.aln
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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

# For numerical comparisons
perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t

# For visual/diff comparisons
head -n 1 snvTable.tsv > snvTable.valid.tsv && grep -P '\tvalid\t' snvTable.tsv >> snvTable.valid.tsv # generate table of only valid SNVPhyl results

meld PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv snvTable.valid.tsv
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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

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

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

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

# Compare results to original alignment used by Gubbins
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table ../../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table no-density-alignment.summary_of_snp_distribution.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > no-density-alignment.ordered.summary_of_snp_distribution.tsv

perl ../../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected no-density-alignment.ordered.summary_of_snp_distribution.tsv --reference-genome ../../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl then Gubbins (with all positions)
=========================================

To test using an alignment from SNVs detected with SNVPhyl from Gubbins, the SNVPhyl results from no SNV density filtering was used to generate an alignment for Gubbins with all alignment positions detected by SNVPhyl.

## Analyze Results

```bash
cd snvphyl-no-filter

# Construct alignment of all SNV alignment columns and invariant sites
# Uses script from https://github.com/phac-nml/snvphyl-tools/blob/ff57489703be4ba716eb5468b6de82c808f556ad/positions2snv_invariant_alignment.pl
positions2snv_invariant_alignment.pl -i snvTable.tsv -o invariant-alignment-all -f fasta --reference-file ../FM211187.fasta --keep-all
cd invariant-alignment-all/
mv "gi|220673408|emb|FM211187.1|.fasta" no-density-alignment.fasta

# Run Gubbins
run_gubbins.py no-density-alignment-all.fasta

# Compare results to original alignment used by Gubbins
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table ../../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv
perl ../../../scripts/gubbinsSnps2Table.pl --snvphyl-table ../snvTable.tsv --gubbins-table no-density-alignment-all.summary_of_snp_distribution.vcf | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > no-density-alignment-all.ordered.summary_of_snp_distribution.tsv

perl ../../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected no-density-alignment-all.ordered.summary_of_snp_distribution.tsv --reference-genome ../../FM211187.fasta --false-detection-output false | column -t
```

SNVPhyl with Phage/Genomic Islands removed
==========================================

To test the case of SNVPhyl results with no SNV density filtering, but phage/genomic islands removed, we perfomed the following.

## Generate SNVPhyl invalid-positions file

1. Genomic islands were downloaded from [IslandViewer](http://www.pathogenomics.sfu.ca/islandviewer/) for FM211187.1 (NC_011900.1) and stored as `NC_011900.1.islandviewer.tsv`
2. Phage locations downloaded from [PHASTER](http://phaster.ca/) for FM211187.1 (NC_011900.1) and stored as NC_011900.1.phaster.txt
3. Locations complied together into `invalid_positions.tsv` file.

## Run SNVPhyl

```bash
snvphyl.py --deploy-docker --fastq-dir fastq/ --reference-file FM211187.fasta --invalid-positions invalid_positions.tsv --min-coverage 10 --filter-density-window 1 --filter-density-threshold 20 --output-dir snvphyl-no-filter-remove-phage-islands
```

## Analyse Results

```bash
cd snvphyl-no-filter-remove-phage-islands

perl ../../scripts/gubbinsSnps2Table.pl --snvphyl-table snvTable.tsv --gubbins-table ../original_gubbins/PMEN1-with-reference.summary_of_snp_distribution.vcf --mark-invalid-alignments | sed -e 's/{chrom}/gi|220673408|emb|FM211187.1|/' > PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv

perl ../../scripts/compare_positions.pl --variants-true PMEN1-with-reference.ordered.summary_of_snp_distribution.tsv --variants-detected snvTable.tsv --reference-genome ../FM211187.fasta --false-detection-output false | column -t
```

Calculating Tree distances
==========================

All phylogenetic trees were re-generated with phyml in order to compute distances between them.

```bash
mkdir tree-distances && cd tree-distances

# Copy alignment files
for i in ../snvphyl-*; do b=`basename $i`; cp $i/snvAlignment.phy $b.snvAlignment.phy; done
cp ../original_gubbins/PMEN1-with-reference.filtered_polymorphic_sites.phylip original_gubbins.phy
cp ../snvphyl-no-filter/invariant-alignment/no-density-alignment.filtered_polymorphic_sites.phylip snvphyl-gubbins.phy
cp ../snvphyl-no-filter/invariant-alignment-all/no-density-alignment-all.filtered_polymorphic_sites.phylip snvphyl-gubbins-all.phy

# change all 'Reference' names in alignment to lower-case
sed -i -e 's/Reference/reference/' *.phy

# Generate trees with phyml (version 20131022)
for i in *.phy; do phyml -i $i -s BEST -m GTR --r_seed 42; done

# Distances (with treedist from phylip package)
mkdir distances
cp original_gubbins.phy_phyml_tree.txt distances/intree
n=2; for i in original_gubbins.phy_phyml_tree.txt snvphyl*phyml_tree.txt; do b=`basename $i .phy_phyml_tree.txt`; echo -e "$n\t$b" >> distances/number_name.tsv; cat $i >> distances/intree2; n=$(($n+1)); done
cd distances
echo -e "2\nL\nV\nY" | phylip treedist | tee settings
for i in `cut -f 1 number_name.tsv`; do j=`grep "^$i" number_name.tsv | sed -e 's/^[0-9]\t//'`; sed -i -e "s/Trees 1 and $i/Trees 1 and $j/" outfile; done
```

Compare trees with [Ktreedist](http://molevol.cmima.csic.es/castresana/Ktreedist.html) version 1.0.

```bash
for i in *_tree.txt; do Ktreedist.pl -rt original_gubbins.phy_phyml_tree.txt -ct $i -t $i.tree.tsv; done
cat *.tree.tsv | sort -ur | column -ts $'\t' # prints table of distances from ktreedist
```

Compiling Results
=================

All results from `all-vs-valid` comparisons to Gubbins for each case were compiled into a table `snvphyl-gubbins.xlsx`.
