SNVPhyl simulated data validation
=================================

This repository contains code for generating read simulations and comparing with SNVPhyl results.  Instructions to reproduce results are as follows.

```bash
sh generate-simulations.sh
sh run-simulations.sh
```

Results will appear in `simulations/e_coli_sakai_w_plasmids/`, in particular the file `simulations/e_coli_sakai_w_plasmids/snvphyl-runs/output-10-0.75/variants_comparison_summary.tsv`, comparing the variants simulated to the variants detected.

To find the misscalled bases and construct Table S3, the following commands were run.

```bash
# Search for base differences (excluding - and N) in the results.
grep -P 'all-vs-all\tFP' simulations/e_coli_sakai_w_plasmids/snvphyl-runs/output-10-0.75/false_variants.tsv

# Search for specific bases for all simulated genomes in list of positions from command above.  E.g.
grep 345916 simulations/e_coli_sakai_w_plasmids/variants-e_coli_sakai_w_plasmids.tsv simulations/e_coli_sakai_w_plasmids/snvphyl-runs/output-10-0.75/snvTable.tsv

# Use these results to fill in Table S3.
```
To determine the copy numbers of the covering regions for each of the above variants the below commands were run.

```
cd simulations/e_coli_sakai_w_plasmids/fastq/genomes
for i in *.fasta; do nucmer --maxmatch --nosimplify --prefix $i ../../e_coli_sakai_w_plasmids.fasta $i; done

# for each position and genome, run show-snps and record copies covering position
for p in 345916 2205013 4159150; do for g in *.fasta.delta; do echo "Position $p, Genome $g"; show-snps -l -r $g | grep "$p  *[ATCG] [ATCG]  *$p"; done; done
```

Scripts
======= 

These simulations make use of some scripts, mainly:

1. [scripts/generate_variant_table.pl](scripts/generate_variant_table.pl): Generates a table of random variants.
2. [scripts/generate_genomes.pl](scripts/generate_genomes.pl): Constructs mutated reference genomes from the variant table and the genome under [references/](references/).  Simulates reads using `art_illumina`.
3. [scripts/compare_positions.pl](scripts/compare_positions.pl): Compares initially constructed variant table and table produced by SNVPhyl.  Counts TP/FP/TN/FN variants.
4. [generate-simulations.sh](generate-simulations.sh): Generates the read simulations.
5. [run-simulations.sh](run-simulations.sh): Runs simulations.

Dependencies
============

* [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
* Perl
* Perl modules: `cpanm Bio::SeqIO Set::Scalar`
* [SNVPhyl command-line-interface](https://github.com/phac-nml/snvphyl-galaxy-cli)
* [Docker](https://www.docker.com/)
* [MUMMer](http://mummer.sourceforge.net/)

