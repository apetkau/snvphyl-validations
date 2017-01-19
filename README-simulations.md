SNVPhyl simulated data validation
=================================

This repository contains code for generating read simulations and comparing with SNVPhyl results.  Instructions to reproduce results are as follows.

```bash
sh generate-simulations.sh
sh run-simulations.sh
```

Results will appear in `simulations/e_coli_sakai_w_plasmids/`, in particular the file `simulations/e_coli_sakai_w_plasmids/snvphyl-runs/output-10-0.75/variants_comparison_summary.tsv`.

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

