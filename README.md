SNVPhyl Validations
===================

This repository contains the code and describes the steps followed for the three validations of the [SNVPhyl](http://snvphyl.readthedocs.io) pipeline as described in:

Petkau A, Mabon P, Sieffert C, Knox N, Cabral J, Iskander M, Iskander M, Weedmark K, Zaheer R, Katz L, Nadon C, Reimer A, Taboada E, Beiko R, Hsiao W, Brinkman F, Graham M, Van Domselaar G. [SNVPhyl: a single nucleotide variant phylogenomics pipeline for microbial genomic epidemiology](http://dx.doi.org/10.1099/mgen.0.000116). 08/06/2017. *M Gen* (online ahead of inclusion in an Issue). doi: [10.1099/mgen.0.000116](https://doi.org/10.1099/mgen.0.000116).

These consist of:

1. [Simulated data](README-simulations.md): Validation of SNVPhyl against simulated data.
2. [SNV density filtering](snv-density-filtering/README.md): Evaluation of SNVPhyl's SNV density filtering against the [Gubbins](https://github.com/sanger-pathogens/gubbins) software.
3. [Parameter optimization](salmonella_heidelberg/README.md): Evaluation and optimization of SNVPhyl's parameters against a real-world dataset.


Most of these instructions made use of the [SNVPhyl command-line interface](https://github.com/phac-nml/snvphyl-galaxy-cli) and assume that `snvphyl.py` is on the `PATH`.
