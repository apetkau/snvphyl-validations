SNVPhyl Validations
===================

This repository contains the code and describes the steps followed for the three validations of the [SNVPhyl](http://snvphyl.readthedocs.io) pipeline as described in the manuscript.  These consist of:

1. [Simulated data](README-simulations.md): Validation of SNVPhyl against simulated data.
2. [SNV density filtering](snv-density-filtering/README.md): Evaluation of SNVPhyl's SNV density filtering against the [Gubbins](https://github.com/sanger-pathogens/gubbins) software.
3. [Parameter optimization](salmonella_heidelberg/README.md): Evaluation and optimization of SNVPhyl's parameters against a real-world dataset.


Most of these instructions made use of the [SNVPhyl command-line interface](https://github.com/phac-nml/snvphyl-galaxy-cli) and assume that `snvphyl.py` is on the `PATH`.
