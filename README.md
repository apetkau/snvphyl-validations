WGS SNVPhyl Whole Genome Phylogeny Simulations
==============================================

Generates simulated variants and reads for validating the SNVPhyl pipeline.

Dependencies
============

* [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
* Perl
* Perl modules: `cpanm Bio::SeqIO Set::Scalar`

Usage
=====

To generate simulation cases please run the following:

```
sh generate-simulations.sh
```

This will randomly mutate genomes under `references/` and generate simulated reads under `simulations/[reference]`.  fastq files are under `simulations/[reference]/fastq.  Please run fastq files through SNVPhyl pipeline.

A table of the true variants is found under `simulations/[reference]/variants.tsv`.  To compare with generated variants table from SNVPhyl please run:

```
perl compare_positions.pl simulations/[reference]/variants.tsv [snvphyl variants.tsv] | column -t
```

To find exact differences, please run:

```
diff simulations/[reference]/variants.tsv [snvphyl variants.tsv]
```
