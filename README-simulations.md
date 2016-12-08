WGS SNVPhyl Whole Genome Phylogeny Simulations
==============================================

Generates simulated variants and reads for validating the SNVPhyl pipeline.  Simulated variants are stored in a variant table identical to the table produced by SNVPhyl for easy comparison:

**variants.tsv:**

```
#Chromosome                    Position  Status  Reference  08-5578-0  08-5578-1
gi|662858600|ref|NC_013766.2|  2280      valid   G          A          G        
gi|662858600|ref|NC_013766.2|  6156      valid   A          T          A       
gi|662858600|ref|NC_013766.2|  8071      valid   T          C          T      
```

The **variants.tsv** file is used to generate a set of variant reference genomes, in `references/` and fastq files are produced from the genomes using the ART read simulator with varying amounts of coverage.  The fastq files, along with the reference genome can be uploaded to Galaxy and the variants table produced by SNVPhyl can be compared to the original variants table.

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

This will randomly mutate genomes under `references/` and generate simulated reads under `simulations/[reference]`.  fastq files are under `simulations/[reference]/fastq`.  Please run fastq files through SNVPhyl pipeline.

A table of the true variants is found under `simulations/[reference]/variants.tsv`.  To compare with generated variants table from SNVPhyl please run:

```
perl compare_positions.pl --variants-true simulations/[reference]/variants.tsv --variants-detected [snvphyl variants.tsv] --core-genome [reference.fasta] | column -t
```

This will generate a table contining summaries of each file, true positives/negatives, along with calculations of specificity and sensitivity for the alignment.  The core genome is used for determining the number of true negatives and for our simulations is assumed to be equivalent to the reference genome.  For example:

```
perl scripts/compare_positions.pl --variants-true variants-true.tsv --variants-detected variants-detected.tsv --core-genome 08-5578.fasta |column -t
Core_Genome_File  Core_Genome  Variants_True_File  Variants_Detected_File  True_Variants  Variants_Detected  TP   FP  TN       FN  Accuracy  Specificity  Sensitivity  Precision  FP_Rate
08-5578.fasta     3032624      variants-true.tsv   variants-detected.tsv   1000           989                989  0   3031635  11  1.0000    1.0000       0.9890       1.0000     0.0000
```

To find exact differences, please run (you may need to sort the detected variants table first with `sort -k 1,1 -k 2,2n`):

```
diff simulations/[reference]/variants.tsv [snvphyl variants.tsv]
```

For example:

```
diff variants-true.tsv variants-detected.tsv 
18c18
< gi|662858600|ref|NC_013766.2| 45111   valid   G       C       G       C       G       C       G       C       G       C       G
---
> gi|662858600|ref|NC_013766.2| 45111   DP4     G       C       G       C       G       -       G       C       G       C       G
...
```

Note, differences will show up in repetitive regions due to masking those regions out in SNVPhyl.

Considerations
==============

A few considerations when using these simulated datasets include:

* Only very simple mutations are simulated.  No indels.
* ART assumes reads produced by a MiSeq all have the same length.  This is not the case for real-world data.
