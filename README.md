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

This will randomly mutate genomes under `references/` and generate simulated reads under `simulations/[reference]`.  fastq files are under `simulations/[reference]/fastq.  Please run fastq files through SNVPhyl pipeline.

A table of the true variants is found under `simulations/[reference]/variants.tsv`.  To compare with generated variants table from SNVPhyl please run:

```
perl compare_positions.pl simulations/[reference]/variants.tsv [snvphyl variants.tsv] | column -t
```

For example:

```
perl scripts/compare_positions.pl variants-true.tsv variants-detected.tsv |column -t
variants-true.tsv  variants-detected.tsv  Intersection  Unique-variants-true.tsv  Unique-variants-detected.tsv
1000               947                    947           53                        0
```

To find exact differences, please run:

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
* Mutations are generated for random positions.  For any positions in repetitive regions, SNVPhyl will not detect and this will show up as missing variants.
* ART assumes reads produced by a MiSeq all have the same length.  This is not the case for real-world data.
