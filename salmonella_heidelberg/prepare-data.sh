#!/bin/sh

# downsample
for i in fastqs/*.fastq; do b=`basename $i`; echo $i; seqtk sample -s 66 $i 0.41 > fastqs-downsampled/$b; done
