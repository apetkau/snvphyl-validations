#!/bin/bash

# estimate coverages
(for i in fastqs/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee coverages.txt
