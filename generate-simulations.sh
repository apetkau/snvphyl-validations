#!/bin/sh

genome=2011C-3609
dir=simulations
reference=references/2011C-3609.fasta
simulation_dir=$dir/$genome
repeats_file=$simulation_dir/repeats.tsv
fastq_dir=$simulation_dir/fastq
variant_table=$simulation_dir/variants.tsv
case_string=`basename $simulation_dir`
echo "Case $case_string"
if [ -e $simulation_dir ];
then
	echo "Simulation dir $simulation_dir already exists, do not re-generate"
else
	mkdir $simulation_dir
	cp $reference $simulation_dir
	echo "Generate repeat regions list"
	perl scripts/find-repeats.pl $reference > $repeats_file
	echo "Generate variants table"
	perl scripts/generate_variant_table.pl --reference $reference --num-genomes 10 --num-substitutions 10000 --num-insertions 100 --num-deletions 100 --random-seed 42 --exclude-positions $repeats_file | sort -k 1,1 -k 2,2n > $variant_table
	echo "Generate reads. See $simulation_dir/generate_genomes.log for details"
	perl scripts/generate_genomes.pl --reference $reference --variants $variant_table --min-cov 30 --max-cov 60 --random-seed 43 --art-parameters '--noALN --len 250 --rndSeed 43 --paired --seqSys MS --sdev 100 --mflen 500' --out-dir $fastq_dir > $simulation_dir/generate_genomes.log
	
	mkdir $fastq_dir/log
	mv $fastq_dir/*.log $fastq_dir/log
	
	echo "Done $case_string"
	echo "Repeats in $repeats_file"
	echo "List of true variants in $variant_table"
	echo "Reads in $fastq_dir/*.fastq"
fi
