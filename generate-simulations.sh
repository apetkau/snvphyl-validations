#!/bin/sh

genome=08-5578
dir=simulations
reference=references/08-5578.fasta
simulation_dir=$dir/$genome
fastq_dir=$simulation_dir/fastq
variant_table=$simulation_dir/variants.tsv
echo "Case $genome"
if [ -e $simulation_dir ];
then
	echo "Simulation dir $simulation_dir already exists, do not re-generate"
	exit 1
fi

mkdir $simulation_dir
echo "Generate variants table"
perl scripts/generate_variant_table.pl --reference $reference --num-genomes 10 --num-variants 1000 --random-seed 42 | sort -k 1,1 -k 2,2n > $variant_table
echo "Generate reads"
perl scripts/generate_genomes.pl --reference $reference --variants $variant_table --min-cov 30 --max-cov 60 --random-seed 42 --art-parameters '--noALN --len 250 --rndSeed 42 --paired --seqSys MS --sdev 50 --mflen 400' --out-dir $fastq_dir | tee $simulation_dir/generate_genomes.log

mkdir $fastq_dir/log
mv $fastq_dir/*.log $fastq_dir/log

echo "Done $genome"
echo "List of true variants in $variant_table"
echo "Reads in $fastq_dir/*.fastq"
