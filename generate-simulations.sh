#!/bin/bash

usage="$0 [galaxy_url] [galaxy_api_key]"

galaxy_url=$1
galaxy_api_key=$2

if [ "$galaxy_url" == "" ];
then
        echo "Error: no galaxy_url found"
        echo $usage
        exit 1
elif [ "$galaxy_api_key" == "" ];
then
        echo "Error: no galaxy_api_key found"
        echo $usage
        exit 1
fi

genome=e_coli_sakai_w_plasmids
dir=simulations
reference=references/e_coli_sakai_w_plasmids.fasta
invalid_positions_file=references/invalid_positions.bed
simulation_dir=$dir/$genome
repeats_file=$simulation_dir/repeats-$genome.tsv
fastq_dir=$simulation_dir/fastq
variant_table=$simulation_dir/variants-$genome.tsv
case_string=`basename $simulation_dir`
echo "Case $case_string"
if [ -e $simulation_dir ];
then
	echo "Simulation dir $simulation_dir already exists, do not re-generate"
else
	mkdir $simulation_dir
	cp $reference $simulation_dir
	echo "Generate repeat regions list"
	perl scripts/find-repeats.pl -l 150 -p 90 $reference > $repeats_file
	echo "Generate variants table"
	perl scripts/generate_variant_table.pl --reference $reference --num-duplicate-reference-genomes 1 --num-variant-genomes 10 --num-substitutions 9800 --num-insertions 100 --num-deletions 100 --random-seed 42 --repeat-positions $repeats_file | sort -k 1,1 -k 2,2n > $variant_table
	echo "Generate reads. See $simulation_dir/generate_genomes.log for details"
	perl scripts/generate_genomes.pl --reference $reference --variants $variant_table --min-cov 30 --max-cov 30 --random-seed 42 --art-parameters '--noALN --len 250 --rndSeed 43 --paired --seqSys MS --sdev 100 --mflen 500' --out-dir $fastq_dir > $simulation_dir/generate_genomes.log
	
	mkdir $fastq_dir/log
	mv $fastq_dir/*.log $fastq_dir/log
	
	echo "Done $case_string"
	echo "Repeats in $repeats_file"
	echo "List of true variants in $variant_table"
	echo "Reads in $fastq_dir/*.fastq"
fi

snvphyl_run_dir=$simulation_dir/snvphyl-runs
if [ -e $snvphyl_run_dir ];
then
        echo "$snvphyl_run_dir already exists, will not run"
else
	echo "Executing workflow for $snvphyl_run_dir"
	mkdir $snvphyl_run_dir

	alternative_allele_ratio=0.75
	min_coverage=10
	run_name="run-$min_coverage-$alternative_allele_ratio"
	output_dir=$snvphyl_run_dir/output-$min_coverage-$alternative_allele_ratio
	log_out=$snvphyl_run_dir/$run_name-log.out
	log_err=$snvphyl_run_dir/$run_name-log.err
	variants_summary=$output_dir/variants_comparison_summary.tsv
	command="run-snphyl.py --galaxy-url $galaxy_url --galaxy-api-key $galaxy_api_key --fastq-dir $fastq_dir --reference-file $reference --invalid-positions-file $invalid_positions_file --run-name $run_name --alternative-allele-ratio $alternative_allele_ratio --min-coverage $min_coverage --output-dir $output_dir"
	date
	echo $command
	echo "Run 'tail -f $log_out $log_err' for more details"
	$command > $log_out 2> $log_err

	perl scripts/compare_positions.pl --variants-true $variant_table --variants-detected $output_dir/$run_name-pseudo-positions.tsv --reference-genome $reference > $variants_summary	
fi
