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
	false_variants=$output_dir/false_variants.tsv
	command="run-snphyl.py --galaxy-url $galaxy_url --galaxy-api-key $galaxy_api_key --fastq-dir $fastq_dir --reference-file $reference --invalid-positions-file $invalid_positions_file --run-name $run_name --alternative-allele-ratio $alternative_allele_ratio --min-coverage $min_coverage --filter-density-window 0 --filter-density-threshold 100 --output-dir $output_dir"
	date
	echo $command
	echo "Run 'tail -f $log_out $log_err' for more details"
	$command > $log_out 2> $log_err

	perl scripts/compare_positions.pl --variants-true $variant_table --variants-detected $output_dir/$run_name-pseudo-positions.tsv --reference-genome $reference --false-detection-output $false_variants > $variants_summary
fi
