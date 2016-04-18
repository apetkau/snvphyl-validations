#!/bin/bash

usage="$0 [galaxy_url] [galaxy_api_key]"

galaxy_url=$1
galaxy_api_key=$2

directory=cases

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

if [ -e $directory ];
then
	echo "Directory $directory already exists, will not re-run"
	exit 1
else
	mkdir $directory
fi

for min_coverage in 1 5 10 20 30
do
	out_name=sh_mincov_${min_coverage}
	out_directory=$directory/$out_name
	echo "Case $out_name"
	run-snvphyl.py --snvphyl-version 0.3 --run-name ${out_name} --galaxy-url $galaxy_url --galaxy-api-key $galaxy_api_key --fastq-history-name sh-input --reference-file reference/S_HeidelbergSL476.fasta --min-coverage $min_coverage --output-dir ${out_directory} 2>&1 | tee $directory/${out_name}.log
done
