# Initial Setup

Directory `S_heildelberg_QC/` should contain a folder for each sample, with fastq files within each folder.

# Data Preparation

1. Re-name all files/link into `fastqs/` directory.

    ```
    cd fastqs/
    for i in `cat ../../../data/strains.txt`; do prev=`echo $i|cut -d ',' -f 1`; curr=`echo $i|cut -d ',' -f 2`; ln -s ../S_heildelberg_QC/$prev/*_R1_001.fastq ${curr}_1.fastq; done
    for i in `cat ../../../data/strains.txt`; do prev=`echo $i|cut -d ',' -f 1`; curr=`echo $i|cut -d ',' -f 2`; ln -s ../S_heildelberg_QC/$prev/*_R2_001.fastq ${curr}_2.fastq; done
    cd ..
    ```

2. Estimate coverages for each sample based on reference genome length.

    ```
    (for i in fastqs/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee coverages.txt
    ```

Minimum coverage is `61`, so downsample accordingly.

3. Downsample all sequence reads using `seqtk` (1.0-r31) so minimum coverage is ~30.

    ```
    for i in fastqs/*.fastq; do b=`basename $i`; echo $i; seqtk sample -s 121 $i 0.5 > fastqs-downsampled/$b; done
    ```

4. Re-check coverages.

    ```
    (for i in fastqs-downsampled/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs-downsampled/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs-downsampled/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee coverages-downsampled.txt
    ```

5. Run SNVPhyl with all default settings initially.  Will upload all fastq files to Galaxy.

    ```
    run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-downsampled/ --run-name initial-snvphyl-run --output-dir experiments/initial-snvphyl-run
    ```

# Experiments

## Minimum Coverage

Run with minimum coverage of 5, 10, 20 (case 15 was run initially).

```
dir=cov
mkdir experiments/$dir
for cov in 5 10 20; do name=cov-${cov}; echo $name; run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-history-name 'snvphyl-S_HeidelbergSL476-2016-02-07-initial-snvphyl-run' --min-coverage $cov --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee min-coverage.log
```

## Alternative Allele Ratio

```
dir=alt
mkdir experiments/$dir
for alt in 0.25 0.5 0.9; do name=alt-${alt}; echo $name; run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-history-name 'snvphyl-S_HeidelbergSL476-2016-02-07-initial-snvphyl-run' --alternative-allele-ratio $alt --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee alt-allele-ratio.log
```

## Sample Coverage

Sample `SH13-001` is at a coverage of ~71 after the initial downsampling.  Downsample further so that it is at coverages **15,20,30**.

```
mkdir fastqs-sample-coverage
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.42 > fastqs-sample-coverage/SH13-001_c30_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.42 > fastqs-sample-coverage/SH13-001_c30_2.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.28 > fastqs-sample-coverage/SH13-001_c20_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.28 > fastqs-sample-coverage/SH13-001_c20_2.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.21 > fastqs-sample-coverage/SH13-001_c15_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.21 > fastqs-sample-coverage/SH13-001_c15_2.fastq

(for i in fastqs-sample-coverage/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs-sample-coverage/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs-sample-coverage/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee fastqs-sample-coverage/coverages.txt
```

Manually upload fastq files to Galaxy and combine with other samples to construct datasets.

## Contamination

# Compile Results


