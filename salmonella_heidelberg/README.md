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

# link previous run to proper directory
cp -r experiments/initial-snvphyl-run/ experiments/cov/cov-15

# Construct files storing titles
for cov in 5 10 15 20; do echo "Minimum Coverage $cov" > experiments/cov/cov-${cov}/title; done
```

## Alternative Allele Ratio

```
dir=alt
mkdir experiments/$dir
for alt in 0.25 0.5 0.9; do name=alt-${alt}; echo $name; run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-history-name 'snvphyl-S_HeidelbergSL476-2016-02-07-initial-snvphyl-run' --alternative-allele-ratio $alt --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee alt-allele-ratio.log

# link previous run to proper directory
cp -r experiments/initial-snvphyl-run/ experiments/alt/alt-0.75

for alt in 0.25 0.5 0.75 0.9; do echo "Alt. Allele Ratio $alt" > experiments/alt/alt-${alt}/title; done
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

Make directories for each cases fastq files and link up appropriate files.

```
mkdir fastqs-sample-coverage/{c30,c20,c15}
pushd fastqs-sample-coverage/c30; ln -s ../../fastqs-downsampled/*.fastq .; popd
pushd fastqs-sample-coverage/c20; ln -s ../../fastqs-downsampled/*.fastq .; popd
pushd fastqs-sample-coverage/c15; ln -s ../../fastqs-downsampled/*.fastq .; popd
rm fastqs-sample-coverage/c*/SH13-001*.fastq

pushd fastqs-sample-coverage/c30; ln -s ../SH13-001_c30*.fastq .; popd
pushd fastqs-sample-coverage/c20; ln -s ../SH13-001_c20*.fastq .; popd
pushd fastqs-sample-coverage/c15; ln -s ../SH13-001_c15*.fastq .; popd
prename 's/_c\d\d//' fastqs-sample-coverage/c*/SH13-001*.fastq
```

Run SNVPhyl on each case using default parameters.

```
dir=scov
mkdir experiments/$dir
for scov in c30 c20 c15; do name=scov-${scov}; echo $name; echo run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-sample-coverage/${scov} --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee sample-coverage.log
```

## Contamination

Pick two samples `SH12-001` from outbreak 1 at a coverage of ~64 and `SH13-001` from outbreak 2 at a coverage of ~71. Downsample each and concatenate files such that the resulting sample is formed from a mixtured of reads at the appropriate ratios, with `SH13-001` being the higher ratio.

```
mkdir fastqs-contamination

# 5% mixture (28.5x for SH13-001, 1.5x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.402 > fastqs-contamination/SH13-001_5p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.402 > fastqs-contamination/SH13-001_5p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.0236 > fastqs-contamination/SH12-001_5p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.0236 > fastqs-contamination/SH12-001_5p_2.fastq

# 10% mixture (27x for SH13-001, 3x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.381 > fastqs-contamination/SH13-001_10p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.381 > fastqs-contamination/SH13-001_10p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.0472 > fastqs-contamination/SH12-001_10p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.0472 > fastqs-contamination/SH12-001_10p_2.fastq

# 20% mixture (24x for SH13-001, 6x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.339 > fastqs-contamination/SH13-001_20p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.339 > fastqs-contamination/SH13-001_20p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.0944 > fastqs-contamination/SH12-001_20p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.0944 > fastqs-contamination/SH12-001_20p_2.fastq

# 50% mixture (15x for SH13-001, 15x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.212 > fastqs-contamination/SH13-001_50p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.212 > fastqs-contamination/SH13-001_50p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.236 > fastqs-contamination/SH12-001_50p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.236 > fastqs-contamination/SH12-001_50p_2.fastq
```

Create directories with other isolates and concatenate the contaminated isolates together.

```
mkdir fastqs-contamination/{50,20,10,5}p
for p in 50p 20p 10p 5p; do pushd fastqs-contamination/$p; ln -s ../../fastqs-downsampled/*.fastq .; popd; done
rm fastqs-contamination/*p/SH13-001*.fastq

# Concatenate files
for p in 50p 20p 10p 5p; do cat fastqs-contamination/*${p}_1.fastq > fastqs-contamination/$p/SH13-001_1.fastq; done
for p in 50p 20p 10p 5p; do cat fastqs-contamination/*${p}_2.fastq > fastqs-contamination/$p/SH13-001_2.fastq; done

# Check coverage of concatenated files
(for i in fastqs-contamination/*p/SH13-001*_1.fastq; do name=`basename $i _1.fastq`; dname=`dirname $i`; forward=`sed -n 2~4p $dname/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p $dname/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$dname\t$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 6,6n | tee fastqs-contamination/coverages.txt
```

Run SNVPhyl on each case.

```
dir=contamination
mkdir experiments/$dir
for case in 50p 20p 10p 5p; do name=contamination-${case}; echo $name; run-snvphyl.py --galaxy-url [URL] --galaxy-api-key [KEY] --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-contamination/${case} --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee contamination.log

# Titles
for p in 50 20 10 5; do echo "${p}% contaminated" > experiments/contamination/contamination-${p}p/title; done
```

# Compile Results

