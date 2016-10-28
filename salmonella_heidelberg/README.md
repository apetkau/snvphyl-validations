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

3. Downsample all sequence reads using `seqtk` (1.1-r92-dirty) so minimum coverage is ~30.

    ```
    for i in fastqs/*.fastq; do b=`basename $i`; echo $i; seqtk sample -s 121 $i 0.5 > fastqs-downsampled/$b; done
    ```

4. Re-check coverages.

    ```
    (for i in fastqs-downsampled/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs-downsampled/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs-downsampled/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee coverages-downsampled.txt
    ```
# Experiments

## Minimum Coverage

Run with minimum coverage of 5, 10, 15, 20.  Using <https://github.com/phac-nml/snvphyl-galaxy-cli> commit `dcebd40a3ddb5b335caa3e4de6ddc6cdd88f7a8b`.

```
dir=cov
mkdir experiments/$dir
for cov in 5 10 15 20; do name=cov-${cov}; echo $name; snvphyl.py --deploy-docker --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-downsampled --min-coverage $cov --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee min-coverage.log

# Construct files storing titles
for cov in 5 10 15 20; do echo "Minimum Coverage $cov" > experiments/cov/cov-${cov}/title; done

# Re-name directory to 05 so it sorts properly for R script
mv experiments/cov/cov-5 experiments/cov/cov-05
```

## Alternative Allele Ratio

```
dir=alt
mkdir experiments/$dir
for alt in 0.25 0.5 0.75 0.9; do name=alt-${alt}; echo $name; snvphyl.py --deploy-docker --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-downsampled/ --min-coverage 10 --alternative-allele-ratio $alt --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee alt-allele-ratio.log

for alt in 0.25 0.5 0.75 0.9; do echo "SNV Abundance Ratio $alt" > experiments/alt/alt-${alt}/title; done
```

## Sample Coverage

Sample `SH13-001` is at a coverage of ~71 after the initial downsampling.  Downsample further so that it is at coverages **10,15,20,25,30**.

```
mkdir fastqs-sample-coverage
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.423 > fastqs-sample-coverage/SH13-001_c30_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.423 > fastqs-sample-coverage/SH13-001_c30_2.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.282 > fastqs-sample-coverage/SH13-001_c20_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.282 > fastqs-sample-coverage/SH13-001_c20_2.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.212 > fastqs-sample-coverage/SH13-001_c15_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.212 > fastqs-sample-coverage/SH13-001_c15_2.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.141 > fastqs-sample-coverage/SH13-001_c10_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.141 > fastqs-sample-coverage/SH13-001_c10_2.fastq

(for i in fastqs-sample-coverage/*_1.fastq; do name=`basename $i _1.fastq`; forward=`sed -n 2~4p fastqs-sample-coverage/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p fastqs-sample-coverage/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 5,5n | tee fastqs-sample-coverage/coverages.txt
```

Make directories for each cases fastq files and link up appropriate files.

```
mkdir fastqs-sample-coverage/{c30,c20,c15,c10}
pushd fastqs-sample-coverage/c30; ln -s ../../fastqs-downsampled/*.fastq .; popd
pushd fastqs-sample-coverage/c20; ln -s ../../fastqs-downsampled/*.fastq .; popd
pushd fastqs-sample-coverage/c15; ln -s ../../fastqs-downsampled/*.fastq .; popd
pushd fastqs-sample-coverage/c10; ln -s ../../fastqs-downsampled/*.fastq .; popd
rm fastqs-sample-coverage/c*/SH13-001*.fastq

pushd fastqs-sample-coverage/c30; ln -s ../SH13-001_c30*.fastq .; popd
pushd fastqs-sample-coverage/c20; ln -s ../SH13-001_c20*.fastq .; popd
pushd fastqs-sample-coverage/c15; ln -s ../SH13-001_c15*.fastq .; popd
pushd fastqs-sample-coverage/c10; ln -s ../SH13-001_c10*.fastq .; popd
prename 's/_c\d\d//' fastqs-sample-coverage/c*/SH13-001*.fastq
```

Run SNVPhyl on each case using default parameters.

```
dir=scov
mkdir experiments/$dir
for scov in c30 c20 c15 c10; do name=scov-${scov}; echo $name; snvphyl.py --deploy-docker --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-sample-coverage/${scov} --min-coverage 10 --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee sample-coverage.log

for scov in 30 20 15 10; do echo "Subsample coverage $scov" > experiments/scov/scov-c${scov}/title; done
```

## Contamination

Pick two samples `SH12-001` from outbreak 1 at a coverage of 63.6 and `SH13-001` from outbreak 2 at a coverage of 70.9. Pick `SH13-001` as being the read set with the higher amount of data and downsample each fastq file such that the resulting dataset is at the appropriate ratios (e.g. 95% and 5%) of the original coverage of `SH13-001` (70.9).  Concatenate each read set to form resulting contaminated sample.

```
mkdir fastqs-contamination

# 5% mixture (67.4x for SH13-001, 3.5x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.95 > fastqs-contamination/SH13-001_05p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.95 > fastqs-contamination/SH13-001_05p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.0550 > fastqs-contamination/SH12-001_05p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.0550 > fastqs-contamination/SH12-001_05p_2.fastq

# 10% mixture (63.8x for SH13-001, 7.1x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.90 > fastqs-contamination/SH13-001_10p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.90 > fastqs-contamination/SH13-001_10p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.111 > fastqs-contamination/SH12-001_10p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.111 > fastqs-contamination/SH12-001_10p_2.fastq

# 20% mixture (56.7x for SH13-001, 14.2x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.80 > fastqs-contamination/SH13-001_20p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.80 > fastqs-contamination/SH13-001_20p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.223 > fastqs-contamination/SH12-001_20p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.223 > fastqs-contamination/SH12-001_20p_2.fastq

# 30% mixture (49.6x for SH13-001, 21.3x for SH12-001)
seqtk sample -s 121 fastqs-downsampled/SH13-001_1.fastq 0.70 > fastqs-contamination/SH13-001_30p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH13-001_2.fastq 0.70 > fastqs-contamination/SH13-001_30p_2.fastq

seqtk sample -s 121 fastqs-downsampled/SH12-001_1.fastq 0.334 > fastqs-contamination/SH12-001_30p_1.fastq
seqtk sample -s 121 fastqs-downsampled/SH12-001_2.fastq 0.334 > fastqs-contamination/SH12-001_30p_2.fastq
```

Create directories with other isolates and concatenate the contaminated isolates together.

```
mkdir fastqs-contamination/{30,20,10,05}p
for p in 30p 20p 10p 05p; do pushd fastqs-contamination/$p; ln -s ../../fastqs-downsampled/*.fastq .; popd; done
rm fastqs-contamination/*p/SH13-001*.fastq

# Concatenate files
for p in 30p 20p 10p 05p; do cat fastqs-contamination/*${p}_1.fastq > fastqs-contamination/$p/SH13-001_1.fastq; done
for p in 30p 20p 10p 05p; do cat fastqs-contamination/*${p}_2.fastq > fastqs-contamination/$p/SH13-001_2.fastq; done

# Check coverage of concatenated files
(for i in fastqs-contamination/*p/SH13-001*_1.fastq; do name=`basename $i _1.fastq`; dname=`dirname $i`; forward=`sed -n 2~4p $dname/${name}_1.fastq|tr -d '\n'|wc -c`; reverse=`sed -n 2~4p $dname/${name}_2.fastq|tr -d '\n'|wc -c`; ref=`bp_seq_length reference/S_HeidelbergSL476.fasta | cut -d ' ' -f 2| tr -d '\n'`; cov=`echo "($forward+$reverse)/$ref"|bc -l`; echo -e "$dname\t$name\t$forward\t$reverse\t$ref\t$cov"; done) | sort -k 6,6n | tee fastqs-contamination/coverages.txt
```

Run SNVPhyl on each case.

```
dir=contamination
mkdir experiments/$dir
for case in 30p 20p 10p 05p; do name=contamination-${case}; echo $name; snvphyl.py --deploy-docker --reference-file reference/S_HeidelbergSL476.fasta --fastq-dir fastqs-contamination/${case} --run-name $name --output-dir experiments/$dir/$name; done 2>&1 | tee contamination.log

# Titles
for p in 30 20 10 05; do echo "${p}% contaminated" > experiments/contamination/contamination-${p}p/title; done
```

# Compile Results

Please run the script `plot_trees.R`.

```
R CMD BATCH plot_trees.R
```

This will generate a file `figure3_trees.pdf` which has the completed figure.
