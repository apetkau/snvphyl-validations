#!/usr/bin/env perl
# Purpose
# Compares two position tables of variants with respect to TP/FP/TN/FN.

use strict;
use warnings;

use Bio::SeqIO;
use Set::Scalar;
use Getopt::Long;

# read positions file to a set of Set::Scalar objects
# Input
#	$file  File to read positions from
#	$position_set_to_remove  Optional set of "$chrom\t$pos" coordinates to use when generating {'columns-valid-removed-positions'}
# Output
#	A hash table of sets of coordinates and lines (with positions/base information) in the format below.
# {
#	'header' => header_line,
# 	'columns-valid' => Set of position lines with 'valid' status
#	'columns-invalid' => Set of position lines which do not have a 'valid' status
#	'columns-all' => Set of all position lines
#	'columns-valid-removed-positions' => Set of position lines minus any lines with coordinates passed in $position_set_to_remove
#	'columns-removed-positions' => Set of position lines that would have been removed by coordinates passed in $position_set_to_remove
#	'positions-valid' => Set of position coordinates with 'valid' status
#	'positions-invalid' => Set of position coordinates which do not have a 'valid' status
# }
sub read_pos
{
	my ($file, $position_set_to_remove) = @_;

	$position_set_to_remove = Set::Scalar->new if (not defined $position_set_to_remove);

	my %results;

	open(my $fh,"<$file") or die "Could not open $file\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	$results{'header'} = $header;
	$results{'columns-valid'} = Set::Scalar->new;
	$results{'columns-valid-removed-positions'} = Set::Scalar->new;
	$results{'columns-removed-positions'} = Set::Scalar->new;
	$results{'columns-invalid'} = Set::Scalar->new;
	$results{'columns-all'} = Set::Scalar->new;
	$results{'positions-valid'} = Set::Scalar->new;
	$results{'positions-invalid'} = Set::Scalar->new;
	while(my $line = readline($fh))
	{
		chomp($line);
		my @tokens = split(/\t/,$line);

		my ($chrom,$position,$status,@bases) = @tokens;
		die "Error with line $line in $file\n" if (not defined $chrom or not defined $position or $position !~ /\d+/);
		die "Error with line $line in $file, status not properly defined\n" if (not defined $status or $status eq '');

		my $line_minus_status = join(' ',$chrom,$position,@bases);
		$results{'columns-all'}->insert($line_minus_status);
		if ($status eq 'valid') {
			$results{'columns-valid'}->insert($line_minus_status);
			$results{'positions-valid'}->insert("$chrom\t$position");

			if ($position_set_to_remove->has("$chrom\t$position"))
			{
				$results{'columns-removed-positions'}->insert($line_minus_status);
			}
			else
			{
				$results{'columns-valid-removed-positions'}->insert($line_minus_status);
			}
		}
		else
		{
			if ($position_set_to_remove->has("$chrom\t$position"))
                        {
				$results{'columns-removed-positions'}->insert($line_minus_status);
			}

			$results{'columns-invalid'}->insert($line_minus_status);
			$results{'positions-invalid'}->insert("$chrom\t$position");
		}
	}
	close($fh);

	return \%results;
}

sub get_comparisons
{
	my ($var_true_pos,$var_detected_pos,$reference_genome_size) = @_;

	# set operations
	my $true_positives_set = $var_true_pos * $var_detected_pos;
	my $false_positives_set = $var_detected_pos - $var_true_pos;
	# See comment for true negatives below
	#my $true_negatives_set;
	my $false_negatives_set = $var_true_pos - $var_detected_pos;
	
	my $true_valid_positives = $var_true_pos->size;
	my $detected_valid_positives = $var_detected_pos->size;
	my $true_positives = $true_positives_set->size;
	my $false_positives = $false_positives_set->size;
	
	# True Negatives are columns in our alignment that have no variant in any genome 
	# That is, the number of columns in the core minus the total variant columns detected.
	# This assumes that our definition of core genome size above is valid.
	my $true_negatives = $reference_genome_size - $var_detected_pos->size;
	
	my $false_negatives = $false_negatives_set->size;
	my $accuracy = sprintf "%0.4f",($true_positives + $true_negatives) / ($true_positives + $false_positives + $true_negatives + $false_negatives);
	my $specificity = sprintf "%0.4f",($true_negatives) / ($true_negatives + $false_positives);
	my $sensitivity = sprintf "%0.4f",($true_positives) / ($true_positives + $false_negatives);
	my $precision = sprintf "%0.4f",($true_positives) / ($true_positives + $false_positives);
	my $fp_rate = sprintf "%0.4f",($false_positives) / ($true_negatives + $false_positives);
	
	return "$true_valid_positives\t$detected_valid_positives\t$true_positives\t$false_positives\t$true_negatives\t$false_negatives\t".
		"$accuracy\t$specificity\t$sensitivity\t$precision\t$fp_rate";
}

my $usage = "$0 --variants-true [variants-true.tsv] --variants-detected [variants-detected.tsv] --reference-genome [reference-genome.fasta]\n".
"Parameters:\n".
"\t--variants-true: The true variants table.\n".
"\t--variants-detected: The detected variants table\n".
"\t--reference-genome: The reference genome in fasta format.  This is used to get the length to calculate the false negative rate.\n".
"Example:\n".
"$0 --variants-true variants.tsv --variants-detected variants-detected.tsv --reference-genome reference.fasta\n\n";

my ($variants_true_file,$variants_detected_file, $reference_genome_file);

if (!GetOptions('variants-true=s' => \$variants_true_file,
		'variants-detected=s' => \$variants_detected_file,
		'reference-genome=s' => \$reference_genome_file))
{
	die "Invalid option\n".$usage;
}

die "--variants-true not defined\n$usage" if (not defined $variants_true_file);
die "--variants-detected not defined\n$usage" if (not defined $variants_detected_file);
die "--reference-genome not defined\n$usage" if (not defined $reference_genome_file);

my $reference_genome_obj = Bio::SeqIO->new(-file=>"<$reference_genome_file", -format=>"fasta");
my $reference_genome_size = 0;
while (my $seq = $reference_genome_obj->next_seq) {
	$reference_genome_size += $seq->length;
}

my $variants_detected = read_pos($variants_detected_file);
my $variants_true = read_pos($variants_true_file, $variants_detected->{'positions-invalid'});

# must have same genomes and same order of genomes
if ($variants_true->{'header'} ne $variants_detected->{'header'})
{
	die "Error: headers did not match\n";
}

print "Reference_Genome_File\t$reference_genome_file\n";
print "Reference_Genome_Size\t$reference_genome_size\n";
print "Variants_True_File\t$variants_true_file\n";
print "Variants_Detected_File\t$variants_detected_file\n";
print "Case\tTrue_Columns\tColumns_Detected\tTP\tFP\tTN\tFN\tAccuracy\tSpecificity\tSensitivity\tPrecision\tFP_Rate\n";
print "all-vs-valid\t".get_comparisons($variants_true->{'columns-all'}, $variants_detected->{'columns-valid'}, $reference_genome_size)."\n";
print "valid-vs-valid\t".get_comparisons($variants_true->{'columns-valid'}, $variants_detected->{'columns-valid'}, $reference_genome_size)."\n";
print "all-minus-detected-invalid-vs-valid\t".get_comparisons($variants_true->{'columns-valid-removed-positions'}, $variants_detected->{'columns-valid'}, $reference_genome_size)."\n";
print "detected-invalid-vs-invalid\t".get_comparisons($variants_true->{'columns-removed-positions'}, $variants_detected->{'columns-invalid'}, $reference_genome_size)."\n";
print "all-vs-all\t".get_comparisons($variants_true->{'columns-all'}, $variants_detected->{'columns-all'}, $reference_genome_size)."\n";
