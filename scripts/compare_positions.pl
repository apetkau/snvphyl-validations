#!/usr/bin/env perl
# Purpose
# Compares two position tables of variants with respect to TP/FP/TN/FN.

use strict;
use warnings;
use FindBin;

use lib $FindBin::Bin.'/lib';

use Bio::SeqIO;
use Set::Scalar;
use Getopt::Long;
use InvalidPositions;

# read positions file a Set::Scalar object
# format
# {
#	'header' => header_line,
# 	'positions-valid' => Set of position lines with 'valid' status
#	'positions-invalid' => Set of position lines with 'invalid' status
# }
sub read_pos
{
	my ($file) = @_;

	my %results;

	open(my $fh,"<$file") or die "Could not open $file\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	$results{'header'} = $header;
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
		if ($status eq 'valid') {
			$results{'positions-valid'}->insert($line_minus_status);
		} else {
			$results{'positions-invalid'}->insert($line_minus_status);
		}
	}
	close($fh);

	return \%results;
}

my $usage = "$0 --variants-true [variants-true.tsv] --variants-detected [variants-detected.tsv] --core-genome [core-genome.fasta] [--exclude-positions positions.tsv]\n".
"Parameters:\n".
"\t--variants-true: The true variants table.\n".
"\t--variants-detected: The detected variants table\n".
"\t--core-genome: The core genome in fasta format.  For most cases this is assumed to be the same as the reference genome\n".
"\t--exclude-positions: Optional positions to exclude from counts of the core genome.  Assumes variant table was run with --exclude-positions\n".
"Example:\n".
"$0 --variants-true variants.tsv --variants-detected variants-detected.tsv --core-genome reference.fasta --exclude-positions repeats.tsv\n\n";

my ($variants_true_file,$variants_detected_file, $core_genome_file, $excluded_positions_file);

if (!GetOptions('variants-true=s' => \$variants_true_file,
		'variants-detected=s' => \$variants_detected_file,
		'core-genome=s' => \$core_genome_file,
		'exclude-positions=s' => \$excluded_positions_file))
{
	die "Invalid option\n".$usage;
}

die "--variants-true not defined\n$usage" if (not defined $variants_true_file);
die "--variants-detected not defined\n$usage" if (not defined $variants_detected_file);
die "--core-genome not defined\n$usage" if (not defined $core_genome_file);

my $core_genome_obj = Bio::SeqIO->new(-file=>"<$core_genome_file", -format=>"fasta");
my $core_genome_size = 0;
while (my $seq = $core_genome_obj->next_seq) {
	$core_genome_size += $seq->length;
}

my $core_genome_size_minus_excluded_positions = $core_genome_size;
if (defined $excluded_positions_file) {
	my $invalid_pos_parser = InvalidPositions->new;
	my $invalid_positions = $invalid_pos_parser->read_invalid_positions($excluded_positions_file);
	my $count_excluded_positions = scalar(keys %$invalid_positions);
	print STDERR "Excluded $count_excluded_positions defined in file $excluded_positions_file\n";
	$core_genome_size_minus_excluded_positions -= $count_excluded_positions;
}

my $variants_true = read_pos($variants_true_file);
my $variants_detected = read_pos($variants_detected_file);

# must have same genomes and same order of genomes
if ($variants_true->{'header'} ne $variants_detected->{'header'})
{
	die "Error: headers did not match\n";
}

my $var_true_valid_pos = $variants_true->{'positions-valid'};
my $var_detected_valid_pos = $variants_detected->{'positions-valid'};

# Count of non variant positions in core genome
my $non_variant_true_positions_count = $core_genome_size_minus_excluded_positions - $var_true_valid_pos->size;

# set operations
my $true_positives_set = $var_true_valid_pos * $var_detected_valid_pos;
my $false_positives_set = $var_detected_valid_pos - $var_true_valid_pos;
# See comment for true negatives below
#my $true_negatives_set;
my $false_negatives_set = $var_true_valid_pos - $var_detected_valid_pos;

my $true_valid_positives = $var_true_valid_pos->size;
my $detected_valid_positives = $var_detected_valid_pos->size;
my $true_positives = $true_positives_set->size;
my $false_positives = $false_positives_set->size;

# True Negatives are positions in our alignment that have no variant in any genome 
# That is, the number of positions in the core minus the total variant positions detected.
# This assumes that our definition of core genome size above is valid.
my $true_negatives = $non_variant_true_positions_count;

my $false_negatives = $false_negatives_set->size;
my $accuracy = sprintf "%0.4f",($true_positives + $true_negatives) / ($true_positives + $false_positives + $true_negatives + $false_negatives);
my $specificity = sprintf "%0.4f",($true_negatives) / ($true_negatives + $false_positives);
my $sensitivity = sprintf "%0.4f",($true_positives) / ($true_positives + $false_negatives);
my $precision = sprintf "%0.4f",($true_positives) / ($true_positives + $false_positives);
my $fp_rate = sprintf "%0.4f",($false_positives) / ($true_negatives + $false_positives);

print "Core_Genome_File\tCore_Genome\tCore_Genome_Ex_Positions\tVariants_True_File\tVariants_Detected_File\tTrue_Variants\tVariants_Detected\tTP\tFP\tTN\tFN\tAccuracy\tSpecificity\tSensitivity\tPrecision\tFP_Rate\n";
print "$core_genome_file\t$core_genome_size\t$core_genome_size_minus_excluded_positions\t$variants_true_file\t$variants_detected_file\t$true_valid_positives\t$detected_valid_positives\t$true_positives\t$false_positives\t$true_negatives\t$false_negatives\t".
      "$accuracy\t$specificity\t$sensitivity\t$precision\t$fp_rate\n";
