#!/usr/bin/env perl
# Purpose
# Compares two position tables of variants with respect to TP/FP/TN/FN.

use strict;
use warnings;

use Set::Scalar;

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

my $usage = "$0 [variants-true.tsv] [variants-detected.tsv]\n";

die "error:\n$usage" if (@ARGV != 2);

my ($variants_true_file,$variants_detected_file) = @ARGV;

my $variants_true = read_pos($variants_true_file);
my $variants_detected = read_pos($variants_detected_file);

# must have same genomes and same order of genomes
if ($variants_true->{'header'} ne $variants_detected->{'header'})
{
	die "Error: headers did not match\n";
}

my $var_true_valid_pos = $variants_true->{'positions-valid'};
my $var_detected_valid_pos = $variants_detected->{'positions-valid'};

# set operations
my $true_positives = $var_true_valid_pos * $var_detected_valid_pos;
my $false_positives = $var_detected_valid_pos - $var_true_valid_pos;
my $true_negatives = Set::Scalar->new;
my $false_negatives = $var_true_valid_pos - $var_detected_valid_pos;

print "$variants_true_file\t$variants_detected_file\tTP\tFP\tTN\tFN\n";
print $var_true_valid_pos->size,"\t",$var_detected_valid_pos->size,"\t",$true_positives->size,
	"\t",$false_positives->size,"\t",$true_negatives->size,"\t",$false_negatives->size,"\n";
