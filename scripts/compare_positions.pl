#!/usr/bin/env perl
# Purpose
# Compares two position tables of variants with respect to TP/FP/TN/FN.

use strict;
use warnings;

use Bio::SeqIO;
use Set::Scalar;
use Getopt::Long;

my $false_detection_fh;

sub is_substitution
{
	return ($_[0] =~ /substitution$/);
}

sub is_insertion
{
	return ($_[0] =~ /insertion$/);
}

sub is_deletion
{
	return ($_[0] =~ /deletion$/);
}

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
sub read_pos_true_variants_table
{
	my ($file, $position_set_to_remove) = @_;

	$position_set_to_remove = Set::Scalar->new if (not defined $position_set_to_remove);

	my %results;

	open(my $fh,"<$file") or die "Could not open $file\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	$results{'header'} = $header;

	# Sets for columns in an alignment in all regions with the appropriate status
	$results{'columns-all'}{'all'} = Set::Scalar->new;
	$results{'columns-all'}{'substitutions'} = Set::Scalar->new;
	$results{'columns-all'}{'deletions'} = Set::Scalar->new;
	$results{'columns-all'}{'insertions'} = Set::Scalar->new;

	# Sets for columns in an alignment in non-repeat regions with the appropriate status
	$results{'columns-valid'}{'all'} = Set::Scalar->new;
	$results{'columns-valid'}{'substitutions'} = Set::Scalar->new;
	$results{'columns-valid'}{'deletions'} = Set::Scalar->new;
	$results{'columns-valid'}{'insertions'} = Set::Scalar->new;

	$results{'positions-all'}{'substitutions'} = Set::Scalar->new;
	$results{'positions-all'}{'insertions'} = Set::Scalar->new;
	$results{'positions-all'}{'deletions'} = Set::Scalar->new;

	$results{'positions-valid'}{'substitutions'} = Set::Scalar->new;
	$results{'positions-valid'}{'insertions'} = Set::Scalar->new;
	$results{'positions-valid'}{'deletions'} = Set::Scalar->new;

	while(my $line = readline($fh))
	{
		my ($chrom,$position,$status,$line_minus_status) = parse_pos_line($line);
		my $chrom_pos = "$chrom\t$position";

		$results{'columns-all'}{'all'}->insert($line_minus_status);
		$results{'columns-all'}{'substitutions'}->insert($line_minus_status) if (is_substitution($status));
		$results{'columns-all'}{'deletions'}->insert($line_minus_status) if (is_deletion($status));
		$results{'columns-all'}{'insertions'}->insert($line_minus_status) if (is_insertion($status));

		$results{'positions-all'}{'substitutions'}->insert($chrom_pos) if (is_substitution($status));
		$results{'positions-all'}{'deletions'}->insert($chrom_pos) if (is_deletion($status));
		$results{'positions-all'}{'insertions'}->insert($chrom_pos) if (is_insertion($status));
		
		if ($status =~ /^valid/)
		{
			$results{'columns-valid'}{'all'}->insert($line_minus_status);
			$results{'columns-valid'}{'substitutions'}->insert($line_minus_status) if (is_substitution($status));
			$results{'columns-valid'}{'deletions'}->insert($line_minus_status) if (is_deletion($status));
			$results{'columns-valid'}{'insertions'}->insert($line_minus_status) if (is_insertion($status));

			$results{'positions-valid'}{'substitutions'}->insert($chrom_pos) if (is_substitution($status));
			$results{'positions-valid'}{'deletions'}->insert($chrom_pos) if (is_deletion($status));
			$results{'positions-valid'}{'insertions'}->insert($chrom_pos) if (is_insertion($status));
		}
	}
	close($fh);

	return \%results;
}

sub intersect_columns_by_position
{
	my ($columns_a, $columns_b) = @_;

	my $intersection = Set::Scalar->new;

	# generate a table which the combined chromosome/position acts as a key to point to the column line
	my %columns_a_table = ();
	for my $e ($columns_a->elements)
	{
		my ($chrom,$pos) = split(/\t/, $e);
		$columns_a_table{"$chrom\t$pos"} = $e;
	}
	my %columns_b_table = ();
	for my $e ($columns_b->elements)
	{
		my ($chrom,$pos) = split(/\t/, $e);
		$columns_b_table{"$chrom\t$pos"} = $e;
	}

	# for every position in table a, add full line to intersection if it exists in table b
	for my $k (keys %columns_a_table)
	{
		$intersection->insert($columns_b_table{$k}) if (exists $columns_b_table{$k});
	}

	return $intersection;
}

sub parse_pos_line
{
	my ($line) = @_;

	chomp($line);
	my @tokens = split(/\t/,$line);

	my ($chrom,$position,$status,@bases) = @tokens;
	die "Error with line $line\n" if (not defined $chrom or not defined $position or $position !~ /\d+/);
	die "Error with line $line, status not properly defined\n" if (not defined $status or $status eq '');
	my $line_minus_status = join("\t",$chrom,$position,@bases);

	return ($chrom,$position,$status,$line_minus_status);
}

sub find_differences_column_with_base_call
{
	my ($base_call, $column_set) = @_;

	my $results = Set::Scalar->new;

	for my $e ($column_set->elements)
	{
		my ($chrom,$pos,$status,@bases) = split(/\t/,$e);
		$results->insert($e) if ((grep {/^$base_call$/} @bases) >= 1);
	}

	return $results;
}

sub read_pos_actual_variants_table
{
	my ($file, $substitutions_pos, $insertions_pos, $deletions_pos) = @_;

	my %results;

	open(my $fh,"<$file") or die "Could not open $file\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	$results{'header'} = $header;
	$results{'columns-valid'}{'all'} = Set::Scalar->new;
	$results{'columns-valid'}{'substitutions'} = Set::Scalar->new;
	$results{'columns-valid'}{'insertions'} = Set::Scalar->new;
	$results{'columns-valid'}{'deletions'} = Set::Scalar->new;
	$results{'columns-invalid'}{'all'} = Set::Scalar->new;
	$results{'columns-all'}{'all'} = Set::Scalar->new;
	$results{'columns-all'}{'substitutions'} = Set::Scalar->new;
	$results{'columns-all'}{'insertions'} = Set::Scalar->new;
	$results{'columns-all'}{deletions} = Set::Scalar->new;
	while(my $line = readline($fh))
	{
		my ($chrom,$position,$status,$line_minus_status) = parse_pos_line($line);
		my $chrom_pos = "$chrom\t$position";

		$results{'columns-all'}{'all'}->insert($line_minus_status);
		$results{'columns-all'}{'substitutions'}->insert($line_minus_status) if ($substitutions_pos->has($chrom_pos));
		$results{'columns-all'}{'insertions'}->insert($line_minus_status) if ($insertions_pos->has($chrom_pos));
		$results{'columns-all'}{'deletions'}->insert($line_minus_status) if ($deletions_pos->has($chrom_pos));

		if ($status eq 'valid')
		{
			$results{'columns-valid'}{'all'}->insert($line_minus_status);
			$results{'columns-valid'}{'substitutions'}->insert($line_minus_status) if ($substitutions_pos->has($chrom_pos));
			$results{'columns-valid'}{'insertions'}->insert($line_minus_status) if ($insertions_pos->has($chrom_pos));
			$results{'columns-valid'}{'deletions'}->insert($line_minus_status) if ($deletions_pos->has($chrom_pos));
		}
		else
		{
			$results{'columns-invalid'}{'all'}->insert($line_minus_status);
		}
	}
	close($fh);

	return \%results;
}

sub get_comparisons
{
	my ($var_true_col,$var_detected_col,$true_nonvariant_columns, $case) = @_;

	# set operations
	my $true_positives_set = $var_true_col->{'all'} * $var_detected_col->{'all'};
	my $false_positives_set = $var_detected_col->{'all'} - $var_true_col->{'all'};
	my $false_negatives_set = $var_true_col->{'all'} - $var_detected_col->{'all'};
	
	my $true_col_positives = $var_true_col->{'all'}->size;
	my $detected_col_positives = $var_detected_col->{'all'}->size;
	my $true_positives = $true_positives_set->size;
	my $false_positives = $false_positives_set->size;

	my $true_negatives = $true_nonvariant_columns - $var_detected_col->{'all'}->size;

	my $true_col_positives_sub = intersect_columns_by_position($var_true_col->{'all'}, $var_true_col->{'substitutions'})->size;
	my $true_col_positives_ins = intersect_columns_by_position($var_true_col->{'all'}, $var_true_col->{'insertions'})->size;
	my $true_col_positives_del = intersect_columns_by_position($var_true_col->{'all'}, $var_true_col->{'deletions'})->size;

	my $detected_col_positives_sub = intersect_columns_by_position($var_detected_col->{'all'}, $var_detected_col->{'substitutions'})->size;
	my $detected_col_positives_ins = intersect_columns_by_position($var_detected_col->{'all'}, $var_detected_col->{'insertions'})->size;
	my $detected_col_positives_del = intersect_columns_by_position($var_detected_col->{'all'}, $var_detected_col->{'deletions'})->size;
	my $detected_col_positives_other = $detected_col_positives - ($detected_col_positives_sub + $detected_col_positives_ins + $detected_col_positives_del);

	my $true_positives_sub = intersect_columns_by_position($true_positives_set, $var_true_col->{'substitutions'})->size;
	my $true_positives_ins = intersect_columns_by_position($true_positives_set, $var_true_col->{'insertions'})->size;
	my $true_positives_del = intersect_columns_by_position($true_positives_set, $var_true_col->{'deletions'})->size;

	my $false_positives_sub_set = intersect_columns_by_position($false_positives_set, $var_detected_col->{'substitutions'});
	my $false_positives_sub = $false_positives_sub_set->size;
	my $false_positives_ins = intersect_columns_by_position($false_positives_set, $var_detected_col->{'insertions'})->size;
	my $false_positives_del = intersect_columns_by_position($false_positives_set, $var_detected_col->{'deletions'})->size;
	my $false_positives_other = $false_positives - ($false_positives_sub+$false_positives_ins+$false_positives_del);

	my $false_negatives_sub_set = intersect_columns_by_position($false_negatives_set, $var_true_col->{'substitutions'});
	my $false_negatives_sub = $false_negatives_sub_set->size;
	my $false_negatives_ins = intersect_columns_by_position($false_negatives_set, $var_true_col->{'insertions'})->size;
	my $false_negatives_del = intersect_columns_by_position($false_negatives_set, $var_true_col->{'deletions'})->size;
	my $false_negatives = $false_negatives_set->size;

	my $false_positives_filtered_base_set = find_differences_column_with_base_call('N',$false_positives_sub_set) + find_differences_column_with_base_call('-',$false_positives_sub_set);
	print $false_detection_fh "$case\tFP\t",join("\n$case\tFP\t",($false_positives_sub_set - $false_positives_filtered_base_set)->elements),"\n";

	my $false_negatives_filtered_base_set = find_differences_column_with_base_call('N',$false_negatives_sub_set) + find_differences_column_with_base_call('-',$false_negatives_sub_set);
	print $false_detection_fh "$case\tFN\t",join("\n$case\tFN\t",($false_negatives_sub_set - $false_negatives_filtered_base_set)->elements),"\n";
	
	my $accuracy = sprintf "%0.4f",($true_positives + $true_negatives) / ($true_positives + $false_positives + $true_negatives + $false_negatives);
	my $specificity = sprintf "%0.4f",($true_negatives) / ($true_negatives + $false_positives);
	my $sensitivity = sprintf "%0.4f",($true_positives) / ($true_positives + $false_negatives);
	my $precision = sprintf "%0.4f",($true_positives) / ($true_positives + $false_positives);
	my $fp_rate = sprintf "%0.4f",($false_positives) / ($true_negatives + $false_positives);
	
	return "$true_col_positives\t$true_nonvariant_columns\t$detected_col_positives\t$true_positives\t$false_positives\t$true_negatives\t$false_negatives\t".
		"$accuracy\t$specificity\t$sensitivity\t$precision\t$fp_rate\n";
}

my $usage = "$0 --variants-true [variants-true.tsv] --variants-detected [variants-detected.tsv] --reference-genome [reference-genome.fasta] --false-detection-output [false detection output]\n".
"Parameters:\n".
"\t--variants-true: The true variants table.\n".
"\t--variants-detected: The detected variants table\n".
"\t--reference-genome: The reference genome in fasta format.  This is used to get the length to calculate the false negative rate.\n".
"\t--false-detection-output: Output file for storing columns with bases that were miscalled\n".
"Example:\n".
"$0 --variants-true variants.tsv --variants-detected variants-detected.tsv --reference-genome reference.fasta\n\n";

my ($variants_true_file,$variants_detected_file, $reference_genome_file, $false_detection_output);

if (!GetOptions('variants-true=s' => \$variants_true_file,
		'variants-detected=s' => \$variants_detected_file,
		'false-detection-output=s' => \$false_detection_output,
		'reference-genome=s' => \$reference_genome_file))
{
	die "Invalid option\n".$usage;
}

die "--variants-true not defined\n$usage" if (not defined $variants_true_file);
die "--variants-detected not defined\n$usage" if (not defined $variants_detected_file);
die "--reference-genome not defined\n$usage" if (not defined $reference_genome_file);
die "--false-detection-output not defined\n$usage" if (not defined $false_detection_output);

my $reference_genome_obj = Bio::SeqIO->new(-file=>"<$reference_genome_file", -format=>"fasta");
my $reference_genome_size = 0;
while (my $seq = $reference_genome_obj->next_seq) {
	$reference_genome_size += $seq->length;
}

my $variants_true = read_pos_true_variants_table($variants_true_file);
my $variants_detected = read_pos_actual_variants_table($variants_detected_file,
	$variants_true->{'positions-all'}{'substitutions'}, $variants_true->{'positions-all'}{'insertions'}, $variants_true->{'positions-all'}{'deletions'});


# must have same genomes and same order of genomes
die "Error: headers did not match\n" if ($variants_true->{'header'} ne $variants_detected->{'header'});

my $true_nonvariant_columns_all = $reference_genome_size - $variants_true->{'columns-all'}{'all'}->size;
my $true_nonvariant_columns_valid = $reference_genome_size - $variants_true->{'columns-valid'}{'all'}->size;

open($false_detection_fh, ">$false_detection_output") or die "Could not open $false_detection_output for writing";

print "Reference_Genome_File\t$reference_genome_file\n";
print "Reference_Genome_Size\t$reference_genome_size\n";
print "Variants_True_File\t$variants_true_file\n";
print "Variants_Detected_File\t$variants_detected_file\n";

print "Case\tTrue_Variant_Columns\tTrue_Nonvariant_Columns\tColumns_Detected\tTP\tFP\tTN\tFN\tAccuracy\tSpecificity\tSensitivity\tPrecision\tFPR\n";
print "all-vs-valid\t",get_comparisons($variants_true->{'columns-all'}, $variants_detected->{'columns-valid'}, $true_nonvariant_columns_all, 'all-vs-valid');
print "valid-vs-valid\t",get_comparisons($variants_true->{'columns-valid'}, $variants_detected->{'columns-valid'}, $true_nonvariant_columns_valid, 'valid-vs-valid');
print "all-vs-all\t",get_comparisons($variants_true->{'columns-all'}, $variants_detected->{'columns-all'}, $true_nonvariant_columns_all, 'all-vs-all');

close($false_detection_fh);
