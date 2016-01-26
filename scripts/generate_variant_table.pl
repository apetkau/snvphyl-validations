#!/usr/bin/env perl
# Purpose
# Given a reference fasta file, number of genomes to generate, and number of positions generates a 'template' variant table to be filled in manually.

use warnings;
use strict;
use FindBin;

use lib $FindBin::Bin.'/lib';

use Bio::SeqIO;
use File::Basename;
use Getopt::Long;
use List::Util qw(shuffle);
use InvalidPositions;

my $usage =
"$0 --reference [reference.fasta] --num-genomes [number of genomes] --num-variants [number of positions] --random-seed [random seed] [--exclude-positions positions.tsv]\n".
"Parameters:\n".
"\t--reference:  A reference genome in FASTA format\n".
"\t--num-variant-genomes: Number of genomes with variation to include in the table\n".
"\t--num-duplicate-reference-genomes: Number of genomes identical to the reference to include in the table\n".
"\t--num-substitutions: Number of substitution positions to generate in table\n".
"\t--num-insertions: Number of insertion positions to generate in table\n".
"\t--num-deletions: Number of deletion positions to generate in table\n".
"\t--random-seed: Random seed for generating mutations\n".
"\t--exclude-positions: A file of positions to exclude when generating random variants (if no parameter is passed select from all positions on reference).\n".
"Example:\n".
"$0 --reference reference.fasta --num-variant-genomes 10 --num-duplicate-reference-genomes 1 --num-substituions 9000 --num-insertions 500 --num-deletions 500 --random-seed 42 --exclude-positions repeats.tsv | sort -k 1,1 -k 2,2n > variants_table.tsv\n";

my $swap_table = {
'A' => ['T','G','C'],
'T' => ['A','G','C'],
'G' => ['A','T','C'],
'C' => ['A','G','T']};

my $reference_table;
my $reference_name;
my @sequence_names;
my $number_sequences;
my $positions_used = {};

sub get_mutated_base
{
	my ($base) = @_;
	return [shuffle @{$swap_table->{uc($base)}}]->[0];
}

# reads all reference sequences into a table structured like
# ref_id => ref_seq
sub read_reference_sequences
{
	my ($reference_file) = @_;
	my %sequence_table;

	my $ref_io = Bio::SeqIO->new(-file=>"< $reference_file",-format=>"fasta");
	die "could not parse reference file $reference_file\n$usage" if (not defined $ref_io);

	while (my $seq = $ref_io->next_seq)
	{
		$sequence_table{$seq->display_id} = $seq;
	}

	return \%sequence_table;
}

sub print_header_line
{
	my ($num_duplicate_reference_genomes,$num_variant_genomes,$reference_name) = @_;

	# print header line
	print "#Chromosome\tPosition\tStatus\tReference";
	# for each duplicate reference genome
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print "\t$reference_name-duplicate-$i";
	}
	# for each variant reference genome
	for (my $i = 0; $i < $num_variant_genomes; $i++)
	{
		print "\t$reference_name-variant-$i";
	}
	print "\n";
}

sub print_reference
{
	my ($base) = @_;

	print "\t$base";
}

sub print_substitution
{
	my ($base) = @_;

	print "\t".get_mutated_base($base);
}

sub print_deletion
{
	print "\t-";
}

sub print_insertion
{
	my ($base) = @_;

	my $newbase = get_mutated_base($base);
	print "\t${base}${newbase}";
}

sub get_unique_position
{
	my ($sequence_name,$pos,$sequence,$length_sequence);

	# generate unique positions
	do 
	{
		# select random sequence
		my $seq_num = int(rand($number_sequences));

		$sequence_name = $sequence_names[$seq_num];
		$sequence = $reference_table->{$sequence_name};
		$length_sequence = $sequence->length;

		# select random position
		$pos = int(rand($length_sequence));
	} while (exists $positions_used->{"${sequence_name}_${pos}"});
	$positions_used->{"${sequence_name}_${pos}"} = 1;

	my $seq_string = $sequence->seq;
	my $ref_base = substr($seq_string,$pos,1);

	return ($sequence_name,$pos,$ref_base);
}

############
### MAIN ###
############
my ($ref_file,$num_variant_genomes,$num_duplicate_reference_genomes,$num_substitutions,$num_deletions,$num_insertions,$random_seed, $excluded_positions_file);

if (!GetOptions('reference=s' => \$ref_file,
                'num-substitutions=i' => \$num_substitutions,
                'num-insertions=i' => \$num_insertions,
                'num-deletions=i' => \$num_deletions,
                'num-variant-genomes=i' => \$num_variant_genomes,
		'num-duplicate-reference-genomes=i' => \$num_duplicate_reference_genomes,
                'random-seed=i' => \$random_seed,
		'exclude-positions=s' => \$excluded_positions_file)) {
        die "Invalid option\n".$usage;
}

die "reference.fasta is not defined\n$usage" if (not defined $ref_file);
die "$ref_file does not exist\n$usage" if (not -e $ref_file);
die "number of variant genomes is not defined\n$usage" if (not defined $num_variant_genomes);
die "number of variant genomes=$num_variant_genomes is not valid\n$usage" if ($num_variant_genomes !~ /^\d+$/);
die "number of duplicate reference genomes is not defined\n$usage" if (not defined $num_duplicate_reference_genomes);
die "number of duplicate reference genomes=$num_duplicate_reference_genomes is not valid\n$usage" if ($num_duplicate_reference_genomes !~ /^\d+$/);
die "number of substitutions=$num_substitutions is not valid\n$usage" if ($num_substitutions !~ /^\d+$/);
die "number of insertions=$num_insertions is not valid\n$usage" if ($num_insertions !~ /^\d+$/);
die "number of deletions=$num_deletions is not valid\n$usage" if ($num_deletions !~ /^\d+$/);

if (not defined $random_seed) {
	$random_seed = 42;
	warn "--random-seed not defined, defaulting to $random_seed\n";
}


srand($random_seed);


if (defined $excluded_positions_file) {
	print STDERR "Will exclude all positions in $excluded_positions_file\n";
	my $invalid_positions_parser = InvalidPositions->new;
	$positions_used = $invalid_positions_parser->read_invalid_positions($excluded_positions_file);
}


# read original reference sequences
$reference_table = read_reference_sequences($ref_file);
$reference_name = basename($ref_file, '.fasta');
@sequence_names = (keys %$reference_table);
$number_sequences = scalar(@sequence_names);

print_header_line($num_duplicate_reference_genomes,$num_variant_genomes,$reference_name);

# substitution positions
for (my $pos_num = 0; $pos_num < $num_substitutions; $pos_num++)
{
	my ($sequence_name,$pos,$ref_base) = get_unique_position();

	# print variant line, +1 to position since positions start with 1, not 0
	print "$sequence_name\t".($pos+1)."\tvalid\t$ref_base";

	# for each duplicate reference genome to generate
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print_reference($ref_base);
	}

	# for each genome to generate
	for (my $i = 0; $i < $num_variant_genomes; $i++)
	{
		print_substitution($ref_base);
	}
	print "\n";
}

# deletions
for (my $pos_num = 0; $pos_num < $num_deletions; $pos_num++)
{
	my ($sequence_name,$pos,$ref_base) = get_unique_position();

	# print variant line, +1 to position since positions start with 1, not 0
	print "$sequence_name\t".($pos+1)."\tdeletion\t$ref_base";

	# for each duplicate reference genome to generate
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print_reference($ref_base);
	}

	# print 1st variant genome with deletion, others with substitutions
	print_deletion();

	# for each genome to generate
	for (my $i = 1; $i < $num_variant_genomes; $i++)
	{
		print_substitution($ref_base);
	}
	print "\n";
}

# insertions
for (my $pos_num = 0; $pos_num < $num_insertions; $pos_num++)
{
	my ($sequence_name,$pos,$ref_base) = get_unique_position();

	# print variant line, +1 to position since positions start with 1, not 0
	print "$sequence_name\t".($pos+1)."\tinsertion\t$ref_base";

	# for each duplicate reference genome to generate
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print_reference($ref_base);
	}

	# print 1st variant genome with an insertion, others with substitutions
	print_insertions($ref_base);

	for (my $i = 1; $i < $num_variant_genomes; $i++)
	{
		print_substitution($ref_base);
	}
	print "\n";
}
