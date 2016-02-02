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
"$0 --reference [reference.fasta] --num-genomes [number of genomes] --num-variants [number of positions] --random-seed [random seed] [--repeat-positions positions.tsv]\n".
"Parameters:\n".
"\t--reference:  A reference genome in FASTA format\n".
"\t--num-variant-genomes: Number of genomes with variation to include in the table\n".
"\t--num-duplicate-reference-genomes: Number of genomes identical to the reference to include in the table\n".
"\t--num-substitutions: Number of all substitution positions to generate in table\n".
"\t--num-insertions: Number of single insertion and other substitution positions to generate in table\n".
"\t--num-deletions: Number of single deletion and all other substitution positions to generate in table\n".
"\t--random-seed: Random seed for generating mutations\n".
"\t--repeat-positions: A file of repeat positions to annotate when generating random variants.\n".
"Example:\n".
"$0 --reference reference.fasta --num-variant-genomes 10 --num-duplicate-reference-genomes 1 --num-substituions 9000 --num-insertions 500 --num-deletions 500 --random-seed 42 --repeat-positions repeats.tsv | sort -k 1,1 -k 2,2n > variants_table.tsv\n";

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
	print "\t${newbase}${base}";
}

sub get_unique_position
{
	my ($repeat_positions) = @_;

	my ($sequence_name,$pos,$sequence,$length_sequence);

	my $total_length = 0;
	my @sequence_start_length = ();
	for my $seq_name (@sequence_names)
	{
		my $sequence = $reference_table->{$seq_name};
		push(@sequence_start_length,$total_length);
		$total_length += $sequence->length;
	}

	# generate unique positions
	do 
	{
		# select random position, scaled by length of each sequence in the file
		my $position = int(rand($total_length));

		# find appropriate sequence in reference file
		my $index = 0;
		while(($index < @sequence_names) and ($position >= $sequence_start_length[$index]))
		{
			$index++;
		}
		my $chosen_index = $index - 1;

		$sequence_name = $sequence_names[$chosen_index];
		$sequence = $reference_table->{$sequence_name};

		$pos = $position - $sequence_start_length[$chosen_index];
	} while (exists $positions_used->{"${sequence_name}_${pos}"});
	$positions_used->{"${sequence_name}_${pos}"} = 1;

	my $seq_string = $sequence->seq;
	my $ref_base = substr($seq_string,$pos,1);
	my $in_repeat_region = exists $repeat_positions->{"${sequence_name}_${pos}"};

	return ($sequence_name,$pos,$ref_base,$in_repeat_region);
}

############
### MAIN ###
############
my ($ref_file,$num_variant_genomes,$num_duplicate_reference_genomes,$num_substitutions,$num_deletions,$num_insertions,$random_seed, $repeat_positions_file);

if (!GetOptions('reference=s' => \$ref_file,
                'num-substitutions=i' => \$num_substitutions,
                'num-insertions=i' => \$num_insertions,
                'num-deletions=i' => \$num_deletions,
                'num-variant-genomes=i' => \$num_variant_genomes,
		'num-duplicate-reference-genomes=i' => \$num_duplicate_reference_genomes,
                'random-seed=i' => \$random_seed,
		'repeat-positions=s' => \$repeat_positions_file)) {
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

my $repeat_positions = {};

if (defined $repeat_positions_file) {
	print STDERR "Will mark all positions in $repeat_positions_file as repeats\n";
	my $repeat_positions_parser = InvalidPositions->new;
	$repeat_positions = $repeat_positions_parser->read_invalid_positions($repeat_positions_file);
}


# read original reference sequences
$reference_table = read_reference_sequences($ref_file);
$reference_name = basename($ref_file, '.fasta');
@sequence_names = sort {$a cmp $b} (keys %$reference_table);
$number_sequences = scalar(@sequence_names);

print_header_line($num_duplicate_reference_genomes,$num_variant_genomes,$reference_name);

# substitution positions
for (my $pos_num = 0; $pos_num < $num_substitutions; $pos_num++)
{
	my ($sequence_name,$pos,$ref_base,$in_repeat_position) = get_unique_position($repeat_positions);

	# print variant line, +1 to position since positions start with 1, not 0
	if ($in_repeat_position)
	{
		print "$sequence_name\t".($pos+1)."\trepeat,substitution\t$ref_base";
	}
	else
	{
		print "$sequence_name\t".($pos+1)."\tvalid,substitution\t$ref_base";
	}

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
	my ($sequence_name,$pos,$ref_base,$in_repeat_position) = get_unique_position($repeat_positions);

	# print variant line, +1 to position since positions start with 1, not 0
	if ($in_repeat_position)
	{
		print "$sequence_name\t".($pos+1)."\trepeat,deletion\t$ref_base";
	}
	else
	{
		print "$sequence_name\t".($pos+1)."\tvalid,deletion\t$ref_base";
	}

	# for each duplicate reference genome to generate
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print_reference($ref_base);
	}

	# pick a random replicate to have a deletion
	my $replicate_with_deletion = int(rand($num_variant_genomes));

	# for each genome to generate
	for (my $i = 0; $i < $num_variant_genomes; $i++)
	{
		if ($i == $replicate_with_deletion)
		{
			print_deletion();
		}
		else
		{
			print_substitution($ref_base);
		}
	}
	print "\n";
}

# insertions
for (my $pos_num = 0; $pos_num < $num_insertions; $pos_num++)
{
	my ($sequence_name,$pos,$ref_base,$in_repeat_position) = get_unique_position($repeat_positions);

	# print variant line, +1 to position since positions start with 1, not 0
	if ($in_repeat_position)
	{
		print "$sequence_name\t".($pos+1)."\trepeat,insertion\t$ref_base";
	}
	else
	{
		print "$sequence_name\t".($pos+1)."\tvalid,insertion\t$ref_base";
	}

	# for each duplicate reference genome to generate
	for (my $i = 0; $i < $num_duplicate_reference_genomes; $i++)
	{
		print_reference($ref_base);
	}

	my $replicate_with_insertion = int(rand($num_variant_genomes));

	for (my $i = 0; $i < $num_variant_genomes; $i++)
	{
		if ($i == $replicate_with_insertion)
		{
			print_insertion($ref_base);
		}
		else
		{
			print_substitution($ref_base);
		}
	}
	print "\n";
}
