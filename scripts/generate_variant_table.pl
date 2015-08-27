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
use InvalidPositions;

my $usage =
"$0 --reference [reference.fasta] --num-genomes [number of genomes] --num-variants [number of positions] --random-seed [random seed] [--exclude-positions positions.tsv]\n".
"Parameters:\n".
"\t--reference:  A reference genome in FASTA format\n".
"\t--num-genomes: Number of genomes to include in the table\n".
"\t--num-variants: Number of variants to mutate in table\n".
"\t--random-seed: Random seed for generating mutations\n".
"\t--exclude-positions: A file of positions to exclude when generating random variants (if no parameter is passed select from all positions on reference).\n".
"Example:\n".
"$0 --reference reference.fasta --num-genomes 5 --num-variants 100 --random-seed 42 --exclude-positions repeats.tsv | sort -k 1,1 -k 2,2n > variants_table.tsv\n";

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

############
### MAIN ###
############
my ($ref_file,$num_genomes,$num_positions,$random_seed, $excluded_positions_file);

if (!GetOptions('reference=s' => \$ref_file,
                'num-variants=i' => \$num_positions,
                'num-genomes=i' => \$num_genomes,
                'random-seed=i' => \$random_seed,
		'exclude-positions=s' => \$excluded_positions_file)) {
        die "Invalid option\n".$usage;
}

die "reference.fasta is not defined\n$usage" if (not defined $ref_file);
die "$ref_file does not exist\n$usage" if (not -e $ref_file);
die "number of genomes is not defined\n$usage" if (not defined $num_genomes);
die "number of genomes=$num_genomes is not valid\n$usage" if ($num_genomes !~ /^\d+$/);
die "number of positions is not defined\n$usage" if (not defined $num_positions);
die "number of positions=$num_positions is not valid\n$usage" if ($num_positions !~ /^\d+$/);

if (not defined $random_seed) {
	$random_seed = 42;
	warn "--random-seed not defined, defaulting to $random_seed\n";
}


srand($random_seed);

my $positions_used = {};
if (defined $excluded_positions_file) {
	print STDERR "Will exclude all positions in $excluded_positions_file\n";
	my $invalid_positions_parser = InvalidPositions->new;
	$positions_used = $invalid_positions_parser->read_invalid_positions($excluded_positions_file);
}

my $swap_table = {
'A' => 'T',
'T' => 'G',
'G' => 'C',
'C' => 'A'};

# read original reference sequences
my $reference_table = read_reference_sequences($ref_file);

my $reference_name = basename($ref_file, '.fasta');
my @sequence_names = (keys %$reference_table);
my $number_sequences = scalar(@sequence_names);

# print header line
print "#Chromosome\tPosition\tStatus\tReference";
# for each genome
for (my $i = 0; $i < $num_genomes; $i++)
{
	print "\t$reference_name-$i";
}
print "\n";

for (my $pos_num = 0; $pos_num < $num_positions; $pos_num++)
{
	my ($seq_num,$pos,$sequence_name,$sequence,$length_sequence);
	# generate unique positions
	do 
	{
		# select random sequence
		$seq_num = int(rand($number_sequences));

		$sequence_name = $sequence_names[$seq_num];
		$sequence = $reference_table->{$sequence_name};
		$length_sequence = $sequence->length;

		# select random position
		$pos = int(rand($length_sequence));
	} while (exists $positions_used->{"${sequence_name}_${pos}"});
	$positions_used->{"${sequence_name}_${pos}"} = 1;

	my $seq_string = $sequence->seq;
	my $ref_base = substr($seq_string,$pos,1);

	# print variant line, +1 to position since positions start with 1, not 0
	print "$sequence_name\t".($pos+1)."\tvalid\t$ref_base";

	# for each genome to generate
	for (my $i = 0; $i < $num_genomes; $i++)
	{
		# cluster genomes
		if ($i % 2 == 0)
		{
			if ($i % 4 == 0) # every 4th
			{
				if ($pos_num % 2 == 0) # 50% of positions
				{
					print "\t".$swap_table->{uc($ref_base)};
				}
				else
				{
					# swap again for other 50%
					print "\t".$swap_table->{$swap_table->{uc($ref_base)}};
				}
			}
			else
			{
				print "\t".$swap_table->{uc($ref_base)};
			}
		}
		else
		{
			print "\t$ref_base";
		}
	}
	print "\n";
}
