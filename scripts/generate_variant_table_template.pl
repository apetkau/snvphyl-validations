#!/usr/bin/env perl
# Purpose
# Given a reference fasta file, number of genomes to generate, and number of positions generates a 'template' variant table to be filled in manually.

use warnings;
use strict;

use Bio::SeqIO;
use File::Basename;

my $usage =
"$0 [reference.fasta] [number of genomes] [number of positions]\n".
"Parameters:\n".
"\treference.fasta:  A reference genome in FASTA format\n".
"\tnumber of genomes: Number of genomes to include in template table\n".
"\tnumber of positions: Number of positions to mutate in table\n".
"Example:\n".
"$0 reference.fasta 5 100 | sort -n -k 2,2 > variants_table.tsv\n";

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
die $usage if (@ARGV != 3);
my ($ref_file,$num_genomes,$num_positions) = @ARGV;
srand(time);

die "reference.fasta is not defined\n$usage" if (not defined $ref_file);
die "$ref_file does not exist\n$usage" if (not -e $ref_file);
die "number of genomes is not defined\n$usage" if (not defined $num_genomes);
die "number of genomes=$num_genomes is not valid\n$usage" if ($num_genomes !~ /^\d+$/);
die "number of positions is not defined\n$usage" if (not defined $num_positions);
die "number of positions=$num_positions is not valid\n$usage" if ($num_positions !~ /^\d+$/);

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
	# select random sequence
	my $seq_num = int(rand($number_sequences));

	my $sequence_name = $sequence_names[$seq_num];
	my $sequence = $reference_table->{$sequence_name};
	my $length_sequence = $sequence->length;

	# select random position
	my $pos = int(rand($length_sequence));

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
