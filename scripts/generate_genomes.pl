#!/usr/bin/env perl
# Purpose
# Given a reference FASTA file and a table of variants, generates the corresponding genomes for these variants.

use warnings;
use strict;

use FindBin;
use lib $FindBin::Bin.'/lib';

use Getopt::Long;
use File::Copy 'move';
use PositionsTable;
use Bio::SeqIO;
use Bio::Seq::SeqFactory;

my $usage =
"$0 --reference [reference.fasta] --variants [variants_table.tsv] --min-cov [minimum coverage] --max-cov [maximum coverage] --art-parameters [art_illumina_parameters] --out-dir [out_dir]\n".
"Parameters:\n".
"\t--reference:	A reference genome in FASTA format\n".
"\t--variants: A tab-deliminited table listing the genomes and variants to insert (generated from generate_variant_table_template.pl)\n".
"\t--min-cov: The minium coverage of any simulated set of reads\n".
"\t--max-cov: The maximum coverage of any simulated set of reads\n".
"\t--random-seed: Random seed for defining coverage\n".
"\t--art-parameters: A quoted string containing the art_illumina parameters, minus --fcov (coverage).\n".
"\t--out-dir: The directory to store all output files\n".
"Example:\n".
"$0 --reference reference.fasta --variants variants.tsv --min-cov 30 --max-cov 60 --random-seed 42 --art-parameters '--noALN --len 250 --rndSeed 42 --paired --seqSys MS --sdev 50 --mflen 400' --out-dir out/\n\n";

my $art_illumina_bin = `which art_illumina 2>/dev/null`;
chomp $art_illumina_bin;

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
my $min_cov_default = 30;
my $max_cov_default = 60;
my $paired = 0;
my ($ref_file,$variants_file,$art_illumina_parameters,$out_dir,$min_coverage,$max_coverage,$random_seed);
if (!GetOptions('reference=s' => \$ref_file,
		'variants=s' => \$variants_file,
		'min-cov=i' => \$min_coverage,
		'max-cov=i' => \$max_coverage,
		'random-seed=i' => \$random_seed,
		'art-parameters=s' => \$art_illumina_parameters,
		'out-dir=s' => \$out_dir)) {
	die "Invalid option\n".$usage;
}

die "error: dependency 'art_illumina' not found on PATH, please install this software and try again" if ((not defined $art_illumina_bin) or not (-e $art_illumina_bin));

die "reference.fasta is not defined\n$usage" if (not defined $ref_file);
die "$ref_file does not exist\n$usage" if (not -e $ref_file);
die "variants_table.tsv is not defined\n$usage" if (not defined $variants_file);
die "$variants_file does not exist\n$usage" if (not -e $variants_file);
die "out_dir is not defind\n$usage" if (not defined $out_dir);
die "$out_dir already exists\n$usage" if (-e $out_dir);
die "art_illumina_parameters not defined\n$usage" if (not defined $art_illumina_parameters);
die "art_illumina_parameters contains --fcov, this is set in code\n$usage" if ($art_illumina_parameters =~ /--fcov[\s\$]/ or $art_illumina_parameters =~ /-f[\s\$]/);

if (not defined $min_coverage) {
	$min_coverage = $min_cov_default;
	warn "--min-cov not set, defaulting to $min_cov_default\n";
}

if (not defined $max_coverage) {
	$max_coverage = $max_cov_default;
	warn "--max-cov not set, defaulting to $max_cov_default\n";
}

if (not defined $random_seed) {
	$random_seed = 42;
	warn "--random-seed not set, defaulting to 42\n";
}

$paired = ($art_illumina_parameters =~ /-p[\s\$]/ or $art_illumina_parameters =~ /--paired[\s\$]/);

my $genomes_out_dir = "$out_dir/genomes";

# read in all positions to mutate
my $variants_parser = PositionsTable->new;
my ($variants_core_table,$variants_core_table_count) = $variants_parser->read_table($variants_file);
my @genome_names = (keys %$variants_core_table);
die "error: no genomes to generate defined in $variants_file" if (@genome_names == 0);

# read original reference sequences
my $reference_table = read_reference_sequences($ref_file);

print "Running $0 with parameters\n";
print "reference_file=$ref_file\n";
print "variants_table=$variants_file\n";
print "art_illumina_parameters=$art_illumina_parameters\n";
print "out_dir=$out_dir\n";
print "paired-end=".($paired ? 'true' : 'false');
print "\n";

mkdir $out_dir;
mkdir $genomes_out_dir;

my $factory = Bio::Seq::SeqFactory->new;

srand($random_seed);
foreach my $genome_name (sort {$a cmp $b} @genome_names)
{
	next if ($genome_name eq 'Reference');
	my $coverage = int(rand($max_coverage-$min_coverage))+$min_coverage;

	my $genome_file = "$genomes_out_dir/$genome_name.fasta";
	my $fastq_file_prefix = ($paired ? "$out_dir/${genome_name}_R" : "$out_dir/${genome_name}_");

	# mutate each new genome, write to file
	my $variant_entry = $variants_core_table->{$genome_name};
	for my $chrom (keys %$variant_entry)
	{
		my $positions = $variant_entry->{$chrom};
		die "error: no chromsome named $chrom in file $ref_file"
			if (not defined $reference_table->{$chrom});

		my $reference_seq = $reference_table->{$chrom};
		my $seq_string = $reference_seq->seq;

		for my $pos (keys %$positions)
		{
			# insert mutation at position
			my $string_pos = $pos - 1; # position in string starts at 0 not 1
			my $alt = $positions->{$pos}->{'alternative'};
			my $ref = $positions->{$pos}->{'reference'};
			die "error: no alt for $genome_name:$chrom:$pos" if (not defined $alt);
			die "error: invalid alt=$alt for $genome_name:$chrom:$pos" if ($alt !~ /^[ACTG]$/i);
			die "error for $genome_name:$chrom:$pos, position out of bounds in file $ref_file" if ($pos > length($seq_string));

			my $real_ref_base = substr($seq_string,$string_pos,1);
			die "error for $genome_name:$chrom:$pos base($real_ref_base) from file $ref_file != base($ref) from file $variants_file" if (lc($real_ref_base) ne lc($ref));

			# for deletion, make sure to delete alt base
			if ($alt eq '-')
			{
				substr($seq_string,$string_pos,1) = ''; # deletion
			}
			else
			{
				substr($seq_string,$string_pos,1) = $alt; # perform mutation (grows sequence string for insertion)
			}
		}

		my $generated_sequence = $factory->create(-seq => $seq_string, -id => $chrom);
		my $seq_writer = Bio::SeqIO->new(-file=>">>$genome_file",-format=>"fasta");
		$seq_writer->write_seq($generated_sequence);
	}
	print "wrote $genome_file\n";

	# generate reads
	my $art_command = "$art_illumina_bin -i $genome_file $art_illumina_parameters --fcov $coverage -o $fastq_file_prefix 1> $out_dir/$genome_name.cov-$coverage.log 2>&1";
	print "running '$art_command'\n";
	system($art_command) == 0 or die "Could not execute '$art_command'\n";
	if ($paired) {
		move("${fastq_file_prefix}1.fq", "${fastq_file_prefix}1.fastq");
		move("${fastq_file_prefix}2.fq", "${fastq_file_prefix}2.fastq");
		print "wrote ${fastq_file_prefix}1.fastq\n";
		print "wrote ${fastq_file_prefix}2.fastq\n";
	} else {
		move("$fastq_file_prefix.fq", "$fastq_file_prefix.fastq");
		print "wrote ${fastq_file_prefix}1.fastq\n";
	}
	print "\n";
}

