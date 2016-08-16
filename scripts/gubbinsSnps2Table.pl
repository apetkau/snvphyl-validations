#!/usr/bin/env perl
# Purpose: Converts Gubbins .summary_of_snp_distribution.vcf to SNVPhyl SNVTable with same order of columns

use strict;
use warnings;

use Bio::SeqIO;
use Set::Scalar;
use Getopt::Long;

sub get_order_snvphyl_table {
	my ($snvphyl_table) = @_;

	open(my $fh,"<$snvphyl_table") or die "Could not open $snvphyl_table\n";
	my $header = readline($fh);
	chomp($header);
	die "Error with header line: $header\n" if ($header !~ /^#/);
	my ($chrom,$position,$status,$reference, @genome_order) = split(/\t/,$header);
	close($fh);

	return @genome_order;
}


my $usage = "$0 --snvphyl-table [variants.tsv] --gubbins-table [name.summary_of_snp_distribution.vcf]\n";

my ($snvphyl_table,$gubbins_table);

if (!GetOptions('snvphyl-table=s' => \$snvphyl_table,
		'gubbins-table=s' => \$gubbins_table)) {
	die "Invalid option\n".$usage;
}

die "--snvphyl-table not defined\n$usage" if (not defined $snvphyl_table);
die "--gubbins-table not defined\n$usage" if (not defined $gubbins_table);

my @snvphyl_order = get_order_snvphyl_table($snvphyl_table);
open (my $gfh, "<$gubbins_table") or die "Could not open $gubbins_table\n";
my $header_line;
for ($header_line = readline($gfh); $header_line !~ /^#CHROM/; $header_line = readline($gfh)) {}
chomp ($header_line);
my (@gubbins_fields) = split(/\t/,$header_line);

# generate map of 'genome => @gubbins_fields position in array'
my %genome_gubbins_field;
for (my $i = 0; $i < @gubbins_fields; $i++) {
	$genome_gubbins_field{$gubbins_fields[$i]} = $i;
}

# validate we have entries for every snvphyl genome
for (my $i = 0; $i < @snvphyl_order; $i++) {
	die "Error: no entry in gubbins for snvphyl genome ".$snvphyl_order[$i] if (not exists $genome_gubbins_field{$snvphyl_order[$i]});
}

print "#Chromosome\tPosition\tStatus\tReference\t".join("\t",@snvphyl_order)."\n";
while (my $line = readline($gfh)) {
	chomp ($line);
	my (@gubbins_fields) = split(/\t/,$line);
	my ($chr,$pos,$id,$ref,$alt) = @gubbins_fields;

	print "{chrom}\t$pos\tvalid\t$ref";
	for my $genome (@snvphyl_order) {
		my $gubbins_base = $gubbins_fields[$genome_gubbins_field{$genome}];
		print "\t$gubbins_base";
	}
	print "\n";
}

close($gfh);
