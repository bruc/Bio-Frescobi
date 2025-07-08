#!/usr/bin/perl

# Simple script to count amino acid occurences in a FASTA file.

use strict;
use warnings;

use Bio::Frescobi::Sequtil;

my $line = <STDIN>;
my %counts;

# Put in zeroes for the known amino acids.

foreach my $aa (split(//, "GPAVLIMSTKRDENQFYWCH")) {
    $counts{$aa} = 0;
}

my $records = 0;
while (defined $line) {
    my ($name, $annotation, $seq) =
	read_1_fasta_sequence(\*STDIN, $line, 0);
    last if not defined $name;
    foreach my $aa (split(//, $seq)) {
	$counts{uc($aa)} += 1;
    }
    if (++$records % 1000 == 0) {
	print STDERR "$records sequences processed.\n";
    }
}

print STDERR "$records sequences processed.\n";
foreach my $aa (sort {$a cmp $b} keys %counts) {
    printf "%s\t%d\n", $aa, $counts{$aa};
}

