#!/usr/bin/perl

# Retrieve the first sequence from the fasta file and print it.

use strict;
use warnings;
use Bio::Frescobi::Sequtil;

open (SEQ, "<$ARGV[0]") ||
    die "Unable to open $ARGV[0]: $!\n";

$_ = <SEQ>;
my ($name, $annotation, $seq) = read_1_fasta_sequence(\*SEQ, $_, 0);
if (not defined $name) {
    die "No sequences in $ARGV[0]\n";
}
if ($annotation ne "") {
    $annotation = " " . $annotation;
}
print ">${name}${annotation}\n";
print format_seq($seq);

