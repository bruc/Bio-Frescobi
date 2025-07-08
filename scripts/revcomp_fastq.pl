#!perl

use strict;
use warnings;

use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;
use Getopt::Long;

my $usage = <<EOF;
fastq_revcomp.pl input-fastq-file output-fastq-file

Reverse complements a FASTQ file.
EOF
    ;

my $count = 0;

if (scalar(@ARGV) != 2) {
    print STDERR $usage;
    exit(1);
}

open (IN, "<$ARGV[0]") ||
    die "Unable to open $ARGV[0]: $!\n";
open (OUT, ">$ARGV[1]") ||
    die "Unable to open $ARGV[1]: $!\n";


while (my ($name, $annotation, $seq, $qual) = read_1_fastq_sequence(\*IN, 'B')) {
    my $new_qual = reverse($qual);
    my $new_seq = reverse_dna_seq($seq);
    print OUT format_fastq($name, $annotation, $new_seq, $new_qual);
    $count++;
}

print STDERR "A total of $count sequences were reversed.\n";

=head1 NAME

fastq_revcomp.pl - Tool for reverse complementing a FASTQ file

=head1 SYNOPSIS

 fastq_revcomp.pl <input-fastq-file> <output-fastq-file>

=head1 DESCRIPTION

The C<fastq_revcomp.pl> script provides a simple tool for reverse
complementing the sequences in a FASTQ file.

=head1 ARGUMENTS AND OPTIONS

The program just takes two operands, the input and output FASTQ files.

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
