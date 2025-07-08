#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::Frescobi::Sequtil;

my $prefix = "";
my $size = 4000000000;

GetOptions("prefix=s" => \$prefix,
	   "size=i" => \$size);

my $usage = <<EOF;
split_fasta_by_sie.pl -prefix=<string> [fasta file]
                      [ -size=<integer> ]  Default: $size
EOF
    ;

if ($prefix eq "" or scalar(@ARGV) > 1) {
    print STDERR $usage;
    exit(1);
}

if (scalar(@ARGV) == 1) {
    open (STDIN, "< $ARGV[0]") ||
	die "Unable to open $ARGV[0]: $!\n";
}

my $bp_count = 0;
my $file_count = 1;

&open_out;

my $line = <STDIN>;
while (defined $line) {
    my ($name, $annotation, $seq) = read_1_fasta_sequence(\*STDIN, $line, 0);
    last if not defined $name;
    if (length($seq) > $size) {
	print STDERR "Sequence $name is bigger than $size.\n";
    }
    $bp_count += length($seq);
    # In the case we get a sequence that's bigger than $size, we will output it.
    if ($bp_count > $size) {
	if ($bp_count != length($seq)) { #  This test avoids a zero length file.
	    &open_out;
	    $bp_count = length($seq);
	}
    }
    printf OUT ">$name%s\n", length($annotation) == 0 ? "" : " $annotation";
    print OUT format_seq($seq);
}

sub open_out {
    my $filename = "${prefix}_${file_count}.fa";
    open (OUT, ">$filename") ||
	die "Unable to open $filename: $!\n";
    $file_count += 1;
}
