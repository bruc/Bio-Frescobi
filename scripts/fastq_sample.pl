#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::Frescobi::Sequtil;

my $sample = -1;
my $seed = 1;
my $output = "";

GetOptions("sample=i" => \$sample,
	   "seed=i" => \$seed,
	   "output=s" => \$output);

my $usage = <<EOF;
fastq_sample  -sample=<integer> -output=<string> fastq-file-1 [fastq-file-2]
            [ -seed=<integer> ]

The -sample option must be positive.
The -output option is a file name prefix.    
EOF
    ;

if (scalar(@ARGV) < 1 or scalar(@ARGV) > 2 or $sample <= 0) {
    print STDERR $usage;
    exit(1);
}

open (FQ1, "<$ARGV[0]") ||
    die "Unable to open $ARGV[0]: $!\n";

if ($ARGV[1]) {
    open (FQ2, "<$ARGV[1]") ||
	die "Unable to open $ARGV[1]: $!\n";
}

srand($seed);

my $n1 = 0;
my $n2 = 0;
while (my ($name1, $annotation1, $seq1, $qual1) = read_1_fastq_sequence(\*FQ1, 'B')) {
    $n1++;
    if ($ARGV[1]) {
	my ($name2, $annotation2, $seq2, $qual2) = read_1_fastq_sequence(\*FQ2, 'B');
	if (not defined $name2) {
	    print STDERR "The two fastq files are not the same length\n";
	    exit(1);
	}
	$n2++;
	my $id1 = (split(/\//, $name1))[0];
	my $id2 = (split(/\//, $name2))[0];
	if ($id1 ne $id2) {
	    print STDERR "For sequence $n1, id's do not match. id1 = $id1  id2 = $id2\n";
	    exit(1);
	}
    }
}

if ($ARGV[1]) {
    if (read_1_fastq_sequence(\*FQ2, 'B')) {
	$n2++;
    }
}

if ($ARGV[1]) {
    if ($n1 != $n2) {
	print STDERR "The two fastq files are not the same length\n";
	exit(1);
    }
}

if ($sample > $n1) {
    print STDERR "Sample size is larger than the FASTQ sequence count.\n";
    exit(1);
}

print STDERR "$n1 sequences (or pairs) read.\n";

my %selected;

my $sel_count = 0;
while ($sel_count < $sample) {
    my $i = int(rand($n1));
    if (not exists($selected{$i})) {
	$selected{$i} = 1;
	$sel_count++;
    }
}

my $output_file1;
my $output_file2;

if ($ARGV[1]) {
    $output_file1 = $output . ".1.fq";
    $output_file2 = $output . ".2.fq";
    open (OUT2, ">$output_file2") ||
	die "Unable to open $output_file2: $!\n";
}
else {
    $output_file1 = $output . ".fq";
}
open (OUT1, ">$output_file1") ||
    die "Unable to open $output_file1: $!\n";

open (FQ1, "<$ARGV[0]") ||
    die "Unable to open $ARGV[0]: $!\n";

if ($ARGV[1]) {
    open (FQ2, "<$ARGV[1]") ||
	die "Unable to open $ARGV[1]: $!\n";
}

$n1 = 0;
$n2 = 0;
while (my ($name1, $annotation1, $seq1, $qual1) = read_1_fastq_sequence(\*FQ1, 'B')) {
    if ($selected{$n1}) {
	print OUT1 format_fastq($name1, $annotation1, $seq1, $qual1);
    }
    if ($ARGV[1]) {
	my ($name2, $annotation2, $seq2, $qual2) = read_1_fastq_sequence(\*FQ2, 'B');
	if ($selected{$n1}) {
	    print OUT2 format_fastq($name2, $annotation2, $seq2, $qual2);
	}
    }
    $n1++;
}


