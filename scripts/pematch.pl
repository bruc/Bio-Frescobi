#!perl

use strict;
use warnings;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::PEMatch;
use Getopt::Long;

my $file1 = "";
my $file2 = "";
my $minoffset = 18;
my $maxoffset = 22;
my $minoverlap = 10;
my $maxd = 5;
my $qualbase = 'B';
my $match_weight = 5;
my $badfile = "";
my $debug = 0;

my $usage = <<EOF;
pematch.pl  -bad=<file> fastq-file1 fastq-file2
           [-minoffset=<integer> ]  Default: $minoffset
           [-maxoffset=<integer> ]  Default: $maxoffset
           [-maxd=<integer>      ]  Default: $maxd
           [-debug=<integer>     ]  Default: $debug
           [-qualbase=<char>     ]  Default: $qualbase
           [-matchweight=<float> ]  Default: $match_weight
           [-minoverlap=<integer>]  Default: $minoverlap
EOF
;
GetOptions("minoffset=i" => \$minoffset,
	   "maxoffset=i" => \$maxoffset,
	   "maxd=i" => \$maxd,
	   "debug=i" => \$debug,
	   "qualbase=s" => \$qualbase,
	   "matchweight=f" => \$match_weight,
	   "minoverlap=i" => \$minoverlap,
	   "bad=s" => \$badfile);

if (scalar(@ARGV) != 2
    or not $badfile) {
    print STDERR $usage;
    exit(1);
}

open (FQ1, "<$ARGV[0]") ||
    die "Unable to open $ARGV[0]: $!\n";
open (FQ2, "<$ARGV[1]") ||
    die "Unable to open $ARGV[1]: $!\n";
open (BAD, ">$badfile") ||
    die "Unable to open $badfile: $!\n";

my $matcher = new Bio::Frescobi::PEMatch;
$matcher->minoffset($minoffset);
$matcher->maxoffset($maxoffset);
$matcher->maxd($maxd);
$matcher->minoverlap($minoverlap);
$matcher->qualbase($qualbase);
$matcher->match_weight($match_weight);
$matcher->debug($debug);

while (my ($name1, $annot1, $seq1, $qual1) = &read_1_fastq_sequence(\*FQ1)) {
    $annot1 = "" if not defined $annot1;
    my ($name2, $annot2, $seq2, $qual2) = &read_1_fastq_sequence(\*FQ2);
    print STDERR "Now processing $name1 and $name2\n" if $debug;
    if (not defined $name2) {
	die "Mismatched number of sequences: too few in $ARGV[1]\n";
    }
    $annot2 = "" if not defined $annot2;
    my ($id1, $dir1) = split(/\//, $name1);
    my ($id2, $dir2) = split(/\//, $name2);
    if ($id1 ne $id2 or
	$dir1 ne '1' or
	$dir2 ne '2') {
	die "Bad sequence id pairs: name1 = $name1  name2 = $name2\n";
    }
    my ($join_seq, $join_qual) = $matcher->match($seq1, $qual1, $seq2, $qual2);
    if (not defined $join_seq) {
	print BAD format_fastq($name1, $annot1, $seq1, $qual1);
	print BAD format_fastq($name2, $annot2, $seq2, $qual2);
    }
    else {
	print format_fastq($id1,
			   $annot1 . $matcher->equal,
			   $join_seq,
			   $join_qual);
    }
}

=head1 NAME

pematch.pl - Command line tool for Paired End Matching.

=head1 SYNOPSIS

  pematch.pl  -bad=<file> fastq-file1 fastq-file2 >merged-fastq-file
             [-minoffset=<integer> ]  Default: $minoffset
             [-maxoffset=<integer> ]  Default: $maxoffset
             [-maxd=<integer>      ]  Default: $maxd
             [-debug=<integer>     ]  Default: $debug
             [-qualbase=<char>     ]  Default: $qualbase
             [-matchweight=<float> ]  Default: $match_weight
             [-minoverlap=<integer>]  Default: $minoverlap

=head1 DESCRIPTION

The C<pematch.pl> script is a command line version of the
L<Bio::Frescobi::PEMatch> module, which implements an algorithm for
joining two reads together from a paired end read, typically from
Illumina sequencing, where the two ends of the read are taken from a
fragment that is shorter than twice the read length. Please refer to that
module for details of the algorithm.

The script requires two FASTQ files which must be in register, so that
each element of a paired end read is in the same order in both
files. The identifiers must match, with the exception of the C<\1> and
C<\2> designations for each half of the paired end.  The output is
written to standard output. Any pairs of reads that cannot be joined
are written to the file specified by the C<-bad> option so you can
investigate further and adjust the options.

It should be noted that matching up overlapping segments of two reads
is not a well defined problem because sequences can be repetitive or
overlaps can result by chance, so you must take care to have a long
enough overlapping region to avoid artifacts.

=head1 ARGUMENTS AND OPTIONS

The options for C<pematch.pl> are taken from directly from the
L<Bio::Frescobi::PEMatch> module. Summarizing:

=over 4

=item C<< fastq-file1 >>

FASTQ formatted file with the first element of each paired end read.

=item C<< fastq-file2 >>

FASTQ formatted file with the second element of each paired end
read. The order must match C<fastq-file1>.

=item C<< -bad=<file> >>

Specify the file name for the file of bad sequence pairs that cannot be matched.

=item C<< -minoffset=<integer> >>

Set the minimum offset for the overlap of the second sequence on the first.

=item C<< -maxoffset=<integer> >>

Set the maximum offset for the overlap. If set to a negative value,
then the option, C<-minoverlap>, sets the maximum offset based on the
sequence length.

=item C<< -minoverlap=<integer> >>

If C<-maxoffset> is negative, the C<-minoverlap> parameter specifies
the minimum allowed overlap to join two sequences. The value for this
parameter should be large enough to ensure that the overlap is
sufficiently unlikely to arise by chance.

=item C<< -maxd=<integer> >>

Sets the maximum allowed mismatch distance in the overlap region. The
mismatch distance is computed by adding the minimum number of
insertions, deletions, and substitutions needed to change the overlap
region of one sequence into the other.

=item C<< -matchweight=<float> >>

Sets the weight of the match distance relative to the length of the
overlap. This script attempts to maximize the following score:

 score = length(overlap) - matchweight * mismatch_distance

=item C<< -debug=<integer> >>

If this option is non-zero, debugging information is written to standard error.

=item C<< -qualbase=<char> >>

Specifies the minimum possible quality score as an ASCII
character. Typical values are 'B' (Ascii code 66) and '#' (Ascii code
35) depending on the version of the Illumina sequencer. The matching
code will substitute bases in the overlap region when the quality
scores have this value.

=back

=head1 SEE ALSO

  Bio::Frescobi::PEMatch(3)

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
