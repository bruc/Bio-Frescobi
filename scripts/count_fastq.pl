#!perl

use strict;
use warnings;

use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;
use Getopt::Long;

my $tmpdir = $ENV{TMPDIR};
my $mincount = 2;
my $usememory = 1;
my $sortmem = "4G";
my $trim = -1;

GetOptions("tmpdir=s" => \$tmpdir,
	   "mincount=i" => \$mincount,
	   "usememory!" => \$usememory,
	   "trim=i" => \$trim,
	   "sortmem=s" => \$sortmem);

my $count = 0;
my %counts;
if ($usememory) {
    while (my ($name, $annotation, $seq, $qual) =
	   read_1_fastq_sequence(\*STDIN, 'b')) {
	if ($trim > 0) {
	    $seq = substr($seq, 0, $trim);
	}
	$counts{$seq}++;
	$count++;
    }
    print STDERR "$count sequences read.\n";
}
else {
    if (not $tmpdir) {
	$tmpdir = "/tmp";
    }
    if (not -d $tmpdir) {
	system_with_check("mkdir -p $tmpdir", 1);
    }
    my $sorted_file = "${tmpdir}/count_fastq.$$.inp";
    open (PIPE,
	  "| sort --buffer-size=${sortmem} " .
	  " --stable --temporary-directory=${tmpdir} --ignore-case " .
	  " >${sorted_file}") ||
	die "Unable to open sort pipe: $!\n";
    while (my ($name, $annotation, $seq, $qual) =
	   read_1_fastq_sequence(\*STDIN, 'b')) {
	if ($trim > 0) {
	    $seq = substr($seq, 0, $trim);
	}
	print PIPE "$seq\n";
	print STDERR "$count reads processed.\n" if ++$count % 10000 == 0;
    }
    print STDERR "$count sequences read.\n";
    close PIPE;
    my $prev_seq = "";
    my $seq_count = 0;
    open (SORT, "<${sorted_file}") ||
	die "Unable to open ${sorted_file}: $!\n";
    while (<SORT>) {
	chomp;
	$_ = uc($_);
	if ($prev_seq eq $_) {
	    $seq_count++;
	}
	else {
	    if ($prev_seq) {
		if ($seq_count >= $mincount) {
		    $counts{$prev_seq} = $seq_count;
		}
	    }
	    $prev_seq = $_;
	    $seq_count = 1;
	}
    }
    close SORT;
    if ($prev_seq) {
	if ($seq_count >= $mincount) {
	    $counts{$prev_seq} = $seq_count;
	}
    }
#    system_with_check("rm ${sorted_file}");
}
foreach my $seq (sort {$counts{$b} <=> $counts{$a} ||
		       $a cmp $b } keys %counts) {
    if ($counts{$seq} >= $mincount) {
	printf "%8d: %s\n", $counts{$seq}, $seq;
    }
}


=head1 NAME

count_fastq.pl - Tool for counting recurrences in a FASTQ file

=head1 SYNOPSIS

 count_fastq.pl [-tmpdir=<directory> ] < FASTQ-file
                [-mincount=<integer> ]
                [-[no]usememory      ]
                [-trim=<integer>     ]

=head1 DESCRIPTION

The C<count_fastq.pl> script provides a simple tool for counting
recurrences in a FASTQ file. It is very handy for identifying repeats
in a sequence library.

For short Fastq files, it can use Perl memory to store sequences, and
for longer files, it can use the system C<sort> program to handle any
size file that can fit in the available storage.

=head1 ARGUMENTS AND OPTIONS

The standard input to the program is a Fastq file, and the standard
output is a list of most frequent sequences sorted for most abundant
to least. Messages about progress and the number of reads are written
to standard error.

The options for C<count_fastq.pl> are as follows:

=over 4

=item C<< -tmpdir=<directory> >>

This option specifies the temporary directory for sorting. If not
specified, then the value of the environment variable, C<TMPDIR>, is
used. If there is no such environment variable, then C</tmp> is used.

=item C<< -mincount=<integer> >>

The C<-mincount> option specifies the minimum number of duplicates
needed for display. The default is 2.

=item C<< -[no]usememory >>

This option controls whether the counting is done in Perl's memory
(C<-usememory>) or by using the Unix C<sort> program
(C<-nousememory>).

=item C<< -trim=<integer> >>

This option specifies that each sequence should be trimmed to the
specified length.  If the option is specified as a negative number,
then the full length of each read is used.

=item C<< -sortmem=<string> >>

=back

=head1 SEE ALSO

  sort(1)

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
