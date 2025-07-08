# This package contains some simple functions for sequence manipulation.

# Copyright 2010 Congenomics, LLC
# This package is free software; you can redistribute it and/or modify
# it under the same terms as Perl itself, either Perl version 5.8.8 or,
# at your option, any later version of Perl 5 you may have available.

package Bio::Frescobi::Sequtil;
use strict;
use warnings;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Carp;
use SHA;
use Bio::Frescobi::Genutil;
@ISA = ('Exporter');
@EXPORT = qw(&seq_checksum
	     &get_checksum_size
	     &set_checksum_size
	     &format_seq
	     &get_format_seq_size
	     &set_format_seq_size
	     &get_checkseq
	     &set_checkseq
	     &format_qual
	     &get_format_qual_size
	     &set_format_qual_size
	     &qual_to_array
	     &array_to_qual
	     &read_1_fasta_sequence
	     &extract_1_sequence
	     &read_1_fastq_sequence
             &read_fasta_file
             &format_fastq
	     &reverse_dna_seq
	     &substr_dna_seq
	     &is_totally_masked
	     &prot_3_to_1
	     &prot_1_to_3
	     &na_2_prot
	     &colorize_seq
             &degenerate_index
             &fastq_qual_to_array
             &fastq_array_to_qual
             &find_ns);

# The following functions are currently broken and will be fixed in a future release.

# 	     &packseq
#	     &unpackseq

$VERSION = '0.06';

# checksum_size is the default length of a checksum.

my $checksum_size = 8;

# format_seq_size is the default format_seq size.

my $format_seq_size = 50;

# format_qual_size is the default format_qual size.

my $format_qual_size = 50;

# checkseq is the mode variable controlling sequence checking.

my $checkseq = 1;

my %packseq = ("a" => "00",
	       "c" => "01",
	       "g" => "10",
	       "t" => "11");

my %unpackseq = reverse %packseq;

my %aa3to1 = ("ALA", "A",
	      "ARG", "R",
	      "ASN", "N",
	      "ASX", "B",
	      "ASP", "D",
	      "CYS", "C",
	      "GLN", "Q",
	      "GLX", "Z",
	      "PCA", "Q",
	      "GLU", "E",
	      "GLY", "G",
	      "HIS", "H",
	      "ILE", "I",
	      "LEU", "L",
	      "LYS", "K",
	      "MET", "M",
	      "PHE", "F",
	      "PRO", "P",
	      "HYP", "P",
	      "SER", "S",
	      "THR", "T",
	      "TRP", "W",
	      "TYR", "Y",
	      "VAL", "V",
	      "XAA", "U",
	      "UNK", "U");

my %aa1to3 = ("A", "ALA",
	      "R", "ARG",
	      "N", "ASN",
	      "D", "ASP",
	      "C", "CYS",
	      "Q", "GLN",
	      "E", "GLU",
	      "G", "GLY",
	      "H", "HIS",
	      "I", "ILE",
	      "L", "LEU",
	      "K", "LYS",
	      "M", "MET",
	      "F", "PHE",
	      "P", "PRO",
	      "S", "SER",
	      "T", "THR",
	      "W", "TRP",
	      "Y", "TYR",
	      "V", "VAL",
	      "*", "STOP");

our %genetic_code = (
    "UUU" => [ 'PHE', 'F' ],
    "UUC" => [ 'PHE', 'F' ],
    "UCU" => [ 'SER', 'S' ],
    "UCC" => [ 'SER', 'S' ],
    "UAU" => [ 'TYR', 'Y' ],
    "UAC" => [ 'TYR', 'Y' ],
    "UGU" => [ 'CYS', 'C' ],
    "UGC" => [ 'CYS', 'C' ],
    "UUA" => [ 'LEU', 'L' ],
    "UCA" => [ 'SER', 'S' ],
    "UAA" => [ '***', '*' ],
    "UGA" => [ '***', '*' ],
    "UUG" => [ 'LEU', 'L' ],
    "UCG" => [ 'SER', 'S' ],
    "UAG" => [ '***', '*' ],
    "UGG" => [ 'TRP', 'W' ],
    "CUU" => [ 'LEU', 'L' ],
    "CUC" => [ 'LEU', 'L' ],
    "CCU" => [ 'PRO', 'P' ],
    "CCC" => [ 'PRO', 'P' ],
    "CAU" => [ 'HIS', 'H' ],
    "CAC" => [ 'HIS', 'H' ],
    "CGU" => [ 'ARG', 'R' ],
    "CGC" => [ 'ARG', 'R' ],
    "CUA" => [ 'LEU', 'L' ],
    "CUG" => [ 'LEU', 'L' ],
    "CCA" => [ 'PRO', 'P' ],
    "CCG" => [ 'PRO', 'P' ],
    "CAA" => [ 'GLN', 'Q' ],
    "CAG" => [ 'GLN', 'Q' ],
    "CGA" => [ 'ARG', 'R' ],
    "CGG" => [ 'ARG', 'R' ],
    "AUU" => [ 'ILE', 'I' ],
    "AUC" => [ 'ILE', 'I' ],
    "ACU" => [ 'THR', 'T' ],
    "ACC" => [ 'THR', 'T' ],
    "AAU" => [ 'ASN', 'N' ],
    "AAC" => [ 'ASN', 'N' ],
    "AGU" => [ 'SER', 'S' ],
    "AGC" => [ 'SER', 'S' ],
    "AUA" => [ 'ILE', 'I' ],
    "ACA" => [ 'THR', 'T' ],
    "AAA" => [ 'LYS', 'K' ],
    "AGA" => [ 'ARG', 'R' ],
    "AUG" => [ 'MET', 'M' ],
    "ACG" => [ 'THR', 'T' ],
    "AAG" => [ 'LYS', 'K' ],
    "AGG" => [ 'ARG', 'R' ],
    "GUU" => [ 'VAL', 'V' ],
    "GUC" => [ 'VAL', 'V' ],
    "GCU" => [ 'ALA', 'A' ],
    "GCC" => [ 'ALA', 'A' ],
    "GAU" => [ 'ASP', 'D' ],
    "GAC" => [ 'ASP', 'D' ],
    "GGU" => [ 'GLY', 'G' ],
    "GGC" => [ 'GLY', 'G' ],
    "GUA" => [ 'VAL', 'V' ],
    "GUG" => [ 'VAL', 'V' ],
    "GCA" => [ 'ALA', 'A' ],
    "GCG" => [ 'ALA', 'A' ],
    "GAA" => [ 'GLU', 'E' ],
    "GAG" => [ 'GLU', 'E' ],
    "GGA" => [ 'GLY', 'G' ],
    "GGG" => [ 'GLY', 'G' ]);

my @seq_colors = ();

our $debug = 1;

1;

=head1 NAME

Bio::Frescobi::Sequtil -- Module containing sequence manipulations functions.

=head1 SYNOPSIS

 use Bio::Frescobi::Sequtil;

=head2 Functions

 &seq_checksum
 &get_checksum_size
 &set_checksum_size
 &format_seq
 &get_format_seq_size
 &set_format_seq_size
 &format_qual
 &get_format_qual_size
 &set_format_qual_size
 &get_checkseq
 &set_checkseq
 &qual_to_array
 &array_to_qual
 &read_1_fasta_sequence	     
 &extract_1_sequence	     
 &read_1_fastq_sequence
 &format_fastq
 &reverse_dna_seq
 &substr_dna_seq
 &is_totally_masked
 &packseq
 &unpackseq
 &prot_3_to_1
 &prot_1_to_3
 &colorize_seq
 &degenerate_index
 &fastq_qual_to_array
 &fastq_array_to_qual

=head1 METHODS

=over 4

=cut
   

=item C<seq_checksum($seq[, $size])>

    Computes a checksum for C<$seq> using the Secure Hashing Algorithm at
    a length of 40 characters (160 bits).  The checksum is truncated to a
    length of C<$size>. If not the $size variable is not provided, a
    default size, stored within this module, is used. The functions,
    C<&get_checksum_size> and C<&set_check_size>, can be used for
    manipulating this default. The initial default size is 8.
    
=cut

sub seq_checksum {
    my ($seq, $size) = @_;
    my ($checksum, $sha);

    $sha = new SHA;
    ($checksum = $sha->hexhash(uc($seq))) =~ s/ //g;
    if (not defined $size) {
	$size = $checksum_size;
    }
    else {
	# Only test if the user passes something in.
	if ($size < 1 or $size > 40) {
	    croak "Bad checksum size ($size) used for seq_checksum";
	}
    }
    return substr($checksum, 0, $size);
}

=item C<get_checksum_size>

Return the current default checksum size.

=cut

sub get_checksum_size {

    return $checksum_size;
}


=item C<set_checksum_size($size)>

Set the current default checksum size to C<$size>.

=cut

sub set_checksum_size {
    my $size = $_[0];
    my $st;

    if (not defined $size or $size < 1 or $size > 40) {
	$st = "($size) " if defined $size;
	croak "set_checksum_size must be called with checksum size value ${st}between 1 and 40";
    }
    $checksum_size = $size;
}


sub format_seq {
    my ($seq, $size) = @_;
    my ($i);

=item C<format_seq($seq[, $size])>

Accept a C<$seq> and  return a multiline string in which the sequence is broken
into more manageable size lines. The size of the output lines is given by C<$size>.
If no C<$size> is specified, then a default value is used. The
functions, C<&set_format_seq_size> and C<&get_format_seq_size>, are used to
manipulate the default size value.

=cut    

    if (not defined $size) {
	$size = $format_seq_size;
    }
    else {
        if ($size < 10) {
	    carp "Size value ($size) for format_seq must be at least 10";
	    $size = 10;
	}
    }
    my @pieces = ();
    for ($i = 0; $i < length($seq); $i+= $size) {
	push (@pieces, sprintf("%s\n", substr($seq, $i, $size)));
    }
    return join("", @pieces);
}

sub get_format_seq_size 
{

=item C<get_format_seq_size>

Return the current default C<format_seq> line size.

=cut

    return $format_seq_size;
}

sub set_format_seq_size {
    my $size = $_[0];

=item C<set_format_seq_size($size)>

Set the current default C<format_seq> size to C<$size>.

=cut

    my $st;

    if (not defined $size or $size < 1) {
	$st = "($size) " if defined $size;
	croak "set_format_seq_size must be called with checksum size value ${st}at least 10";
    }
    $format_seq_size = $size;
}

sub get_checkseq {
    
=item C<get_checkseq>

Return the current default C<checkseq> setting.

=cut

    return $checkseq;
}

sub set_checkseq {
    my $mode= $_[0];

=item C<set_checkseq($mode)>

Set the current default C<checkseq> setting to C<$mode>.

=cut


    if (not defined $mode) {
        croak "set_checkseq must be called with checksum mode value.\n";
    }
    $checkseq = $mode;
}

sub format_qual {
    my ($qual, $size) = @_;

=item C<format_qual($qual[, $size])>

Input a sequence of quality scores in Phred format (blank separated integers)
and return a multiline string in which the quality scores
are broken into smaller lines.
The size of the output lines is given by C<$size>.
If no C<$size> is specified, then a default value is used. The
functions, C<&set_format_qual_size> and C<&get_format_qual_size>, are used to
manipulate the default size value.

=cut    

    my ($i, $l);

    if (not defined $size) {
	$size = $format_qual_size;
    }
    else {
        if ($size < 10) {
	    carp "Size value ($size) for format_qual must be at least 10";
	    $size = 10;
	}
    }

    my @pieces = ();
    $i = 0;
    while ($i < length($qual)) {
	$l = $size;
	while ($i < length($qual) and substr($qual, $i, 1) eq " ") {
	    $i += 1;
	}
	last if $i >= length($qual);
	while ($i + $l - 1 < length($qual) and substr($qual, $i + $l, 1) ne " ") {
	    $l += 1;
	}
	push(@pieces, sprintf("%s\n", substr($qual, $i, $l)));
	$i += $l;
    }
    return join("", @pieces);
}

sub get_format_qual_size {
    
=item C<get_format_qual_size>

Return the current default C<format_qual> line size.

=cut

    return $format_qual_size;
}

sub set_format_qual_size {
    my $size = $_[0];

=item C<set_format_qual_size($size)>

Set the current default C<format_qual> size to C<$size>.

=cut

    my $st;

    if (not defined $size or $size < 1) {
	$st = "($size) " if defined $size;
	croak "set_format_qual_size must be called with checksum size value ${st}at least 10";
    }
    $format_qual_size = $size;
}

sub qual_to_array {
    my ($qual) = @_;

=item C<qual_to_array($qual)>

Convert a string of Phred style quality scores into an array of
integers, and return the array.

=cut
    
    my (@ret);
    local ($_);

    $qual =~ s/^\s+//;
    $qual =~ s/\s+$//;
    @ret = split(/\s+/, $qual);
    foreach (@ret) {
	$_ += 0;
    }
    return @ret;
}

sub array_to_qual {

=item C<array_to_qual($quality1, $quality2, ...)>

Convert an array of quality scores into a string with each value is
separated by a blank. No checks on the quality scores are made.

=cut

    my ($ret, $st);
    local ($_);

    return join(" ", @_);
}

sub read_1_fasta_sequence {

=item C<read_1_fasta_sequence($fh, $line[, $is_qual])>

Function for reading sequences or quality scores from a FASTA
file. The use of this function is somewhat complicated because of the
need of the calling program to maintain the last line read. The
filehandle for access to the FASTA file is passed via the C<$fh>
parameter. The last line read is passed via the C<$line> parameter,
and must be initialized to the first line of the file prior to the
first call of this function. This variable is modified in the calling
program. The optional C<$is_qual> parameter specifies whether the
FASTA file contains quality scores.

This script will work with DOS formatted files, i.e., it will handle
the CR LF line terminator properly.

On return, this function returns an array of three values, C<($name,
$annotation, $seq)>.  The C<$name> is the word following the E<gt>
character in the FASTA header line.  If the FASTA file is at end of
file, then C<$name> will be returned with an undefined value.  The
C<$annotation> will be returned with the remaining text on the header
line after all leading and trailing white space removed, tabs changed
to spaces, and multiple spaces converted to single spaces.  The
C<$seq> will be returned with the assembled sequence.  For regular
sequences, all white space will be removed prior to the return and the
sequence will be checked for non-alphanumeric characters being
included (if the C<checkseq> mode is turned on).
For quality score sequences, all tabs will be converted to
spaces and multiple spaces will be converted to single spaces, and
leading and trailing blanks will be removed.

The following code segment illustrates the use of this function:

    open(SEQ, "<$seq_file") || die "Unable to open $seq_file: $!\n";
    my $seq_line = <SEQ>;
    
    while (defined $seq_line) {
	my ($name, $annotation, $seq) =
            read_1_fasta_sequence(\*SEQ, $seq_line, 0);
	last if not defined $name;

        # process sequence
    }

=cut

    my ($fh, $line, $is_qual) = @_;
    my ($name, $annotation, $sequence, $st);
    local ($_);

    if (not defined $is_qual) {
	$is_qual = 0;
    }
    $_ = $line;
    $name = undef;
    $annotation = undef;
    $sequence = undef;
    &dos_chomp($_);
    while (defined $_ and /^\s*$/) {
	$_ = <$fh>;
    }
    if (substr($_, 0, 1) ne ">") {
	warn "Sequence start missing. Line:\n$_";
	while ($_ = <$fh>) {
	    last if substr($_, 0, 1) eq ">";
	}
    }
    else {
	if (/^>(\S+)\s+(.*)/) { # /
	    ($name, $annotation) = ($1, $2);
	}
	else {
	    $name = substr($_, 1);
	    $annotation = "";
	}
    }
    my @pieces = ();
    while ($_ = <$fh>) {
	last if substr($_, 0, 1) eq ">";
	dos_chomp($_);
	if ($is_qual) {
	    push(@pieces, " " . $_ . " ");
	}
	else {
	    s/\r//g;
	    s/\s+//g;
	    if ($checkseq) {
		if (m/[^A-Za-z0-9*]/) {
		    $st = $_;
		    $st =~ s/[A-Za-z0-9*]+//g;
		    carp "Sequence named $name has non alphanumeric characters: $st\n";
		    s/[^A-Za-z0-9*]//g;
		}
	    }
	    push (@pieces, $_);
	}
    }
    $sequence = join("", @pieces);
    &squeeze_blanks($annotation);
    &squeeze_blanks($sequence) if $is_qual;
    $_[1] = $_;
    return ($name, $annotation, $sequence);
}

sub extract_1_sequence {

=item C<extract_1_sequence($file, $id[, $is_qual])>

Function for reading one specific sequence or quality score from a
FASTA file. This function uses C<read_1_fasta_sequence> in this
library.  The C<$file> argument specifies the FASTA file name, and the
C<$id> argument specifies the identifier for the sequence being
extracted.  The optional C<$is_qual> parameter specifies whether the
FASTA file contains quality scores.

N.B.: This program does a sequential scan of the file, and therefore,
it is only meant for small sequence files.

On return, this function returns an array of two values,
C<($annotation, $seq)>.  The C<$annotation> will be returned with the
remaining text on the header line after all leading and trailing white
space removed, tabs changed to spaces, and multiple spaces converted
to single spaces.  The C<$seq> will be returned with the assembled
sequence.  For regular sequences, all white space will be removed
prior to the return and the sequence will be checked for
non-alphanumeric characters being included (provided the C<checkseq> mode is
turned on).  For quality score sequences, all tabs will be converted
to spaces and multiple spaces will be converted to single spaces, and
leading and trailing blanks will be removed. If the sequence with the
searched for C<$id> is not found, then both elements in the array will
be returned as C<undef>.
;

=cut

    my ($file, $id, $is_qual) = @_;
    my ($ret_annotation, $ret_seq);
    local (*FILE);

    if (not open(FILE, "<$file")) {
	carp "Unable to open $file: $!\n";
    }
    else {
	my $line = <FILE>;
	while (defined $line) {
	    my ($name, $annotation, $sequence) = read_1_fasta_sequence(\*FILE, $line, $is_qual);
	    last if not defined $name;
	    if ($name eq $id) {
		$ret_annotation = $annotation;
		$ret_seq = $sequence;
		last;
	    }
	}
	close FILE;
    }
    return ($ret_annotation, $ret_seq);
}

sub read_1_fastq_sequence {

=item C<read_1_fastq_sequence($fh[, $offset])>

Function for reading sequences and quality scores from a FASTQ file.
The filehandle for access to the FASTQ file is passed via the C<$fh>
parameter. Unlinke the C<read_1_fasta_sequence> function, this
function does not need to store the previous line. Just check to see
if any of the returned parameters are defined or not.

On return, this function returns an array of four values, C<($name,
$annotation, $seq, $qual)>.  If the FASTQ file is at end of
file, then C<$name> will be returned with an undefined value.  The
C<$annotation> will be returned with the remaining text on the header
line after all leading and trailing white space removed, tabs changed
to spaces, and multiple spaces converted to single spaces.  The
C<$seq> will be returned with the assembled sequence.

The sequence will be checked for non-alphanumeric characters being
included (if the C<checkseq> mode is turned on).

The quality score will just be returned as is. The C<offset> is
currently ignored.

The following code segment illustrates the use of this function:

    open(SEQ, "<$seq_file") || die "Unable to open $seq_file: $!\n";
    while (my ($name, $annotation, $seq, $qual) =
	   read_1_fastq_sequence(\*SEQ, 'b')) {

        # process sequence
    }

=cut

    my ($fh, $offset) = @_;
    my ($name, $annotation, $sequence, $qual);

    if (not defined $offset) {
	$offset = 0;
    }

    my @lines;
    for (my $i = 0; $i < 4; $i++) {
        $lines[$i] = <$fh>;
	if (not defined($lines[$i])) {
	    return ();
	}
        &dos_chomp($lines[$i]);
    }
    if (substr($lines[0], 0, 1) ne '@') {
	croak "First sequence line missing \@ character.\n";
    }
    if (substr($lines[2], 0, 1) ne '+') {
	croak "Third sequence line missing + character.\n";
    }
    ($name, $annotation) = split(/\s+/, substr($lines[0], 1), 2);
    $sequence = $lines[1];
    $qual = $lines[3];
    if ($checkseq) {
	if ($sequence =~ m/[^A-Za-z0-9]/) {
	    carp "Sequence named $name has non alphanumeric characters: $sequence\n";
	    $sequence =~ s/[^A-Za-z0-9]/N/g;
	}
    }
    return ($name, $annotation, $sequence, $qual);
}

=item C<read_fasta_file($file)>

Function for all the sequences in a
FASTA file. This function uses C<read_1_fasta_sequence> in this
library.  The C<$file> argument specifies the FASTA file name. The return is a hash
with the sequence name being the key and the sequence being the value.

The file must have unique names for each sequence.

=cut

sub read_fasta_file {

    my $file = shift;

    my %ret;
    local (*FILE);

    if (not open(FILE, "<$file")) {
	carp "Unable to open $file: $!\n";
    }
    else {
	my $line = <FILE>;
	while (defined $line) {
	    my ($name, $annotation, $sequence) = read_1_fasta_sequence(\*FILE, $line, 0);
	    last if not defined $name;
	    if (exists($ret{$name})) {
		print STDERR "Duplicate sequence name ($name) found in $file\n";
		print STDERR "Only the first sequence will be kept.\n";
	    }
	    else {
		$ret{$name} = $sequence;
	    }
	}
	close FILE;
    }
    return %ret;
}

sub format_fastq {

=item C<format_fastq($name, $annotation, $seq, $qual)>

Generate an output string containing four lines of FASTQ output.

=cut

    my ($name, $annotation, $seq, $qual) = @_;
    $annotation = "" if not defined $annotation;

    if (length($seq) != length($qual)) {
        croak sprintf("Sequence length (%d) does not match quality length (%d)\n",
		      length($seq),
		      length($qual));
    }
    my @lines;
    push (@lines, sprintf("\@%s%s%s",
			  $name,
			  $annotation eq "" ? "" : " ",
			  $annotation));
    push (@lines,
	  $seq,
	  "+",
	  $qual);
    return join("\n", @lines) . "\n";
}

sub reverse_dna_seq {

=item C<reverse_dna_seq($seq)>

Reverse complement a DNA sequence, passed in C<$seq>, including all
ambiguity codes and return the result.

=cut

    # reverse a DNA sequence.
    my ($seq) = @_;
    my ($i, $j, $ret);

    $ret = reverse($seq);
    $ret =~ tr/ATCGNMRWSYKVHDBatcgnmrwsykvhdb/ /cd;
    $ret =~ tr/ATCGNMRWSYKVHDBatcgnmrwsykvhdb/TAGCNKYWSRMBDHVtagcnkywsrmbdhv/;
    if (index($ret, " ") >= 0) {
	my $seqcopy = $seq;
	$seqcopy =~ tr/ATCGNMRWSYKVHDBatcgnmrwsykvhdb/ /cd;
	croak "Sequence $seq had unknown codes: $seqcopy\n";
    }
    return $ret;
}

sub substr_dna_seq {

=item C<substr_dna_seq($seq, $start, $stop, $is_plus)>

Return a substring of a DNA sequence passed as C<$seq> taking into account polarity,
provided in C<$is_plus> and with position identified by C<$start> and C<$stop>.
The count of position starts at 1 for start and stop. If C<$is_plus> is false,
then the extracted subsequence is reverse complemented.

=cut

    # Avoid copying the sequence, so we just reference it when needed. For big sequences,
    # this makes a huge performance difference.
    
    my $start = $_[1]; 
    my $stop = $_[2];
    my $is_plus = $_[3];
    my $seq;

    if ($start < 1) {
	carp "Start position less than 1 (= $start)\n";
	$start = 1;
    }
    if ($stop > length($_[0])) {
	carp "Stop position past end of string (start = $start  length = " . length($_[0]) . "\n";
	$stop = length($_[0]);
    }
    $seq = substr($_[0], $start - 1, $stop - $start + 1);
    if (defined $is_plus and not $is_plus) {
	$seq = reverse_dna_seq($seq);
    }
    return $seq;
}

sub is_totally_masked {

=item C<is_totally_masked($seq)>

Returns true or false if C<$seq> is composed purely of masked out characters,
i.e., X or N.    

=cut

    my ($seq) = @_;
    return $seq =~ m/^[xXnN]*$/;
}

sub packseq { }
#     my $seq = lc(shift(@_));

# =item C<packseq($seq)>

# Given a nucleotide sequence, pack it into a binary string using 2 bits for
# each nucleotide.

# =cut

#     my ($i, $template, $bitstream, $code);

#     $bitstream = "  " x length($seq);
#     for ($i = 0; $i < length($seq); $i++) {
# 	$code = $packseq{substr($seq, $i, 1)};
# 	die "Bad sequence i= $i seq = $seq\n" if length($code) != 2;
# 	substr($bitstream, 2 * $i, 2) = $code;
#     }
#     return pack(sprintf("B%d", length($bitstream)), $bitstream);
# }


sub unpackseq { }
#     my ($bits, $length) = @_;

# =item C<$unpackseq($bits, $length)>

# Inverse of C<&packseq>. The length must be passed because there is no way
# to tell when the bitstring has ended.

# =cut

#     my ($bitstream, $i, $seq, $ch);


#     $bitstream = unpack(sprintf("B%d", $length * 2), $bits);
#     $seq = " " x $length;
#     if ($length * 2 != length($bitstream)) {
# 	die "unpackseq logic error for length\n";
#     }
#     for ($i = 0; $i < length($bitstream) - 1; $i += 2) {
# 	$ch = $unpackseq{substr($bitstream, $i, 2)};
# 	die "Bad bit pattern\n" if $ch eq "";
# 	substr($seq, $i/2, 1) = $ch;
#     }
#     return $seq;
# }

sub prot_3_to_1 {

=item C<prot_3_to_1($res1, $res2, ...)>
    
Translate an array of protein residues into a string of one letter codes.
A warning message will be generated if a residue is unidentified.

=cut

    my @args = @_;
    my @trans = ();
    my $errs = "";
    my $count = 0;

    foreach my $res (@args) {
	$count += 1;
	$res = uc($res);
	if (exists($aa3to1{$res})) {
	    push (@trans, $aa3to1{$res});
	}
	elsif ($res eq "TER") {
	    if ($count != scalar(@args)) {
		$errs .= "TER(not at end) ";
		push (@trans, "U");
	    }
	}
	else {
	    $errs .= "$res ";
	    push (@trans, "U");
	}
    }
    if ($errs) {
	carp "Unable to interpret protein residue(s): $errs\n";
    }
    return join("", @trans);
}

sub prot_1_to_3
{

# =item C<prot_1_to_3($aa_string)>
	
# Translate a string of one letter codes into an array of protein residues.
# A warning message will be generated if a residue is unidentified.
# A stop letter, the asterisk, will be translated to "STOP".
	
# =cut

    my $aa_string = shift;
    my @ret;    
    my @args = @_;
    my $errs = "";
    my $count = 0;

    for (my $i = 0; $i < length($aa_string); $i++) {
	my $aa1 = substr($aa_string, $i, 1);
	if (exists($aa1to3{$aa1})) {
	    push (@ret, $aa1to3{$aa1});
	}
	else {
	    $errs .= "Unable to translate $aa1 ";
	    push (@ret,'UNK');
	}
    }
    if ($errs) {
	carp "Unable to interpret protein residue(s): $errs\n";
    }
    return @ret;
}

sub na_2_prot {

=item C<na_2_prot($seq[, $one_letter])>

Translate a DNA or RNA sequence into a protein sequence. Returns an array.
If $one_letter is defined and true, then return one letter codes.

=cut

    my @ret = ();
    my $na_seq = uc(shift);
    $na_seq =~ tr/T/U/;
    my $one_letter = shift;
    my $lookup;
    if ($one_letter) {
	$lookup = 1;
    }
    else {
	$lookup = 0;
    }
    my $div_err = 0;
    if (length($na_seq) % 3 != 0) {
	carp "na_2_prot passed a sequence that is not divisible by 3.\n" if $debug;
	$div_err = 1;
	$na_seq .= 'UUU'; # Add extra junk that will be fixed later.
    }
    my $codon_err = 0;
    for (my $i = 0; $i < length($na_seq) - 2; $i += 3) {
	my $codon = substr($na_seq, $i, 3);
	my $trans;
	if (not exists($genetic_code{$codon})) {
	    if (not $codon_err) {
		carp "na_2_prot was passed an invalid nucleic sequence: '$codon'.\n" if $debug;
		$codon_err = 1;
	    }
	    $trans = [ '???', '?' ];
	}
	else {
	    $trans = $genetic_code{$codon};
	}
	push @ret, $trans->[$lookup];
    }
    if ($div_err) {
	pop @ret;
    }
    return @ret;
}

sub colorize_seq {

=item C<colorize_seq($seq)>

Return an HTML formatted string that will colorize the sequence, which
is presumed to be a protein sequence. The colors were tweeked to make
them all distinguishable.

=cut

    my $seq = shift;
    my @ret;

    if (scalar(@seq_colors) == 0) {
	&set_colors;
    }
    for (my $i = 0; $i < length($seq); $i++) {
	my $ch = substr($seq, $i, 1);
	my $code = ord(uc($ch)) - ord('A');
	if ($code < 0 or $code > 25) {
	    $code = 26;
	}
	push(@ret,
	     sprintf('<font color="%s">%s</font>', $seq_colors[$code], $ch));
    }
    return join("", @ret);
}

sub set_colors {
    # This algorithm was tweaked to generate at least 27 colors that
    # are all distinguishable from each other.
    my @choices = (0x0, 0x55, 0xaa, 0xff);
    my $i = 0;

    foreach my $c1 (@choices) {
	foreach my $c2 (@choices) {
	    foreach my $c3 (@choices) {
		my $sum = $c1 + $c2 + $c3;
		next if ($sum <= 0x99 or $sum >= 2 * 0xFF);
		next if $c2 == 0xff;
		push (@seq_colors, '#' . join("", map { sprintf("%02x", $_); } ($c1, $c2, $c3)));
	    }
	}
    }
    srand(1);
    @seq_colors = sort { &color_sum($a) <=> &color_sum($b) } @seq_colors;
    for (my $i = 1; $i < 1000; $i++) {
 	my $j = min(scalar(@seq_colors) - 1, int(rand(scalar(@seq_colors))));
	my $k = min(scalar(@seq_colors) - 1, int(rand(scalar(@seq_colors))));
 	my $t = $seq_colors[$j];
 	$seq_colors[$j] = $seq_colors[$k];
 	$seq_colors[$k] = $t;
    }
}


sub color_sum {
    my $code = shift;
    return hex(substr($code, 1, 2)) + hex(substr($code, 3, 2)) + hex(substr($code, 5, 2));
}

sub color_table {
    # A debugging function.
    my $q = shift;
    
    my @rows = ($q->TR($q->th("Code"),
		       $q->th("Color")));
    foreach my $color (@seq_colors) {
	push (@rows, $q->TR($q->td($color),
			    $q->td($q->font({-color => $color}, "ABCD"))));
    }
    return $q->table( {-style => 'font-family: monospace; font-size: larger',
		       -border => 1},
		      @rows);
}

sub degenerate_index {

=item C<degenerate_index($source, $target)>

A string indexing function for the source string can have an N which
is effectively a wildcard when matching the target. If there is no
match, then a -1 is returned. For examples,

 degenerate_index("ATCGATCG", "CGA") == 2
 degenerate_index("ATCGNTCG", "GCT") == 3

=cut
    
    my $source = lc(shift);
    my $target = lc(shift);

    if (length($target) > length($source)) {
	return -1;
    }
    elsif (index($source, 'n') == -1) {
	return index($source, $target);
    }
    else {
	for (my $i = 0; $i < length($source) - length($target) + 1; $i++) {
	    my $match = 1;
	    for (my $j = 0; $j < length($target); $j++) {
		my $sch = substr($source, $i + $j, 1);
		if ($sch ne 'n' and
		    $sch ne substr($target, $j, 1)) {
		    $match = 0;
		    last;
		}
	    }
	    if ($match) {
		return $i;
	    }
	}
	return -1;
    }
}


=item C<fastq_qual_to_array($qual[, $offset])>

Converts a FASTQ quality string, C<$qual>, to an array of phred
quality scores. The C<$offset> is the zero quality character, and
is optional. If not specified, then it defaults to 'B'.

=cut

sub fastq_qual_to_array {
    my $qual = shift;
    my $offset = shift;

    my $offset_i = defined($offset) ? ord($offset) : ord('B');
    my @ret = map { $_ -= $offset_i; } unpack("C*", $qual);
    return @ret;
}


=item C<fastq_array_to_qual($qual_array[, $offset])>

Converts an array of phred quality scores pointed to C<$qual_array> to
a FASTQ quality string which is returned.  The C<$offset> is the zero
quality character, and is optional. If not specified, then it defaults
to 'B'.

=cut

sub fastq_array_to_qual {
    my $qual_array = shift;
    my $offset = shift;

    my $offset_i = defined($offset) ? ord($offset) : ord('B');
    my @quals = map { $_ = chr($_ + $offset_i); } @{$qual_array};
    return join("", @quals);
}

=item C<find_ns($dna_seq)>

Searches the $dna_seq for N's and returns a array containing arrays of
size 2 with start and end positions of N's. The positions are 1-based.

=cut

sub find_ns {
    my $seq = uc(shift);
    my @ret;

    my $pos = 0;
    while ((my $ind = index($seq, 'N', $pos)) > -1) {
	my $start = $ind + 1;
	my ($n_run) = (substr($seq, $ind) =~ m/^(N+)/);
	if (not defined $n_run) {
	    print STDERR "Error in find_ns. \$n_run not defined.\n";
	    print STDERR "ind = $ind. Substr Seq = ", substr($seq, $ind, 100), "\n";
	    last;
	}
	$ind += length($n_run);
	my $stop = $ind;
	push (@ret, [ $start, $stop ]);
	last if $ind >= length($seq);
	$pos = $ind;
    }
    return @ret;
}
    
=back

=head1 SEE ALSO

Bio::Frescobi::Genutil(3)

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
