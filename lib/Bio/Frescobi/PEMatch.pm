package Bio::Frescobi::PEMatch;

use 5.008;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION $AUTOLOAD);
use warnings;
use Carp;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::Sequtil;
use Text::LevenshteinXS qw(distance);

require Exporter;

@ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration
# use Bio::Frescobi::PEMatch ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.

our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw();

our $AUTOLOAD;

my %component =
    map { ($_ => 1) }
        qw(MINOFFSET MAXOFFSET MINOVERLAP MAXD QUALBASE MATCH_WEIGHT DEBUG);

our $VERSION = '0.01';

1;

=head1 NAME

Bio::Frescobi::PEMatch - Perl extension for Paired End Matching. 

=head1 SYNOPSIS

  use Bio::Frescobi::PEMatch;
  my $pematch = new Bio::Frescobi::PEMatch;
  $pematch-><components>([new_value]);
  ($joined_seq, $joined_qual) = $pematch->match($seq1, $qual1,
                                                $seq2, $qual2);
  $pematch->equal;
  $pematch->best_offset;

=head2 Constructor

  $pematch = new Bio::Frescobi::PEMatch;

  # You can also use $obj = Bio::Frescobi::PEMatch->new();

=head2 Object Methods

=head3 Parameter Setting Methods

  $pematch->minoffset([$minoffset]);
  $pematch->maxoffset([$maxoffset]);
  $pematch->maxd([$maxd]);
  $pematch->minoverlap([$minoverlap]);
  $pematch->qualbase([$qualbase]);
  $pematch->match_weight([$match_weight]);
  $pematch->debug([$debug]);

=head3 Output Methods

  ($joined_seq, $joined_qual) = $pematch->match($seq1, $qual1,
                                                $seq2, $qual2);
  $pematch->equal;
  $pematch->best_offset;

=head1 DESCRIPTION

The Bio::Frescobi::PEMatch module provides a mechanism for joining
paired end reads together if the ends of the reads overlap. The module
will search over a range of possible overlap values looking for a good
overlap. The programmer can specify certain aspects of the search and
scoring of the overlaps based on the nature of the sequencing
project. Quality scores are also taken into account in an I<ad hoc>
way.

The search for possible overlaps begins by taking the reverse
complement of the second sequence.  The range of possible offsets is
specified by a minimum offset, and either a maximum offset or minimum
allowed overlap.

For each possible offset, the overlap region is identified. If any
unknown bases are present, i.e. 'N', then the base from the other
sequence is substituted, i.e. N's become free matches. If bases at any
position are different and the quality for a base is less than or
equal to the
lowest possible quality score (parameter C<qualbase>), then the lower
quality base is substituted with the higher quality one.

The sequences in the overlap region are compared using the Levenshtein
algorithm (see L<Text::LevenshteinXS>) and a distance (d) is
computed. The distance is the number of insertions, deletions, and
substitutions needed to transform one sequence into the other. If the
distance is greater than the C<maxd> parameter, then the method
continues with the next possible offset. Otherwise, a score for the
overlap is generated using the formula,

score = length(overlap) - match_weight x d

The C<match_weight> parameter is used to balance the length of the
overlap against the quality of the match. The length must be taken
into account, because short sequences can overlap by chance
alone. Some judgment and possible trial and error is necessary to set
the parameters.

The program identifies the offset with the highest score and uses that
offset, and the substitutions above to construct a merged sequence and
quality score that it returns. If no acceptable offset is found, then
undefined values are returned.

The most reliable setting for C<maxd> would be zero. Otherwise, the
program can get confused when the offset is close to correct, because
incorrect offset by just one base can be adjusted by an
insertion. Thus, the C<match_weight> parameter should be kept greater
than one so that the overlap dominates.

=head1 METHODS

The following methods are provided:

=over 4

=item C<new()>

Creates a new Bio::Frescobi::PEMatch object. There are no parameters to the method.

=item C<minoffset([$minoffset_value])>

Retrieves the current setting for minimum offset used in the
search. If a argument is passed to this method, then the minimum
offset is set in the object. Its effect will not occur until another
search is done. The default is 0.

=item C<maxoffset([$maxoffset_value])>

Retrieves the current setting for maximum offset used in the
search. If a argument is passed to this method, then the mamimum
offset is set in the object. Its effect will not occur until another
search is done. The default is -1.

If the maximum offset is set to a negative number, then the maximum
offset is determined by the minimum overlap parameter, see below, and
the length of the sequences being compared.

=item C<minoverlap([$minoverlap_value])>

Retrieves the current setting for minimum overlap used in the
search. If a argument is passed to this method, then the minimum
overlap is set in the object. Its effect will not occur until another
search is done. The default is 10.

This parameter has an effect only with the C<maxoffset> parameter is
negative.

=item C<maxd([$maxd_value])>

Retrieves the current setting for maximum allowed mismatch distance
for comparing two overlapped subsequences.  If a argument is passed to
this method, then the maximum allowed mismatch distance is set in the
object. Its effect will not occur until another search is done. The
default is 1.

=item C<qualbase([$qualbase_value])>

Retrieves the quality score base value.
If a argument is passed to this method, then the minimum
offset is set in the object. Its effect will not occur until another
search is done.

The quality score base value depends on the sources of the sequencing
and is expressed as a character. The binary code of the character is
the base value. The default is C<'B'>.

=item C<match_weight([$match_weight_value])>

Retrieves the current setting for the match weight used in the search,
see above. If a argument is passed to this method, then the match weight
is set in the object. Its effect will not occur until another
search is done.

=item C<debug([$debug_value])>

Retrieves the current setting for debugging output used in the
search. If a argument is passed to this method, then the debug setting
will be changed in the object. Currently, a value of 2 gives the most
information. Output is sent to C<STDERR>.

=item C<match($seq1, $qual1, $seq2, $qual2)>

The C<match> method actually does the work of matching two sequences
as described above. If the matching is successful, the method returns
a list of two values, the joined sequence and joined quality
scores. If the match is unsuccessful, then two C<'undef'> values are
returned.

=item C<best_offset>

Returns the best offset from the previous search. It is an
error to call this method before a match is performed.

=item C<equal>

Returns either the string, C<"Equal"> or C<"Unequal">, if the
overlapping sequences match exactly or not, respectively.  It is an
error to call this method before a match is performed.

=back

=head1 SEE ALSO

pematch.pl

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut


sub new {
    my $pkg;
    my $class = shift;
    eval {($pkg) = caller(0);};
    if ($class ne $pkg) {
	unshift @_, $class;
    }
    my $self = {};
    bless $self;
    $self->{MINOFFSET} = 0;
    $self->{MAXOFFSET} = -1;
    $self->{MINOVERLAP} = 10;
    $self->{MAXD} = 1;
    $self->{DEBUG} = 0;
    $self->{MATCH_WEIGHT} = 2;
    $self->{QUALBASE} = 'B';
    $self->{EQUAL} = undef;
    $self->{BESTOFFSET} = undef;
    $self->{STATE} = "setup";
    return $self;
}

sub DESTROY {
    
    my $self = shift;
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self);
    my $pkg = __PACKAGE__;
    my ($pos, $new_val);

    confess "$self is not an object. AUTOLOAD = $AUTOLOAD\n" if not $type;
    confess "$pkg AUTOLOAD function fails on $type\n"
	if $type !~ m/^$pkg$/;
    my $name = uc($AUTOLOAD);
    $name =~ s/^.*:://;
    if (not exists($component{$name})) {
	confess "$name is not a valid method for $type.\n";
    }
    elsif (not exists($self->{$name})) {
	confess "$name is not a valid method for $type\n";
    }
    elsif (defined $_[0]) {
	$self->{$name} = shift;
	$self->{STATE} = "setup";
    }
    if ($name =~ m/^(MINOFFSET|MINOVERLAP|MAXD)$/ and
	$self->{$name} < 0) {
	confess ucfirst($name) . " parameter must be non-negative.\n";
    }
    if ($name eq "QUALBASE" and
	length($self->{QUALBASE}) != 1) {
	confess "Qualbase parameter must be a single characater.\n";
    }
    if ($name eq "MATCH_WEIGHT" and
	$self->{MATCH_WEIGHT} <= 0.0) {
	confess "Match_weight parameter must be positive.\n";
    }
    return $self->{$name};
}

sub match {
    my $self = shift;
    my ($seq1, $qual1, $seq2, $qual2) = @_;

    if (not defined $seq1 or
	not defined $seq2 or
	$seq1 eq "" or
	$seq2 eq "") {
	confess "Both sequences for matching must be defined and not empty.\n";
    }
    if (length($seq1) != length($seq2)) {
	confess "Both sequences for matching must be the same length.\n";
    }
    if (not defined $qual1 or
	$qual1 eq "") {
	$qual1 = chr(ord($self->{QUALBASE}) + 30) x length($seq1);
    }
    if (not defined $qual2 or
	$qual2 eq "") {
	$qual2 = chr(ord($self->{QUALBASE}) + 30) x length($seq2);
    }
    if (length($qual1) != length($qual2)) {
	confess "Both quality strings for matching must be the same length.\n";
    }
    my $rev_seq2 = reverse_dna_seq($seq2);
    my $rev_qual2 = reverse($qual2);
    my $best_offset = -1;
    my $best_score = 0;
    my $maxoffset = $self->{MAXOFFSET};
    if ($maxoffset < 0) {
	$maxoffset = length($seq1) - $self->{MINOVERLAP};
	print STDERR "maxoffset = $maxoffset\n" if $self->{DEBUG};
    }
    for (my $offset = $self->{MINOFFSET}; $offset <= $maxoffset; $offset++) {
	my $f1 = substr($seq1, $offset);
	my $r2 = substr($rev_seq2, 0, length($rev_seq2) - $offset);
	my $f1qual = substr($qual1, $offset);
	my $r2qual = substr($rev_qual2, 0, length($rev_qual2) - $offset);
	for (my $i = 0; $i < length($f1); $i++) {
	    my $ch1 = substr($f1, $i, 1);
	    my $ch2 = substr($r2, $i, 1);
	    my $q1 = substr($f1qual, $i, 1);
	    my $q2 = substr($r2qual, $i, 1);
	    if ($ch1 eq 'N' xor $ch2 eq 'N') {
		if ($ch1 eq 'N') {
		    substr($f1, $i, 1) = $ch2;
		}
		else {
		    substr($r2, $i, 1) = $ch1;
		}
	    }
	    elsif (($q1 le $self->{QUALBASE} xor $q2 le $self->{QUALBASE}) and
		   $ch1 ne $ch2) {
		if ($q1 le $self->{QUALBASE}) {
		    substr($f1, $i, 1) = $ch2;
		}
		else {
		    substr($r2, $i, 1) = $ch1;
		}
	    }   
	}
	my $d = distance($f1, $r2);
	if ($self->{DEBUG} > 1) {
	    print STDERR "Matching (d=$d) :\n";
	    print STDERR "$f1\n$r2\n";
	}
	if ($d <= $self->{MAXD}) {
	    # The MATCH_WEIGHT parameter avoids a longer overlap counting
	    # more than a shorter exact match.
	    my $score = length($f1) - $self->{MATCH_WEIGHT} * $d;
	    print STDERR "Score = $score\n" if $self->{DEBUG} > 1;
	    if ($score > $best_score) {
		if ($self->{DEBUG} > 1) {
		    print STDERR " Updating best_score from $best_score to $score. Offset = $offset\n";
		}
		$best_score = $score;
		$best_offset = $offset;
	    }
	}
    }
    my @ret;
    $self->{BEST_OFFSET} = $best_offset;
    if ($best_offset == -1) {
	push(@ret, undef, undef);
    }
    else {
	my $f1 = substr($seq1, $best_offset);
        my $r2 = substr($rev_seq2, 0, length($rev_seq2) - $best_offset);
	if ($f1 eq $r2) {
	    if ($best_offset == 0) {
		@ret = ($seq1, &_best_qual($qual1, $qual2));
		$self->{EQUAL} = "Equal";
	    }
	    else {
		my $join_seq = substr($seq1, 0, $best_offset) .
		    $f1 . substr($rev_seq2, -$best_offset);
		my $join_qual = substr($qual1, 0, $best_offset) .
		    &_best_qual(substr($qual1, $best_offset),
				substr($rev_qual2, 0, length($rev_qual2) - $best_offset)) .
		    substr($rev_qual2, -$best_offset);
		@ret = ($join_seq, $join_qual);
		$self->{EQUAL} = "Equal";
	    }
	}
	else {
	    my $f1qual = substr($qual1, $best_offset);
	    my $r2qual = substr($rev_qual2, 0, length($rev_qual2) - $best_offset);
	    my $join_seq = substr($seq1, 0, $best_offset);
	    my $join_qual = substr($qual1, 0, $best_offset);
	    for (my $i = 0; $i < length($f1); $i++) {
		my $base1 = substr($f1, $i, 1);
		my $base2 = substr($r2, $i, 1);
		my $q1 = ord(substr($f1qual, $i, 1));
		my $q2 = ord(substr($r2qual, $i, 1));
		if ($base1 eq $base2) {
		    $join_seq .= $base1;
		    $join_qual .= chr(max($q1, $q2));
		}
		else {
		    $join_qual .= chr(int(($q1 + $q2) / 2));
		    if ($q2 > $q1) {
			$join_seq .= $base2;
		    }
		    else {
			$join_seq .= $base1;
		    }
		}
	    }
	    $join_seq .= substr($rev_seq2, -$best_offset);
	    $join_qual .= substr($rev_qual2, -$best_offset);
	    @ret = ($join_seq, $join_qual);
	    $self->{EQUAL} = "Unequal";
	}
    }
    $self->{STATE} = "computed";
    return @ret;
}

sub _best_qual {
    my ($qual1, $qual2) = @_;
    die "Qual length mismatch in best_qual:\n $qual1\n $qual2\n"
	if length($qual1) != length($qual2);
    my $ret = " " x length($qual1);
    for (my $i = 0; $i < length($qual1); $i++) {
	my $q1 = substr($qual1, $i, 1); 
	my $q2 = substr($qual2, $i, 1); 
	if (ord($q1) > ord($q2)) {
	    substr($ret, $i, 1) = $q1;
	}
	else {
	    substr($ret, $i, 1) = $q2;
	}
    }
    return $ret;
}

sub best_offset {
    my $self = shift;
    if ($self->{STATE} eq 'setup') {
	croak "best_offset called before any matching was done.\n";
    }
    return $self->{BEST_OFFSET};
}

sub equal {
    my $self = shift;
    if ($self->{STATE} eq 'setup') {
	croak "best_offset called before any matching was done.\n";
    }
    return $self->{EQUAL};
}
