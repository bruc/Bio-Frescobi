# This package contains general standalone utilities.

# Copyright 2010 Congenomics, LLC
# This package is free software; you can redistribute it and/or modify
# it under the same terms as Perl itself, either Perl version 5.8.8 or,
# at your option, any later version of Perl 5 you may have available.

package Bio::Frescobi::Genutil;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Carp;
use POSIX;
use Digest::MD5;
@ISA = ('Exporter');
@EXPORT = qw(&system_with_check
	     &system_with_echo
	     &file2array
	     &file2string
	     &scanarray
	     &incrarray
	     &dated_mesg
	     &set_line_buffering
	     &max &min &sum &abs
	     &maxst &minst
	     &piper
	     &trim &squeeze_blanks
	     &permute
	     &add_to_where
	     &sprint1f_undef
	     &string_truth
	     &lc_hash_keys
	     &bool2logical
	     &logical2bool
	     &quotify
	     &find_smallest
	     &commify
	     &date
	     &clean_undef
	     &mark_undef
	     &find_num_in_array
	     &find_string_in_array
	     &dos_chomp
	     &md5sum
             &parse_ms_config
             &display_variables
	     &display_environment
             &random_string
             &range_overlap);

$VERSION = '0.08';

#Guarantee reproducibility unless the programmer calls srand themselves.

BEGIN {
    srand(1);
}


=head1 NAME

Bio::Frescobi - Free Relational Eclectic System for the Collection Of Biological Information
    
Bio::Frescobi::Genutil - General utility functions    

=head1 SYNOPSIS

 use Bio::Frescobi::Genutil;

=head2 Functions

 system_with_check
 system_with_echo
 file2array
 file2string
 scanarray
 incrarray
 dated_mesg
 set_line_buffering
 max
 min
 abs
 sum
 maxst
 minst
 piper
 trim
 squeeze_blanks
 permute
 add_to_where
 sprint1f_undef
 string_truth
 lc_hash_keys
 bool2logical
 logical2bool
 quotify
 find_smallest
 commify
 date
 clean_undef
 mark_undef
 dos_chomp
 md5sum
 parse_ms_config
 display_variables
 display_environment
 random_string
 range_overlap
    
=head1 DESCRIPTION

The C<Bio::Frescobi::Genutil> modules provides a number of general utility
functions for Perl programming. 

=head2 Function Descriptions

=over 4

=cut

sub system_with_check {

=item C<< system_with_check($command[, $echo]) >>

Issue a system command, but check the result. If an error is detected,
then the program croaks. If the $echo parameter is specified and is
true, then the command is echoed to STDERR along with the time prior to execution.

=cut
    
    my($command, $echo) = @_;
    my($status);

    if (defined $echo and $echo) {
	&dated_mesg($command);
    }
    $status = system("$command");
    if ($status != 0) {
	croak sprintf("Shell command: $command returned $status error code. Current directory is %s\n",
		      `pwd`);
    }
}

sub system_with_echo {

=item C<system_with_echo($command[, $echo])>

Issue a system command with optional echoing. If the $echo parameter is specified and is
true, then the command is echoed to STDERR along with the time prior to execution.

=cut

    my($command, $echo) = @_;
    my($status);

    if (defined $echo and $echo) {
	&dated_mesg($command);
    }
    $status = system("$command");
    return $status;
}

sub file2array {

=item C<file2array($file)>

Read a file and return an array containing all the lines. Each of the
lines will be C<chomp>ed, so this
function is sensitive to the input record separator, C<$/>.

=cut

    my ($file) = @_;
    local (*IN, $_);
    my (@ret);

    open (IN, "<$file") || croak "Unable to open $file: $!";
    @ret = ();
    while (<IN>) {
	chomp;
	push (@ret, $_);
    }
    close IN;
    return @ret;
}

sub file2string {

=item C<file2string($file)>

Read a file and return a string containing all the lines concatenated together.

=cut

    my ($file) = @_;
    local (*IN, $_);
    my @pieces = ();
    open (IN, "<$file") || croak "Unable to open $file: $!";
    while (<IN>) {
	push(@pieces, $_);
    }
    close IN;
    return join("", @pieces);
}

sub scanarray {

=item C<scanarray($arrayref, $cursor, $pat, $croak_if_not_found)>

Scan an array, referenced by C<$arrayref>, for the next occurrence of
a pattern, C<$pat>. This function is
intended for use in a context where the programmer wishes to scan for
all occurrences of a pattern and take action knowing the index of the
array element which contains the string.

Prior to the first call of this function, the variable, C<$cursor>,
should be set to 0, or to the first subscript where you wish to start
the scan. This variable will be modified by this function to the
position of the next element which contains the pattern.

If the pattern is not found and C<$croak_if_not_found> is true, then
the function will C<croak>. Otherwise, C<$cursor> will be set to -1.

=cut

    
    my ($arrayref, $cursor, $pat, $croak_if_not_found) = @_;
    my ($i, $len);

    croak "arrayref argument is not an array\n" if not ref($arrayref) eq "ARRAY";
    $len = scalar(@{$arrayref});
    if ($cursor < 0 || $cursor >= $len) {
	croak sprintf("Cursor argument ($cursor) is out of bounds (0..%d)\n",
		      $len - 1);
    }
    for ($i = $cursor; $i < $len; $i++) {
	if ($arrayref->[$i] =~ m/$pat/) {
	    $_[1] = $i;
	    return;
	}
    }
    if (defined $croak_if_not_found and
	$croak_if_not_found) {
	croak "Unable to find $pat\n";
    }
    $_[1] = -1;
    return;
}

sub incrarray {

=item C<incrarray($arrayref, $cursor, $increment)>

Increment the C<$cursor> variable by the C<$increment> for the array
referenced by C<$arrayref>. If the C<$cursor> is bumped beyond the end of
the array, then C<croak>.

=cut

    my ($arrayref, $cursor, $increment) = @_;
    my $len = scalar(@{$arrayref});

    $cursor += $increment;
    croak "incrarray: cursor ($cursor) off end\n" if $cursor >= $len or $cursor < 0;
    $_[1] = $cursor;
}
 
sub dated_mesg {

=item C<dated_mesg($message)>

Print onto C<STDERR> the C<$message>, with a date string at the
beginning. A newline is added to the message if not already present.

=cut

    # Print a dated message to STDERR
    my ($mesg) = @_;
    local $/ = "\n";

    my $date = &date;
    if (substr($mesg, -1, 1) ne "\n") {
	$mesg .= "\n";
    }
    print STDERR "At $date: $mesg";
}
     
sub set_line_buffering {

=item C<set_line_buffering>

Turn on line buffering for C<STDERR> and C<STDOUT>.

=cut
    
    select((select(STDOUT), $| = 1)[$[]);
    select((select(STDERR), $| = 1)[$[]);
}

sub max {

=item C<max($arg1[, $arg2]...)>

Find the numerical maximum of all the arguments and return it. If
there are no arguments, C<undef> will be returned. All arguments must
be numeric.

=cut

    my ($ret);
    local ($_);

    return undef if scalar(@_) == 0;
    $ret = shift(@_);
    foreach (@_) {
	if ($_ > $ret) {
	    $ret = $_;
	}
    }
    return $ret;
}

sub min {

=item C<min($arg1[, $arg2]...)>

Find the numerical minimum of all the arguments and return it. If
there are no arguments, C<undef> will be returned.

=cut

    my ($ret);
    local ($_);

    return undef if scalar(@_) == 0;
    $ret = shift(@_);
    foreach (@_) {
	if ($_ < $ret) {
	    $ret = $_;
	}
    }
    return $ret;
}

sub abs {

=item C<abs($arg)>

Find the absolute value of the argument, C<$arg>.

=cut

    my ($arg) = @_;
    if ($arg < 0.0) { $arg = -$arg; }
    return $arg;
}

sub sum {

=item C<sum($arg1[, $arg2]...)>

Find the sum of all the arguments and return it. If
there are no arguments, C<undef> will be returned.

=cut

    my $ret = 0;
    local ($_);

    return undef if scalar(@_) == 0;
    foreach (@_) {
	$ret += $_;
    }
    return $ret;
}

sub minst {

=item C<minst($arg1[, $arg2]...)>

Find the string minimum of all the arguments and return it. If
there are no arguments, C<undef> will be returned.

=cut

    my (@args) = @_;
    my ($i, $minst);

    confess "Zero arguments passed to minst" if $#args < 0;
    $minst = $args[0];
    for ($i = 1; $i <= $#args; $i++) {
	if ($args[$i] lt $minst) {
	    $minst = $args[$i];
	}
    }
    return $minst;
}

sub maxst {

=item C<maxst($arg1[, $arg2]...)>

Find the string maximum of all the arguments and return it. If
there are no arguments, C<undef> will be returned.

=cut

    my (@args) = @_;
    my ($i, $maxst);

    confess "Zero arguments passed to maxst\n" if $#args < 0;
    $maxst = $args[0];
    for ($i = 1; $i <= $#args; $i++) {
	if ($args[$i] gt $maxst) {
	    $maxst = $args[$i];
	}
    }
    return $maxst;
}

    
sub piper {

=item C<piper($cmd, $lines_array_ref)>

Runs a command in double ended pipe.
Two arguments must be provided:

=over 2

=item 1

a string containing the command to be executed.

=item 2

an array reference containing the lines to be passed  to the command.

=back

An array is returned. The first element is a number
containing the status encoded as follows:

=over 2

0: normal execution

>0: exit status of command.

<0: exec failure. In this case, the second element of the array
will contain the error message.

=back

If the status >= 0, the remaining array elements contain the
output of the command including errors.

=cut

    my ($command, $lineref) = @_;
    my ($pid, $waitpid);
    my ($line, $status);
    my (@log);
    local (*PIPER_IN, *PIPER_OUT, *PIPER_COM);

#   First, setup a pipe and turn off buffering.

    pipe(PIPER_IN, PIPER_OUT);
    select((select(PIPER_IN), $| = 1)[$[]);
    select((select(PIPER_OUT), $| = 1)[$[]);

#   Now fork. The child handles passing data to the command
#   and parent process reads the results.

    if (($pid = fork()) != 0) {
	# parent reads.
	close(PIPER_OUT);
	while (<PIPER_IN>) {
	    push (@log, $_);
	}
	close PIPER_IN;
	$status = waitpid($pid, 0);
	if ($status == -1) {
	    $status = 0;
	}
	else {
	    $status = $? >> 8;
	}
	unshift (@log, $status);
    }
    elsif (defined $pid) {
	# Here the child executes. STDOUT and STDERR are connected
	# to the pipe, and then we run the command
	# using inputs from the user.
	# Standard input is used to avoid any problems with quoting
	# parameters to the shell.
	close(PIPER_IN);
	select((select(STDOUT), $| = 1)[$[]);
	select((select(STDERR), $| = 1)[$[]);
	open (STDOUT, ">&PIPER_OUT");
	open (STDERR, ">&PIPER_OUT");
	if (! open(PIPER_COM, "| $command")) {
	    exit 5;
	}
	foreach my $line (@$lineref) {
	    chomp $line;
	    print PIPER_COM $line, "\n";
	}
	close PIPER_COM;
	$status = $? >> 8;
	close PIPER_OUT;
	close STDOUT;
	exit $status;
    }
    else {
	@log = (-1, "Unable to exec. Reason: " . $!);
    }
    return @log;
}

sub trim {

=item C<trim($st)>

Remove leading and trailing white space from C<$st>. The argument,
C<$st>, is modified.

=cut
    
    $_[0] =~ s/^\s+//;
    $_[0] =~ s/\s+$//;
    return $_[0];
}

sub squeeze_blanks {

=item C<squeeze_blanks($st)>

Remove leading and trailing white space from C<$st> and condense all
white space to single blanks. The argument, C<$st>, is modified.

=cut

    trim($_[0]);
    $_[0] =~ s/\t/ /g;
    $_[0] =~ s/\s{2,}/ /g;
}

sub permute {

=item C<permute($alphabet, $splitfactor)>

Generates all possible permutations of the C<$alphabet> with a length
of C<$splitfactor>. They are returned as an array of strings.

This function is implemented recursively.

=cut

    my ($alphabet, $splitfactor) = @_;
    my ($ch, $st, @ret);

    @ret = ();
    if ($splitfactor == 0) {
	return @ret;
    }
    elsif ($splitfactor == 1) {
	@ret = split(//, $alphabet);
	return @ret;
    }
    else {
	for $ch (split(//, $alphabet)) {
	    foreach $st (&permute($alphabet, $splitfactor - 1)) {
		push (@ret, $ch . $st);
	    }
	}
	return @ret;
    }
}

sub add_to_where {

=item C<add_to_where($where, $addition, $conjunction)>

A convenience function for SQL command generation, specifically C<'WHERE'> clauses.

The first argument is the where clause under construction. It should
be initialized to a blank string. This argument is modified by this
function.

The second argument is the additional expression, which is intended to
be a primary expression, e.g. C<"seqid = 'N3_5'">. The third argument
is the conjunction to be used, typically, C<AND>.

The function checks to see if the where clause is non-blank. If so, it
adds the conjunction.

Next, it adds the additional expression to the where clause and returns.

=cut

    my ($where, $addition, $conjunction) = @_;

    return if $addition eq "";
    if ($where eq "") {
	$_[0] = $addition;
    }
    else {
	$_[0] = $where . " " . $conjunction . " " . $addition;
    }
}

sub sprint1f_undef {

=item C<sprint1f_undef($fmt, $val)>

This function is a convenience function for handling the case of
printing a scalar variable which might not be defined. It accepts two
arguments, a format suitable for C<sprintf> for displaying the variable
and a possibly undefined variable.

If the value is undefined, then the function retrieves the length of
the format by looking for a number preceded by a percent sign, and
returns a blank string of that length. If no such number is found,
then a single blank is returned. Otherwise, it returns the value of
C<sprintf($fmt, $val)>.

=cut

    my ($fmt, $val) = @_;
    my ($size) = ($fmt =~ m/%(\d+)/);
    $size = 1 if not defined $size;
    if (defined $val) {
	return sprintf($fmt, $val);
    }
    else {
	return sprintf("%${size}s", "");
    }
}

sub string_truth {
    my $arg = shift;
    my $default = shift;
    my $name = shift;

=item C<string_truth($arg[, $default[, $name]])>

Converts yes, no, on, off, 0, or 1 to
0 or 1 truth value. If the argument is not defined, then the default argument
is used. If no default argument is specified, then 0 is returned.
The name of the option can optionally be passed in as the third argument.

=cut

    $default = 0 if not defined $default;
    if (not defined $arg) {
	return $default;
    }
    elsif ($arg =~ m/^(yes|on|1)$/i) {
	return 1;
    }
    elsif ($arg =~ m/^(no|off|0)$/i) {
	return 0;
    }
    else {
	if (defined $name) {
	    carp "Bad $name option: $arg\n";
	}
	else {
	    carp "Bad option: $arg\n";
	}
	return undef;
    }
}

sub lc_hash_keys {
    my %hash = ();

=item C<lc_hash_keys(%hash)>

Returns a hash where all the keys have been converted to lower case.

=cut

    my ($key, $value);

    $key = shift;
    while (defined $key) {
	$hash{lc($key)} = shift;
	$key = shift;
    }
    return %hash;
}

sub bool2logical {

=item C<bool2logical($arg)>

Converts a PostgreSQL boolean (C<t> or C<f>) to a Perl logical value
(1 or 0).

=cut

    my ($arg) = @_;

    if (not defined $arg) {
	return 0;
    }
    elsif ($arg eq "t" or $arg eq 1) {
	return 1;
    }
    elsif ($arg eq "f" or $arg eq 0) {
	return 0;
    }
    else {
	carp "Unexpected boolean argument ($arg). Returning false.\n";
	return 0;
    }
}

sub logical2bool {

=item C<logical2bool($arg)>

Converts a Perl logical value (1 or 0) to a PostgreSQL boolean (C<t> or C<f>).

=cut
    my ($arg) = @_;

    if ($arg) {
	return 't';
    }
    else {
	return 'f';
    }
}

sub quotify {

=item C<quotify($st)>

Converts a character string, which may contain quotes or backslashes,
into a PostgreSQL single quoted string by doubling up the single
quotes and the backslashes.

=cut

    my ($arg) = @_;

    if (not defined $arg) {
	return 'NULL';
    }
    $arg =~ s/\\/\\\\/g;
    $arg =~ s/'/''/g; # ' Make emacs happy.
    return "'" . $arg . "'";
}

sub find_smallest {

=item C<find_smallest(@args)>

Return the index into the arguments of the smallest numeric value. If
the list is empty, the result is undefined.

=cut

    my ($i, $smallest, $ret);

    $ret = undef;
    for ($i = 0; $i < scalar(@_); $i++) {
	if (not defined $ret) {
	    $ret = $i;
	    $smallest = $_[$i];
	}
	elsif ($_[$i] < $smallest) {
	    $ret = $i;
	    $smallest = $_[$i];
	}
    }
    return $ret;
}


sub commify {

=item C<commify($digit_string>

Put commas every three digits in the string.

Code taken from the Perl Cookbook.

=cut
    
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub date {

=item C<date>

Return a well formatted date like Mon Feb 07 14:30:06 EST 2011

=cut
    
    return strftime("%a %b %d %T %Z %Y", localtime);
}

sub clean_undef {

=item C<clean_undef(@args)>

Replace all undefined elements in the argument list with null strings,
and return the new list. The original list is unmodified.

=cut
    
    my (@ret, $arg);

    foreach $arg (@_) {
	push (@ret, defined $arg ? $arg : "");
    }
    return @ret;
}

sub mark_undef {

=item C<mark_undef(@args)>

Replace all undef in the argument list with "undef" strings
and return the new list. The original list is unmodified.

=cut
    
    my (@ret, $arg);

    foreach $arg (@_) {
	push (@ret, defined $arg ? $arg : "undef");
    }
    return @ret;
}

sub find_num_in_array {

=item C<find_num_in_array($target, $arrayp)>

Find the number in C<$target> in the array referenced by C<$arrayp>
and return the index where it is first found. If it is not found,
return C<-1>.

=cut
    
    my $target = shift;
    my $arrayp = shift;

    for (my $i = 0; $i < scalar(@{$arrayp}); $i++) {
	if ($target == $arrayp->[$i]) {
	    return $i;
	}
    }
    return -1;
}

sub find_string_in_array {

=item C<find_string_in_array($target, $arrayp)>

Find the string in C<$target> in the array referenced by C<$arrayp>
and return the index where it is first found. If it is not found,
return C<-1>.

=cut

    my $target = shift;
    my $arrayp = shift;

    for (my $i = 0; $i < scalar(@{$arrayp}); $i++) {
	if ($target eq $arrayp->[$i]) {
	    return $i;
	}
    }
    return -1;
}

=item C<dos_chomp($line)>

Remove the CR LF from the end of a line. Only removes characters
if they match the expected LF first and then CR from the end. The
C<$line> varible is modified.

=cut

sub dos_chomp {
    if (substr($_[0], -1) eq chr(10)) {
	$_[0] = substr($_[0], 0, length($_[0]) - 1);
	if (substr($_[0], -1) eq chr(13)) {
	    $_[0] = substr($_[0], 0, length($_[0]) - 1);
	}
    }
    return $_[0];
}

=item C<md5sum($file)>

Compute the MD5 checksum of a file. If the file does not exist, return undef.

=cut

sub md5sum {

    my $file = shift;
    local (*FILE);
    open (FILE, $file) ||
	return undef;
    binmode(FILE);
    my $md5 = Digest::MD5->new->addfile(*FILE)->hexdigest;
    return $md5;
}


=item C<parse_ms_config($file)>

Very simple parser of Microsoft configuration files (.INI) files.
Returns an array of sections where each section is just a pointer
to an array of lines, except for the section name, which is parsed
out of the square brackets. Lines are dos_chomped.

=cut

sub parse_ms_config {
    my $file = shift;
    local (*FILE, $_);
    
    my @ret = ();
    if (not open (FILE, $file)) {
	carp "Unable to open $file: $!\n";
	return @ret;
    }
    my @cur_section;
    my $state = 'init';
    while (<FILE>) {
	dos_chomp($_);
	if (m/^\s*$/) {
	    if (scalar(@cur_section) == 0) {
		next;
	    }
	}
	if (m/^\[([^\]]+)\]/) {
	    my $section_name = $1;
	    if (scalar(@cur_section) == 0) {
		push (@cur_section, $section_name);
	    }
	    else {
		push (@ret, [ @cur_section ]);
		@cur_section = ($section_name);
	    }
	}
	else {
	    push (@cur_section, $_);
	}
    }
    if (scalar(@cur_section)) {
	push (@ret, [ @cur_section ]);
    }
    return @ret;
}

sub display_variables {

=item C<display_variables($q)>

Print HTML to display form variables for an HTML query. The argument
is CGI object from the CGI Perl module of Lincoln Stein.

Taken from SEEBUGS.

=cut

    my ($q) = @_;
    my ($param, $val);
    
    print "<pre>\n";
    printf "query_string = %s\n", $q->query_string();
    foreach $param ($q->param()) {
	foreach $val ($q->param($param)) {
	    printf "%s=%s\n", $param, $val;
	}
    }
    print "</pre>\n";
}

    
sub display_environment {

=item C<display_environment($q)>

Print HTML to display the shell environment. The argument is CGI
object from the CGI Perl module of Lincoln Stein. It will also display
CGI form variables if the function C<display_variables> exists.

=cut

    my ($q) = @_;
    local ($_);
    my ($line);

    print "<pre>\n";
    print "Environment:<br>\n";
    foreach (sort keys %ENV) {
	printf "%s = \"%s\"\n", $_, $ENV{$_};
    }
    print "</pre>\n<p>\n";
    if (defined &display_variables) {
	&display_variables($q);
    }
    print "<p>Current day and time is ", &date, "</p>\n";
}
    
sub random_string {

=item C<random_string([$length])>

Make a random alphanumeric string of the given length. If no length is
specified, then a length of 6 is used.

=cut

    my $length = shift || 6;
    confess "Length ($length) must be greater than 0 for random_string.\n"
	if $length < 1;
    my $ret = " " x $length;
    my $list = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

    for (my $i = 0; $i < $length; $i++) {
	my $ir = int(rand(length($list)));
	substr($ret, $i, 1) = substr($list, $ir, 1);
    }
    return $ret;
}

    
sub range_overlap {

=item C<range_overlap($start1, $end1, $start2, $end2)>

Return true (1) if coordinate ranges overlap
and false (0) if not. The ranges are presumed to be closed.

=cut

    my ($start1, $end1, $start2, $end2) = @_;

    confess "Undefined parameters for range_overlap\n"
	if not (defined $start1 and defined $end1 and defined $start2 and defined $end2);
    confess "start1 ($start1) > end1 ($end1) in range_overlap\n"
	if $start1 > $end1;
    confess "start2 ($start2) > end2 ($end2) in range_overlap\n"
	if $start2 > $end2;
    if ($end2 < $start1 or
	$end1 < $start2) {
	return 0;
    }
    else {
	return 1;
    }
}

=back

=head1 SEE ALSO

Bio::Frescobi::Sequtil(3)

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
