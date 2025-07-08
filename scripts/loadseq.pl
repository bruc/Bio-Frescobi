#!perl
# Load a bunch of sequences into the sequence database.
# Usage: loadseq fasta-sequence-file [fasta-quality-file]

use warnings;
use strict;

use SHA;

use Getopt::Long;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::Dbseq;
use Bio::Frescobi::Frescobi;
use Bio::Frescobi::Config;

my $usage = <<EOF;
loadseq.pl -dbname=name -library=name fasta-sequence-file [fasta-quality-file]
           [-driver=string]
           [-libdesc=string]
           [-repository=name]
           [-strip=string]...
           [-poststrip=string]...
           [-name_type=string]
           [-[no]mergeids]
           [-merge_string=string]
           [-taxon=number|'annotation']
           [-prefix=string]
           [-replace]
           [-minlength=integer]
           [-uniquify]
           [-uniquify_separator=string]
           [-[no]bars]
           [-seqtype=string]
           [-[no]skipdupl]
           [-[no]concatid]
           [-seqid_head=string]
           [-[no]checkseq]
           [-[no]notefile]
           [-[no]usedisk]
EOF
    
my $library = "";
my $libdesc = "";
my $repository = "";
my $name_type = "";
my $mergeids = 0;
my $merge_string = ':';
my $taxon = undef;
my $taxon_from_annotation = 0;
my $prefix = "";
my $replace = 0;
my $minlength = 30;
my $uniquify = 0;
my $uniquify_separator = undef;
my $bars = 0;
my $seq_type = "nucleic";
my $skipdupl = 0;
my $seqid_head = "";
my @strip = ();
my @poststrip = ();
my $concatid = 0;
my $checkseq = 1;
my $notefile = 0;
my $usedisk = 1;
my $dbname = "";
my $driver = 'Pg';

&GetOptions("dbname=s" => \$dbname,
	    "library=s" => \$library,
	    "driver=s" => \$driver,
	    "libdesc=s" => \$libdesc,
	    "repository=s" => \$repository,
	    "strip=s" => \@strip,
	    "poststrip=s" => \@poststrip,
	    "name_type=s" => \$name_type,
	    "mergeids!" => \$mergeids,
	    "merge_string=s" => \$merge_string,
	    "taxon=s" => \$taxon,
	    "prefix=s" => \$prefix,
	    "replace!" => \$replace,
	    "minlength=i" => \$minlength,
	    "uniquify!" => \$uniquify,
	    "uniquify_separator=s" => \$uniquify_separator,
	    "bars!" => \$bars,
	    "seqtype=s" => \$seq_type,
	    "skipdupl!" => \$skipdupl,
	    "concatid!" => \$concatid,
	    "seqid_head=s" => \$seqid_head,
	    "checkseq!" => \$checkseq,
	    "notefile!" => \$notefile,
	    "usedisk!" => \$usedisk);

if (scalar(@ARGV) < 1 or scalar(@ARGV) > 2) {
    print STDERR $usage;
    die;
}

if ($dbname eq "") {
    die "The database name, (-dbname), must be specified.\n";
}

if (defined $taxon and $taxon eq 'annotation') {
    $taxon = undef;
    $taxon_from_annotation = 1;
}
if ($replace and $uniquify) {
    print STDERR "The replace and uniquify options are in conflict. Only one may be specified.\n";
    die;
}

if ($concatid and $uniquify) {
    print STDERR "The concatid and uniquify options are in conflict. Only one may be specified.\n";
    die;
}
if ($library eq "") {
    print STDERR "A library name must be specified\n";
    die;
}

if ($seqid_head eq "") {
    if ($seq_type eq "nucleic") {
	$seqid_head = "N";
    }
    elsif ($seq_type eq "protein") {
	$seqid_head = "P";
    }
    else {
	$seqid_head = "U";
    }
}
my $units;
if ($seq_type eq "nucleic") {
    $units = "bp";
}
elsif ($seq_type eq "protein") {
    $units = "AA";
}
else {
    $units = "letters";
}
set_checkseq($checkseq);

my ($seq_file, $qual_file) = @ARGV;

$qual_file = "" if not defined $qual_file;

my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => 0,
				  default_echo => "yes");
# my %names_found = ();
# my %seqid_for_name = ();
my %quals_for_name = ();
my %qual_annotations_for_name = ();
if ($qual_file) {
    open (QUAL, "<$qual_file") ||
	die "Unable to open quality score file: $!\n";
    my $qual_line = <QUAL>;
    while (defined $qual_line) {
	my ($name, $annotation, $seq) = read_1_fasta_sequence(\*QUAL, $qual_line, 1);
	last if not defined $name;
	if (not $bars) {
	    $name =~ s/\|/_/g;
	}
	&strip_name($name);
	if (exists($quals_for_name{$name})) {
	    print STDERR "Sequence identifier ($name) duplicated in $qual_file\n";
	}
	$quals_for_name{$name} = $seq;
	$qual_annotations_for_name{$name} = $annotation;
    }
    close QUAL;
}

my ($seqobj, $curtime); 
foreach my $mode (("scan", "load")) {
    next if $mode eq "scan" and $uniquify;
    next if $mode eq "scan";
#     if ($mode eq "scan") {
# 	$result = $pg->query(PGRES_TUPLES_OK,
# 			     "select name, seqid from raw_seqs " .
# 			     " where library = " . quotify($library));
# 	for ($i = 0; $i < $result->ntuples; $i++) {
# 	    ($name, $seqid) = $result->fetchrow;
# 	    $seqid_for_name{$name} = $seqid;
# 	}
#     }
    if ($mode eq "load") {
	$seqobj = Bio::Frescobi::Dbseq::new($pg);
	$seqobj->set_seqid_head($seqid_head .
				$pg->nextval('seqid_serial') .
				"_");

	$curtime = $pg->curtime();
	if ($concatid) {
	    $seqobj->attribute_for_name("concatid");
	}
	else {
	    $seqobj->attribute_for_name("name");
	}
	$seqobj->setup_table("raw_seqs");
	$seqobj->replace($replace);
	$seqobj->uniquify_names($uniquify);
	$seqobj->usedisk($usedisk);
	if (defined $uniquify_separator) {
	    $seqobj->uniquify_separator($uniquify_separator);
	}
	$seqobj->default_seq_type($seq_type);
	$seqobj->checkseq($checkseq);
	$seqobj->begin;
    }
    
    open(SEQ, "<$seq_file") || die "Unable to open $seq_file: $!\n";
    my $seq_line = <SEQ>;
    
    while (defined $seq_line) {
	my ($name, $annotation, $seq) = read_1_fasta_sequence(\*SEQ, $seq_line, 0);
	last if not defined $name;
	if (not $bars) {
	    $name =~ s/\|/_/g;
	}
	&strip_name($name);
	my $qual = "";
	if (exists($quals_for_name{$name})) {
	    $qual = $quals_for_name{$name};
	}
	next if $seq eq "";
	if (length($seq) < $minlength) {
	    printf STDERR "Sequence $name is too short (%d $units)\n", length($seq) if $mode eq "load";
	    next;
	}
	if ($prefix ne "") {
	    $name = $prefix . $name;
	}
	if (length($annotation) > 1000000 and $mode eq "load") {
	    print STDERR "Annotation for $name is too long. Truncating...\n";
	    $annotation = substr($annotation, 0, 1000000) . "... (truncated)";
	}
	if ($notefile) {
	    $annotation .= " Read from ${seq_file}";
	}
	if ($mergeids) {
	    my @pieces = ();
	    push(@pieces, $repository) if $repository;
	    push(@pieces, $library) if $library;
	    push(@pieces, $name_type) if $name_type;
	    push(@pieces, $name);
	    $name = join($merge_string, @pieces);
	}
	if ($mode eq "load") {
# 	    if (exists($names_found{$name})) {
# 		print STDERR "Sequence named $name was already present.\n";
# 	    }
# 	    else { 
		my $id = "";
		if ($concatid) {
		    $id = "$repository:$library:$name_type:$name";
		}
		if ($taxon_from_annotation) {
		    $taxon = undef;
		    if ($annotation =~ m/\/taxid=(\d+)\s/) {
			$taxon = $1;
		    }
		}
		$seqobj->add_data({ "name" => $name,
				    "seq" => $seq,
				    "qual" => $qual,
				    "annotation" => $annotation,
				    "taxon" => $taxon,
				    "repository" => $repository,
				    "library" => $library,
				    "name_type" => $name_type,
				    "filtered" => "f",
				    "create_time" => $curtime},
				    "concatid" => $id);
# 	    }
	}
	else {
# 	    if (exists($seqid_for_name{$name})) {
# 		($db_seq, $db_qual) = retrieve_sequence_data($pg, $seqid_for_name{$name});
# 		if (uc($db_seq) eq uc($seq) and $db_qual eq $qual) {
# 		    $names_found{$name} = 1;
# 		}
# 	    }
	}
    }
    close SEQ;
    close QUAL;
    if ($mode eq "load") {
	$seqobj->finish;
	$seqobj = undef;
    }
}
update_duplicates($pg, "raw_seqs") unless $skipdupl;
system_with_check("update_libraries.pl -dbname=\"$dbname\" -driver=$driver", "echo");
if ($libdesc ne "") {
    $pg->command("update libraries set description = " . &quotify($libdesc) .
		 " where library = ". &quotify($library),
		 1);
}
    
sub strip_name {
    # Remove any prefixes in the argument.
    my $name = $_[0];

    if (defined $name) {
    	foreach my $st (@strip) {
	    if ($name =~ s/^$st//) {
		last;
	    }
	}
    	foreach my $st (@poststrip) {
	    if ($name =~ s/$st$//) {
		last;
	    }
	}
    }
    $_[0] = $name;
}


=head1 NAME

loadseq.pl - Load sequences into the database.

=head1 SYNOPSIS

 loadseq.pl -dbname=name -library=name fasta-sequence-file [fasta-quality-file]
           [-libdesc=string]
           [-repository=name]
           [-strip=string]...
           [-poststrip=string]...
           [-name_type=string]
           [-[no]mergeids]
           [-merge_string=string]
           [-taxon=number|'annotation']
           [-prefix=string]
           [-driver=string]
           [-replace]
           [-minlength=integer]
           [-uniquify]
           [-uniquify_separator=string]
           [-[no]bars]
           [-seqtype=string]
           [-[no]skipdupl]
           [-[no]concatid]
           [-seqid_head=string]
           [-[no]checkseq
           [-[no]notefile]
           [-[no]usedisk]

=head1 DESCRIPTION

The C<loadseq.pl> script is the sequence loader. Sequences must be loaded into
the database before they can be annotated or otherwise manipulated. This scripts takes care of identifying identical sequences between the existing sequences and the ones being loaded to ensure that identical sequences get identical sequence identifiers (seqid's).
Much of the work for loading sequences is handled by Bio::Frescobi::Dbseq Perl module.

=head1 ARGUMENTS AND OPTIONS

=over 4

=item C<< -dbname=name >>

Specifies the database name into which the sequences will be loaded. There are
currently two database systems that Frescobi can use, PostgreSQL and
SQLite. In the case of PostgreSQL, you must specify the name of the
database in this option, and you must configure environment variables
so the script can login. In the case of SQLite, the C<-dbname> option
is used to specify the filename of the SQLite database.

This option is mandatory.

=item C<< -driver=[Pg|SQLite] >>

Specifies the database engine to be used for data retrieval and
storage. See the C<-dbname> option above.

=item C<< -library=name >>

Mandatory option specifying the library name for the sequences to be
loaded. Library names can be specified in multiple C<loadseq.pl>
runs. These serve to help organize the sequences in the database.

=item C<< fasta-sequence-file >>

Specifies a file of sequences to be loaded using the FASTA format.

=item C<< fasta-quality-file >>

Specifies a file of quality scores for the sequences to be
loaded. This file is optional, and was originally meant for use with
Sanger sequence data. These quality scores must be specified as
integers separated by white space. The identifiers of the sequences
must be matched to the identifiers used in the C<fasta-sequence-file>
above.

=item C<< -libdesc=string >>

Specifies a description for the library. If a library description is
already available, use of this option will replace it.

=item C<< -repository=name >>

Specifies the repository for the sequences. The repository is just
another piece of metadata that can be used to describe sequences. It
is optional.

=item C<< -strip=string ... >>

This option specifies a prefix to be removed from the beginning of
sequence identifiers in the Fasta file. Multiple prefixes can be
specified by using this option multiple times.

=item C<< -poststrip=string].. >>

This option specifies a suffix to be removed from the end of
sequence identifiers in the Fasta file. Multiple suffixes can be
specified by using this option multiple times.

=item C<< -name_type=string >>

If the C<-mergeids> option is used, then this option specifies the
third component in the merged ID.

=item C<< -[no]mergeids >>

If specified, C<loadseq.pl> will merge the C<repository>, C<library>,
C<name_type>, and Fasta identifier separated by the C<merge_string>
option into a larger identifier for storage in the database.

=item C<< -merge_string=string >>

When C<mergeids> is used, this option specifies the string to separate
the components when they are concatenated together. The default is a
colon (":").

=item C<< -taxon=number|'annotation' >>

Specifies the taxon to be saved with the sequences being loaded. This
can be done in two ways. If this option is specified as "annotation",
then C<loadseq.pl> will look at the sequence annotation for the
string, "taxid=", followed by a number and white space, and use that
number for the txaon. Otherwise, it will just use the value of this
option.

=item C<< -prefix=string >>

If specified, this option specifies a prefix to be concatenated to the
beginning of the Fasta sequence identifiers with no extra delimited.

=item C<< -replace >>

If specified, any sequence in the database which matches the name will
be replaced by the sequence in the sequence file. If not specified,
then C<loadseq.pl> will use the C<-uniquify> option to decide what to
do if a new sequence has the same name as an existing sequence. The
options, C<-replace> and C<-uniquify>, are mutually exclusive.

=item C<< -minlength=integer >>

Specifies the minimum length of sequences to be added to the
database. Any sequence shorter than this option will be ignored. The
default is 30.

=item C<< -uniquify >>

If this option is specified, then the sequence loader will generate
unique names when new sequences have the same names as previous
sequences. This is done by generating random strings that are checked
for uniqueness and added with the C<uniquify_separator> string
separating the original name and the random string.

The options, C<-replace> and C<-uniquify>, are mutually exclusive.

=item C<< -uniquify_separator=string >>

See the C<-uniquify> option above for the use of this option. The
default is an underscore character ("_").

=item C<< -[no]bars >>

When specified, any vertical bars ("|") in the sequence name are
replaced with underscores ("_"). By default this option is true, so if
you don't want the vertical bar replacement, you must specify
C<-nobars>.

=item C<< -seqtype=string >>

Specify the type of sequence. Only two choices are allowed: "nucleic"
or "protein". The default is "nucleic".

=item C<< -[no]skipdupl >>

If this option is off, then the loader will update a table of
duplicate raw sequences where there are multiple identical sequences
with different names. If this option is on, then no update of
duplicate sequences will be performed.

=item C<< -[no]concatid >>

If this option is on, then an additional field (C<concatid>) will be
loaded with the concatenation of the C<repository>, C<library>,
C<name_type>, and Fasta identifier separated by colons. This option is
incompatiable with the C<-uniquify> option above.

=item C<< -seqid_head=string >>

The construction of the sequence identifier uses three components, the
C<seqid_head>, a serial number generated by each run of C<loadseq.pl>
or the C<Bio::Frescobi::Dbseq> Perl module, and serial number for each
new sequence found in a sequence loading run. By default, the
C<seqid_head> is "N" for nucleic acid sequences and "P" for protein
sequences. The C<seqid_head> can be changed using the C<-seqid_head>
option.

=item C<< -[no]checkseq >>

If this option is turned on, the loader will check that all sequence
characters are either alphanumerica or the asterisks. Anything else
will be deleted. This option is turned on by default.

=item C<< -[no]notefile >>

If this option is turned on, then the following note, C< Read from
${seq_file}>, is added to the annotation for each sequence. The
C<${seq_file}> is the sequence file specified on the command line.

=item C<< -[no]usedisk >>

This is option must be left on which is the default.

=back

=head1 SEE ALSO

Bio::Frescobi::Dbseq(3pm)

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
