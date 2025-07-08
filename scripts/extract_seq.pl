#!perl

# Simple script to extract all sequences from a database.

use strict;
use warnings;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Config;
use Getopt::Long;
use Bio::Frescobi::BigQuery;
use Bio::Frescobi::PGLock;
use Bio::Frescobi::PGLoad;

my $usage = <<EOF;
extract_seq.pl -dbname=dbname         [output-file]
             [ -driver=driver ]  Default: Pg
             [ -sql=sql-statement ] 
             [ -type=protein|nucleic ]
             [ -tablepat=tablepat ]
             [ -[no]quality ]
             [ -[no]sequences ]
             [ -id=sql-statement ]
             [ -annotlim=integer ]
             [ -prefix=string ]
             [ -header=sql-statement ... ]
             [ -headers=sql-statement ]
             [ -linesize=integer ]
             [ -[no]lock ]        Default: -lock
             [ -[no]cleanseq ]    Default: -cleanseq
             [ -[no]cleanheader ] Default: -cleanheader
             [ -[no]verbose ]     Default: -noverbose
             [ -help
	     
    Default annotlim (annotation length limit) = 500
EOF

my $dbname = "";
my $driver = "Pg";
my $select_sql = "";
my $type = "";
my $tablepat = "";
my $do_sequences = 1;
my $do_quality = 0;
my @header_sql = ();
my $headers_sql = "";
my $id_sql = "";
my $annotlim = 500;
my $prefix = "";
my $help = 0;
my $lock = 0;
my $linesize = undef;
my $cleanseq = 1;
my $cleanheader = 1;
my $verbose = 0;
GetOptions("dbname=s" => \$dbname,
	   "driver=s" => \$driver,
	   "sql=s" => \$select_sql,
	   "type=s" => \$type,
	   "tablepat=s" => \$tablepat,
	   "quality!" => \$do_quality,
	   "sequences!" => \$do_sequences,
	   "annotlim=i" => \$annotlim,
	   "prefix=s" => \$prefix,
	   "header=s" => \@header_sql,
	   "headers=s" => \$headers_sql,
	   "id=s" => \$id_sql,
	   "linesize=i" => \$linesize,
	   "lock!" => \$lock,
	   "cleanseq!" => \$cleanseq,
	   "cleanheader!" => \$cleanheader,
	   "verbose!" => \$verbose,
	   "help!" => \$help);

if ($help) {
    print STDERR $usage;
    exit 0;
}

if ($dbname eq "") {
    die "The database name, (-dbname), must be specified.\n";
}

if ($type && $type !~ m/^(protein|nucleic)$/) {
    die "Sequence type must be one of \"protein\" or \"nucleic\"\n";
}

if ($do_quality and $do_sequences) {
    die "Only quality or sequences can be specified.\n";
}

if (not ($do_quality or $do_sequences)) {
    die "One of quality or sequences must be specified.\n";
}

if (defined $linesize) {
    set_format_seq_size($linesize);
}

if (scalar(@header_sql) > 0 and $headers_sql) {
    die "Only one of -header or -headers can be specified.\n";
}

if (scalar(@ARGV) == 1) {
    open(STDOUT, ">$ARGV[0]") || die "Unable to open $ARGV[0]: $!\n";
}
elsif (scalar(@ARGV) > 1) {
    print STDERR $usage;
    exit(1);
}

my $miscreants = '/\|()[{^$*+?.-';
my $header_allowed = "";
for (my $i = 32; $i <= 126; $i++) {
    my $ch = chr($i);
    if (index($miscreants, $ch) > -1) {
	$header_allowed .= "\\";
    }
    $header_allowed .= chr($i);
}
eval 'sub clean_annotation { my $annotation = shift; $annotation =~ tr/' . $header_allowed . '//dc; return $annotation; }';
die "Unable to generate clean_annotation subroutine: $@\n" if $@;
    
my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => "no");
if ($verbose) {
    $pg->default_echo("on");
}
my $all_seqs = 1;
$all_seqs = 0 if $select_sql ne "" or $tablepat ne "";

&Bio::Frescobi::PGLock::grab_lock($pg, "seq_data") if $lock;
END {
    &Bio::Frescobi::PGLock::free_lock($pg, "seq_data") if $lock;
}

my $temp_table1 = "extract_seq1_$$";
my $temp_table2 = "extract_seq2_$$";
my $temp_table3 = "extract_seq3_$$";
my $sql;
my $piece;
$piece = "seq_piece" if $do_sequences;
$piece = "qual_piece" if $do_quality;

$pg->command("create temporary table $temp_table1 (seqid text)", 1);
if ($all_seqs) {
    $pg->command("insert into $temp_table1 " .
		 "select d.seqid " .
		 "  from seq_data d ");
}
else {
    if ($select_sql) {
	$pg->command("insert into $temp_table1 $select_sql", 1);
    }
    if ($tablepat) {
	my @tables = grep {m/$tablepat/} $pg->get_all_tables;
	foreach my $table (@tables) {
	    $pg->command("insert into $temp_table1 " .
			 "select seqid from $table",
			 1);
	}
    }
}
$pg->command("create temporary table $temp_table2 " .
	     " as " .
	     "select distinct t.seqid " .
	     "  from $temp_table1 t, seq_data d " .
	     " where t.seqid = d.seqid " .
	     ($type ? "and d.seq_type = '$type'" : ""),
	     1);
if ($headers_sql) {
    $pg->command("create temporary table $temp_table3 ( " .
		 "       seqid      text, " .
		 "       annot      text)");
    my $loader = new Bio::Frescobi::PGLoad($pg, $temp_table3);
    my $query = new Bio::Frescobi::BigQuery($pg, $headers_sql);
    while (my $row = $query->next) {
	my @data = @{$row};
	my $seqid = shift(@data);
	my $annot = join(" ", @data);
	$loader->add_data({seqid => $seqid,
			   annot => $annot});
    }
    $loader->finish;
    printf STDERR "A total of %d annotations written to $temp_table3\n", $loader->record_count;
    $sql = qq( select t.seqid,
                      p.segnum,
                      p.$piece,
                      string_agg(h.annot, ' ') as annot
	         from seq_pieces p,
                      $temp_table2 t,
                      $temp_table3 h
	        where t.seqid = p.seqid
                  and t.seqid = h.seqid
                group by t.seqid, p.segnum, p.$piece
                order by t.seqid, p.segnum );
}
else {
    $sql = qq( select t.seqid, p.segnum, p.$piece, '' as annot
	     from seq_pieces p, $temp_table2 t
	    where t.seqid = p.seqid 
            order by p.seqid, p.segnum );
}
$pg->default_echo(1);
my $query = new Bio::Frescobi::BigQuery($pg, $sql, 10000);
$pg->default_echo($verbose);

my $row_count = 0;
my $row = $query->next;
while (defined $row) {
    $row_count += 1;
    if (not $verbose) {
	if ($row_count % 10000 == 0) {
	    print STDERR "At ", &date, " $row_count records retrieved.\n";
	}
    }
    my ($seqid, $segnum, $piece, $annot) = @{$row};
    my @pieces = ();
    push (@pieces, $piece);
    while (defined ($row = $query->next)) {
	last if $seqid ne $row->[0];
	if ($row->[1] != $segnum + 1) {
	    print STDERR "For seqid $seqid, segnum ($row->[1] $segnum) out of order\n";
	}
	$segnum = $row->[1];
	$piece = $row->[2];
	push(@pieces, $piece) if $piece;
    }
    my $data = join("", @pieces);
    if ($cleanseq) {
	if ($do_quality) {
	    if ($data =~ m/[^ 0-9]/) {
		print STDERR "For seqid $seqid, quality data has junk characters\n";
		$data =~ tr/ 0-9//dc;
	    }
	}
	elsif ($do_sequences) {
	    $data =~ s/\?/X/g;
	    if ($data =~ m/[^\*a-zA-Z]/) {
		print STDERR "For seqid $seqid, sequence data has junk characters\n";
		$data =~ tr/[^\*a-zA-Z]//dc;
	    }
	}
    }
    my @annotations = ();
    if (scalar(@header_sql) > 0) {
	foreach my $header_sql (@header_sql) {
	    my $sql = $header_sql;
	    $sql =~ s/\$seqid/$seqid/g;
	    my $hrows = $pg->get_all_rows($sql);
	    my $l = 0;
	    foreach my $hrow (@{$hrows}) {
		my $piece = " " . join(" ", clean_undef(@{$hrow}));
		$l += length($piece);
		push(@annotations, $piece);
		last if $l > $annotlim + 1;
	    }
	    last if $l > $annotlim + 1;
	}
    }
    else {
	push(@annotations, " ", $annot);
    }
    my $annotation = join("", @annotations);
    my $ctrl_a = chr(1);
    $annotation =~ s/$ctrl_a/ /g;
    if ($cleanheader) {
	my $old_length = length($annotation);
	$annotation = &clean_annotation($annotation);
	if (length($annotation) < $old_length) {
	    print STDERR "Non-ASCII characters found in annotation for seqid $seqid\n";
	}
    }
    if (length($annotation) > $annotlim + 1) {
	$annotation = substr($annotation, 0, $annotlim + 1);
    }
    # It appears that a single blank annotation causes serious problems for critica,
    # so delete such annotations.
    if ($annotation =~ m/^\s+$/) {
	$annotation = "";
    }
    my $id = $seqid;
    if ($id_sql ne "") {
	my $tsql = $id_sql;
	$tsql =~ s/\$seqid/$seqid/g;
	$id = $pg->get_single_value($tsql);
	$id = $seqid if not $id;
    }
    print ">${prefix}${id}$annotation\n";
    print format_seq($data) if $do_sequences;
    print format_qual($data) if $do_quality;
}

print STDERR "At ", &date, " a total of $row_count records retrieved.\n";

=head1 NAME

extract_seq.pl - Extract sequences in Fasta format

=head1 SYNOPSIS

 extract_seq.pl -dbname=dbname         [output-file]
              [ -driver=driver ]  Default: Pg
              [ -sql=sql-statement ] 
              [ -type=protein|nucleic ]
              [ -tablepat=tablepat ]
              [ -[no]quality ]
              [ -[no]sequences ]
              [ -id=sql-statement ]
              [ -annotlim=integer ]
              [ -prefix=string ]
              [ -header=sql-statement ... ]
              [ -headers=sql-statement ]
              [ -linesize=integer ]
              [ -[no]lock ]        Default: -lock
              [ -[no]cleanseq ]    Default: -cleanseq
              [ -[no]cleanheader ] Default: -cleanheader
              [ -[no]verbose ]     Default: -noverbose
              [ -help

=head1 DESCRIPTION

The C<extract_seq.pl> command extracts sequences into a Fasta file
with great flexibility. It uses SQL statements to retrieve annotation
text.

Because Frescobi associates sequence identifiers (C<seqid's>) with
sequences, any redundant sequences can have multiple annotations and
names. This multiplicity can complicate the generation of names and
annotations in the sequence output.

=head1 ARGUMENTS AND OPTIONS

=over 4

=item C<< -dbname=dbname >>

Specifies the database name from where the sequences will be
extracted.  There are currently two database systems that Frescobi can
use, PostgreSQL and SQLite. In the case of PostgreSQL, you must
specify the name of the database in this option, and you must
configure environment variables so the script can login. In the case
of SQLite, the C<-dbname> option is used to specify the filename of
the SQLite database.

This option is mandatory.

=item C<< output-file >>

If an C<output-file> is optionally specified, then the sequences will
be written to this file. If no C<output-file> is specified, then the
output goes to standard output.

=item C<< -driver=[Pg|SQLite]  >>

Specifies the database engine to be used for data retrieval and
storage. See the C<-dbname> option above.

=item C<< -sql=sql-statement ] >>

Specifies an SQL statement that will select the sequences to be
output. The SQL statement must have a C<seqid> field which is used to
retrieve the sequences.  If this option and the C<-tablepat> option
are not used, then all sequences in the database will be output.

=item C<< -type=protein|nucleic  >>

Specifies the type of sequence to be extracted. If no C<-type> is
specified, then all sequences will be output.

=item C<< -tablepat=tablepat  >>

Specifies a Perl pattern to select tables in the database that have
C<seqid's> for extraction. If this option and the C<-sql> option
are not used, then all sequences in the database will be output.

=item C<< -[no]quality  >>

If specified, the quality scores are output. Only one of quality
scores or sequences can be output in one run.

=item C<< -[no]sequences  >>

If specified, the sequences are output. This is the default. Only one
of quality scores or sequences can be output in one run.

=item C<< -id=sql-statement  >>

This option specifies an SQL statement that will determine what
identifier will be output for each sequence. This SQL statement is
executed for every sequence being output, and it must contain the
string, C<$seqid>, in the statement. When a sequence is output, the
string, C<$seqid>, gets replaced with the actual sequence identifier
being output. The SQL statement must return a single value for that
sequence identifier.

This option can be tricky to debug so use the C<-verbose> option to
see the actual SQL statement used for each sequence identifier.

=item C<< -annotlim=integer  >>

Specifies a length limit for the annotation of each sequence. The
default is 500 characters.

=item C<< -prefix=string  >>

If specified, this string is added to the beginning of each sequence identifier.

=item C<< -header=sql-statement ...  >>

This option specifies multiple SQL statements that will be used to retrieve annotation information
for each sequence. This option is mutually exclusive with the C<-headers> option below.

These SQL statements are executed for every sequence being output, and
it must contain the string, C<$seqid>, in the statement. When a
sequence is output, the string, C<$seqid>, gets replaced with the
actual sequence identifier being output. These SQL statements can
return multiple fields per sequence identifier, and all the results
are concatenated together. Undefined fields are ignored. Many special
characters are also removed and Ctrl-A characters are changed into
spaces. Ctrl-A characters are typically found in NCBI sequence
database annotations. The C<-annotlim> also applies to limit the
length of the annotation being output.

=item C<< -headers=sql-statement  >>

This option specifies one SQL statement that will be executed once to
accumulate all the annotation information in a single query to the
database, in contrast to the SQL query per sequence method used by the
C<-header> option above. This option is mutually exclusive with
C<-header> option.

The SQL statement for this option must return only two fields, the
C<seqid> first and the annotation string second. The concatenation
operator or function in the database engine can be very useful in
generating an informative annotation for each sequence.

=item C<< -linesize=integer  >>

Specifies the line size for the sequence data in the output file. The
default is specified in the Bio::Frescobi::Sequtil module.

=item C<< -[no]lock ]     >>

If specified, the sequence data is locked while the output is being
generated. This option is useful if there could be parallel processes
loading sequence data at the same time. The default is not to lock the
sequence data.

=item C<< -[no]cleanseq ] >>

If specified, the sequence data is cleaned of any non alphanumeric
characters or an asterisk ("*"). The default is to clean the
sequences.

=item C<< -[no]cleanheader  >>

If specified, the header is cleaned of any unprintable or otherwise undesirable characters.
The default is to clean the header.

=item C<< -[no]verbose ]    >>

If specified, turn on more logging of SQL statements.

=item C<< -help >>

Prints the short usage text and exits.

=back

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
