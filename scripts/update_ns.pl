#!perl

# Simple script to update the table of N's. 

use warnings;
use strict;

use Bio::Frescobi::Genutil;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::CgPg;
use Getopt::Long;
use Bio::Frescobi::PGLoad;

my $dbname;
my $driver = "Pg";
my $help = 0;
my $verbose = 0;
my $replace = 1;
my $collection = "";

my $usage = <<EOF;
update_ns.pl [ -dbname=dbname       ] <fasta file>
             [ -driver=driver       ]
             [ -[no]replace         ]
             [ -[no]verbose         ] Default: -noverbose
             [ -[no]help            ]
             [ -collection=<string> ] 
EOF
    ;

GetOptions("dbname=s" => \$dbname,
	   "driver=s" => \$driver,
	   "collection=s" => \$collection,
	   "replace!" => \$replace,
	   "verbose!" => \$verbose,
	   "help!" => \$help);

if ($help) {
    print STDERR $usage;
    exit 0;
}

if (scalar(@ARGV) == 1) {
    open(SEQ, "<$ARGV[0]") || die "Unable to open $ARGV[0]: $!\n";
}
else {
    print STDERR $usage;
    exit(1);
}

my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => "no");
if ($verbose) {
    $pg->default_echo("on");
}

if (not $pg->table_exists("seq_ns")) {
    $pg->command("create table seq_ns (" .
		 "       collection text, ". 
		 "       name       text, " .
		 "       start      integer, " .
		 "       stop       integer)");
    if ($pg->driver() eq 'Pg') {
	$pg->command("create index seq_ns_name on seq_ns using btree(name, collection)");
    }
    else {
	$pg->command("create index seq_ns_name on seq_ns (name, collection)");
    }
}

my $loader = new Bio::Frescobi::PGLoad($pg, "seq_ns");

my $line = <SEQ>;
while (defined $line) {
    my ($name, $annotation, $seq) =
	read_1_fasta_sequence(\*SEQ, $line, 0);
    last if not defined $name;
    if ($replace) {
	$pg->command("delete from seq_ns " .
		     " where name = " . quotify($name) .
		     "   and collection = " . quotify($collection));
    }
    elsif ($pg->get_single_value("select count(*) from seq_ns where name = " . quotify($name)) != 0) {
	print STDERR "$name already has data.\n" if $verbose;
	next;
    }
    foreach my $pair (find_ns($seq)) {
	$loader->add_data({collection => $collection,
			   name => $name,
			   start => $pair->[0],
			   stop => $pair->[1]});
    }
}
$loader->finish;
printf STDERR "%d N blocks added to $dbname.\n", $loader->record_count;
