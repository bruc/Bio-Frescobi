#!perl
# Update the library codes.

use warnings;
use strict;

use Getopt::Long;

use Bio::Frescobi::CgPg;
use Bio::Frescobi::PGLock;
use Bio::Frescobi::Frescobi;
use Bio::Frescobi::Config;
use Bio::Frescobi::Genutil;

my $usage = <<EOF;
update_libraries.pl -dbname=name [-driver=string]
EOF

my $driver = "Pg";
my $dbname = "";

&GetOptions("dbname=s" => \$dbname,
	    "driver=s" => \$driver);

if (scalar(@ARGV) != 0) {
    print STDERR $usage;
    die;
}

if ($dbname eq "") {
    die "The database name, (-dbname), must be specified.\n";
}

my %library_code;
my $max_library_code = 0;
my %library_description;

my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => "no");
$pg->default_echo(1);
&Bio::Frescobi::PGLock::grab_lock($pg, "libraries");
&get_libraries($pg);

my @data = @{$pg->query("select count(*),library from raw_seqs group by library")};
my %found = ();

foreach my $rowp (@data) {
    my $library = $rowp->[1];
    $found{$library} = 1;
}

foreach my $library (keys %library_code) {
    next if $library_code{$library} == 0;
    if (not exists($found{$library})) {
	$pg->command("delete from libraries where library = " . quotify($library));
    }
}

foreach my $rowp (@data) {
    my ($count, $library) = @{$rowp};
    if (not exists($library_code{$library})) {
	$max_library_code += 1;
	my $q_library = &quotify($library);
	$pg->command("insert into libraries (library, description, code, seqcount) " .
		     " values($q_library, $q_library, $max_library_code, $count)");
    }
    else {
	my $current_count = $pg->get_single_value("select seqcount from libraries " .
						  " where library = " . quotify($library));
	if ($count != $current_count) {
	    $pg->command("update libraries set seqcount = $count " .
			 " where library = " . quotify($library));
	}
    }
}

@data = @{$pg->query("select create_time, number, library from contig_run " .
		     " where completed " .
		     " order by create_time")};
my %run_number = ();
foreach my $rowp (@data) {
    my ($create_time, $number, $library) = @{$rowp};
    $run_number{$library} = $number;
}

foreach my $library (keys %run_number) {
    my $number = $run_number{$library};
    if ($pg->table_exists("contig_$number") and
	$pg->table_exists("singlets_$number")) {
	my $count = $pg->get_single_value("select count(*) from contig_$number") +
	    $pg->get_single_value("select count(*) from singlets_$number");
	my $current_count = $pg->get_single_value("select assembly_count from libraries " .
					       " where library = " . quotify($library));
	$current_count = -1 if not $current_count;
	if ($count != $current_count) {
	    $pg->command("update libraries " .
			 "   set assembly_count = $count " .
			 " where library = " . quotify($library));
	}
    }
    else {
	print STDERR "No contig or singlet found for library $library, run_number = $number\n";
    }
}
&Bio::Frescobi::PGLock::free_lock($pg, "libraries");

sub get_libraries {
    my $pg = shift;
    # Get the available libraries from Frescobi into the global variables,
    # %library_description, %library_code, and $max_library_code.
    
    my ($i, $library, $description, $code);

    $max_library_code = 0;
    foreach my $rowp (@{$pg->query("select library, description, code from libraries", 1)}) {
	($library, $description, $code) = @{$rowp};
	$library_description{$library} = $description;
	$library_code{$library} = $code;
	$max_library_code = &max($max_library_code, $code);
    }
}
