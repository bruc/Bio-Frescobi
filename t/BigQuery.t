# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl BigQuery.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 56);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::BigQuery');
};
use Bio::Frescobi::Genutil;
use Data::Dumper;

#########################

my $dbname = $ENV{"PG_TEST_DB"} || "regression";
my $pg = new Bio::Frescobi::CgPg(dbname => $dbname,
				 cgi => 0,
				 default_echo => 1);
ok($pg, "open connection");
BAIL_OUT("Unable to open connection to $dbname. Aborting test.\n") if not $pg;
ok($pg->command("drop table if exists bigquery_test_table"), "Drop table");
ok($pg->command("create table bigquery_test_table (" .
		"   i     integer, " .
		"   f     float)"),
   "Create table");
my $test_size = 2000;
$pg->autocommit(0);
for (my $i = 1; $i <= $test_size; $i++) {
    $pg->command("insert into bigquery_test_table (i, f) values ( " .
		 join(", ", quotify($i), quotify($i*$i + 0.5)) . ")",
		0);
}
$pg->autocommit(1);
ok($pg->get_single_value("select count(*) from bigquery_test_table") == $test_size,
   "Fill test table");

my $query = new Bio::Frescobi::BigQuery($pg,
					"select i, f from bigquery_test_table",
					$test_size / 4);
ok($query->verbose_mode(1),
   "Turn on verbose mode");

my $err = 0;
while (my $rowp = $query->next) {
    my ($i, $fuzzy_square) = @{$rowp};
    if ($fuzzy_square - 0.5 != $i * $i) {
	print STDERR "Mismatch found for $i and $fuzzy_square\n";
	$err = 1;
    }
}

ok($err == 0, 
   "Readback of table is wrong");

__END__
