# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl PGLoad.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 56);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::PGLoad');
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
$pg->command("drop table if exists load_test");
$pg->command("create table load_test (tag text, value text, prob double precision)");
$pg->command("create index load_test_tag on load_test using btree (tag)");

my $loader = Bio::Frescobi::PGLoad->new($pg, "load_test");
ok($loader->drop_indexes(1) == 1, "drop index setting");
my @data = ( [ "height", "6 feet", 0.5 ],
	     [ "weight", "190 lbs", 0.7 ],
	     [ "age", "45", 0.9 ],
	     [ "sex", "occasionally", 0.5 ] );
foreach my $datum (@data) {
    $loader->add_data({ "tag" => $datum->[0],
			"prob" => $datum->[2],
			"value" => $datum->[1] });
}
$loader->finish;
ok($loader->record_count == scalar(@data),
   "record count");
foreach my $datum (@data) {
    my ($tag, $value, $prob) = @{$datum};
    my ($ret_value, $ret_prob) = $pg->get_single_row("select value, prob from load_test where tag = " . quotify($tag));
    ok($ret_value eq $value,
       "check value");
    ok($ret_prob eq $prob,
       "check probability");
}
__END__
