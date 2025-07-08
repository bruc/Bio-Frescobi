# Run the regression test script, regress.sh

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 38);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::Genutil');
};

if (not chdir('t')) {
    BAIL_OUT("Unable to cd to test directory, t\n");
}
$ENV{PATH} = "../blib/script:" . $ENV{PATH};
use lib "../blib/lib", "../blib/arch" ;
#########################

foreach my $driver (('Pg', "SQLite")) {
    diag("About to run the regression test using the $driver driver. This will take a while.\n");
    my $status = system("./regress.sh $driver > regress.$driver.log 2>&1");
    ok ($status == 0, "./regress.sh $driver");
    if ($status == 0) {
	ok(system("diff regress.$driver.log ref/regress.$driver.ref.log") == 0,
	   "Diff regress.$driver.log");
    }
}
