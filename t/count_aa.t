# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl count_fastq.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More (tests => 2);
# use Test::More qw(no_plan);

if (not chdir('t')) {
    BAIL_OUT("Unable to cd to test directory, t\n");
}
$ENV{PATH} = "../blib/script:" . $ENV{PATH};
use lib "../blib/lib", "../blib/arch" ;

#########################

ok(0 == system("count_aa.pl <exonerate_prot_db > count_aa.out"),
   "Run count_aa.pl");
ok(system("diff count_aa.out ref >count_aa.dif") == 0, "Diff count_aa.out");

