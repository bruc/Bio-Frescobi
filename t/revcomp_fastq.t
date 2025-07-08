# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl count_fastq.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 38);
use Test::More qw(no_plan);

if (not chdir('t')) {
    BAIL_OUT("Unable to cd to test directory, t\n");
}
$ENV{PATH} = "../blib/script:" . $ENV{PATH};
use lib "../blib/lib", "../blib/arch" ;

#########################

ok(0 == system("revcomp_fastq.pl sample_pe1.fastq revcomp_sample_pe1.fastq"),
   "fastq_recomp.pl");
my $test_file = "revcomp_sample_pe1.fastq";
ok(system("diff $test_file ref >${test_file}.dif") == 0, "Diff $test_file");

