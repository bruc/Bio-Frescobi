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

ok(0 == system("fastq_sample.pl -seed=20 -sample=20 -output=random_sample sample_pe1.fastq sample_pe2.fastq"),
   "Sample reads");
my @test_files = qw(random_sample.1.fq random_sample.2.fq);

foreach my $test_file (@test_files) {
    ok(system("diff $test_file ref >${test_file}.dif") == 0, "Diff $test_file");
}

