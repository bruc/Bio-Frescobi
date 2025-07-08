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

if (not -d "tmp") {
    if (system("mkdir tmp") != 0) {
	BAIL_OUT("Unable to create tmp directory.\n");
    }
}
ok(0 == system("count_fastq.pl -tmpdir=tmp " .
	       "               -mincount=3 " .
	       "               -usememory " .
	       "               -trim=5 " .
	       "   <sample_pe1.fastq >tmp/sample_pe1.count.1"),
  "Use memory run");
ok(0 == system("count_fastq.pl -tmpdir=tmp " .
	       "               -mincount=3 " .
	       "               -nousememory " .
	       "               -trim=5 " .
	       "               -sortmem=10m " .
	       "   <sample_pe1.fastq >tmp/sample_pe1.count.2"),
  "Unix sort run");
my @test_files = qw(sample_pe1.count.1 sample_pe1.count.2);

foreach my $test_file (@test_files) {
    ok(system("diff tmp/$test_file ref >${test_file}.dif") == 0, "Diff $test_file");
}

