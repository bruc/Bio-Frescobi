# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pematch.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 38);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::Sequtil');
    use_ok('Bio::Frescobi::Genutil');
    use_ok('Bio::Frescobi::PEMatch');
};

exit;
if (not chdir('t')) {
    BAIL_OUT("Unable to cd to test directory, t\n");
}
$ENV{PATH} = "../blib/script:" . $ENV{PATH};
use lib "../blib/lib", "../blib/arch" ;

#########################

ok(my $pematch = new Bio::Frescobi::PEMatch, "object creation");
my $dna = "CGACAGCTTGAAATAGATGTTAATCTTTGTTATTGTCGTATGTCATGTATCACCACATTGATCACGTGTTCAAAGACAGGATGTTGGCCGGGTGTCGTGG";
my $revdna = reverse_dna_seq($dna);
my $qual = "a" x 70;
my ($joinseq, $joinqual) = $pematch->match(substr($dna, 0, 70),
					   $qual,
					   substr($revdna, 0, 70),
					   $qual);
ok($joinseq eq $dna, "Single test");

system_with_echo("pematch.pl -minoffset=10 " .
		 "           -maxoffset=90 " .
		 "           -bad=sample_pe.d0.bad " .
		 "           -maxd=0 " .
		 "           -debug=2 " .
		 "           sample_pe1.fastq sample_pe2.fastq " .
		 "           >sample_pe.d0.out 2>sample_pe.d0.err", 1);
system_with_echo("pematch.pl -minoffset=10 " .
		 "           -maxoffset=90 " .
		 "           -bad=sample_pe.d5.bad " .
		 "           -maxd=5 " .
		 "           -debug=2 " .
		 "           sample_pe1.fastq sample_pe2.fastq " .
		 "           >sample_pe.d5.out 2>sample_pe.d5.err", 1);
system_with_echo("pematch.pl -minoffset=10 " .
		 "           -minoverlap=10 " .
		 "           -maxoffset=-1 " .
		 "           -bad=sample_pe.cond3.bad " .
		 "           -maxd=1 " .
		 "           -debug=2 " .
		 "           -matchweight=1 " .
		 "           sample_pe1.fastq sample_pe2.fastq " .
		 "           >sample_pe.cond3.out 2>sample_pe.cond3.err", 1);

my @test_files;
foreach my $tag (('d0', 'd5', 'cond3')) {
    push (@test_files,
	  "sample_pe.${tag}.bad",
	  "sample_pe.${tag}.out",
	  "sample_pe.${tag}.err");
}

foreach my $test_file (@test_files) {
    ok(system("diff $test_file ref >${test_file}.dif") == 0, "Diff $test_file");
}

