# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Sequtil.t'

use warnings;
use strict;
use Bio::Frescobi::Genutil;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More (tests => 53);
# use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::Sequtil');
};
use Data::Dumper;

#########################

my $test_seq =
    "GATCCTCCCCAGGCCCCTACACCCAATGTGGAACCGGGGTCCCGAATGAA" .
    "AATGCTGCTGTTCCCTGGAGGTGTTTTCCTGGACGCTCTGCTTTGTTACC" .
    "AATGAGAAGGGCGCTGAATCCTCGAAAATCCTGACCCTTTTAATTCATGC" .
    "TCCCTTACTCACGAGAGATGATGATCGTTGATATTTCCCTGGACTGTGTG" .
    "GGGTCTCAGAGACCACTATGGGGCACTCTCGTCAGGCTTCCGCGACCACG" .
    "TTCCCTCATGTTTCCCTATTAACGAAGGGTGATGATAGTGCTAAGACGGT" .
    "CCCTGTACGGTGTTGTTTCTGACAGACGTGTTTTGGGCCTTTTCGTTCCA" .
    "TTGCCGCCAGCAGTTTTGACAGGATTTCCCCAGGGAGCAAACTTTTCGAT" .
    "GGAAACGGGTTTTGGCCGAATTGTCTTTCTCAGTGCTGTGTTCGTCGTGT" .
    "TTCACTCACGGTACCAAAACACCTTGATTATTGTTCCACCCTCCATAAGG" .
    "CCGTCGTGACTTCAAGGGCTTTCCCCTCAAACTTTGTTTCTTGGTTCTAC" .
    "GGGCTG";

my $test_qual =
    "26 25 26 27 27 27 24 27 27 22 25 28 20 24 27 32 30 22 12 2 " .
    "22 27 27 23 26 27 27 27 27 26 25 27 22 32 30 17 23 27 27 24 " .
    "26 27 25 16 24 26 27 19 31 26 17 25 27 25 32 29 18 2 30 25 " .
    "26 27 23 24 26 30 24 26 27 25 27 27 26 26 24 27 6 23 27 25 " .
    "24 27 25 31 26 26 26 32 30 19 4 26 27 27 31 26 24 31 25 30 " .
    "24";

set_checksum_size(16);
ok(get_checksum_size() == 16, "get_checksum_size");
ok(seq_checksum($test_seq, 8) eq '89fc701b', "seq_checksum 1");
ok(seq_checksum($test_seq) eq '89fc701badb98186', "seq_checksum 2");
ok(substr(format_seq($test_seq, 20), 0, 20) eq 'GATCCTCCCCAGGCCCCTAC', "format_seq 1");
ok(substr(format_seq($test_seq, 20), 42, 20) eq 'CCCGAATGAAAATGCTGCTG', "format_seq 2");
set_format_seq_size(100);
ok(get_format_seq_size() == 100, "get_format_seq_size");
ok(substr(format_seq($test_seq), 80, 20) eq 'GGACGCTCTGCTTTGTTACC', "format_seq 3");
set_checkseq(1);
ok(get_checkseq() == 1, "get_checkseq");
set_format_qual_size(75);
ok(get_format_qual_size() == 75, "get_format_qual_size");
ok((split(/\n/, format_qual($test_qual)))[2] eq '27 25 32 29 18 2 30 25 26 27 23 24 26 30 24 26 27 25 27 27 26 26 24 27 6 23', "format_qual");
my @quals = qual_to_array($test_qual);
ok($quals[14] eq '27', "qual_to_array 1");
ok(scalar(@quals) == 101, "qual_to_array 2");
ok(array_to_qual(@quals) eq $test_qual, "array_to_qual");
{
    local (*SEQ, *QUAL);
    
    open(SEQ, "<t/dr.seq.1.fasta") ||
	fail("Unable to open t/dr.seq.1.fasta: $!");
    open(QUAL, "<t/dr.seq.1.fasta.qual") ||
	fail("Unable to open t/dr.seq.1.fasta.qual: $!");
    my $line = <SEQ>;
    my ($name, $annotation, $sq) = read_1_fasta_sequence(\*SEQ, $line, 0);
    ok($name eq 'read_1/1', "read_1_fasta_sequence 1");
    ok($annotation eq 'annotation for sequence-1', "read_1_fasta_sequence 2");
    ok(seq_checksum($sq) eq '098c271260385eda', "read_1_fasta_sequence 3");
    my $qline = <QUAL>;
    my ($qname, $qannot, $q) = read_1_fasta_sequence(\*QUAL, $qline, 1);
    ok($qname eq 'read_1/1', "read_1_fasta_sequence 4");
    ok($qannot eq 'annotation for quality-1', "read_1_fasta_sequence 5");
    ok(seq_checksum($q) eq '27f3dae107ee7988', "read_1_fasta_sequence 6");
}
{
    my ($annotation, $seq) = extract_1_sequence("t/dr.seq.4.fasta", "read_8/1");
    ok(seq_checksum($seq) eq '945b474b28569ca7', "extract_1_sequence");
}
{
    local (*SEQ);
    open(SEQ, "<t/sample.fastq") ||
	fail("Unable to open t/sample.fastq: $!");
    my ($name, $annotation, $seq, $qual) = read_1_fastq_sequence(\*SEQ, 33);
    ($name, $annotation, $seq, $qual) = read_1_fastq_sequence(\*SEQ, 33);
    my @qual = fastq_qual_to_array($qual, chr(33));
    ok ($qual[6] == 27, "fastq_qual_to_array 1");
    ok(fastq_array_to_qual(\@qual, chr(33)) eq $qual, "fastq_qual_to_array 2");
    ok($name eq 'SRR022650.2', "read_1_fastq_sequence 1");
    ok($annotation eq "E3MFGYR01D8GNX length=95", "read_1_fastq_sequence 2");
    ok(seq_checksum($seq) eq "6ddd4205cf146ced", "read_1_fastq_sequence 3");
    ok(seq_checksum($qual) eq '95c9500cddcdfa8e', "read_1_fastq_sequence 4");
    ok(format_fastq($name, "new annotation", $seq . "N", $qual . "B") eq <<'EOF', "format_fastq");
@SRR022650.2 new annotation
AAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCGGGAATGCTAGCAAAGATCTGAAGGTGAAGAGGATCACACCAAGGCACACAGGGGATAGGN
+
B;====<B;<=<=C<<<=<<=<<C<C;<B;;=A9<EA2A9:<=<;=<C>*=<:;<<?:B<8=C=8:B<<9;<;:C=B<B<<=<2;3EA4%;9;B;B
EOF
}
ok(reverse_dna_seq("aattccgg") eq "ccggaatt", 'reverse_dna_seq');
ok(substr_dna_seq("GAGTGGCCATTACGGCG", 2, 5, 1) eq "AGTG", "substr_dna_seq 1");
ok(substr_dna_seq("GAGTGGCCATTACGGCG", 7, 10, 0) eq "ATGG", "substr_dna_seq 2");
ok(is_totally_masked("xnXX"), "is_totally_masked 1");
ok(is_totally_masked("AnXX") == 0, "is_totally_masked 2");
# my $p = packseq('GAGTGGCCATTACGGCCG');
# print STDERR Dumper(\$p);
# printf STDERR "\n%x\n", packseq('GAGTGGCCATTACGGCCG');
# ok(packseq('GAGTGGCCATTACGGCCG') == 37485790870, "packseq");
# ok(unpackseq(37485790870, 18) eq 'GAGTGGCCATTACGGCCG', "unpackseq");
ok(prot_3_to_1("ALA", "phe", 'ser', 'gly', 'tyr',
	       'pro', 'asn', 'asp', 'gln', 'GLU',
	       'LYS', 'ARG', 'TRP', 'HIS', 'CYS',
	       'MET', 'LEU', 'ILE', 'VAL', 'THR') eq 'AFSGYPNDQEKRWHCMLIVT',
   "prot_3_to_1");
ok(join(",",
	prot_1_to_3('AFSGYPNDQEKRWHCMLIVT')) eq
   join(",",
	"ALA", "PHE", 'SER', 'GLY', 'TYR',
	'PRO', 'ASN', 'ASP', 'GLN', 'GLU',
	'LYS', 'ARG', 'TRP', 'HIS', 'CYS',
	'MET', 'LEU', 'ILE', 'VAL', 'THR'),
   "prot_1_to_3");

ok(join("", na_2_prot("GACGCA", 1)) eq 'DA', "na_2_prot 1");
# Test debug variable in Sequtil

diag("\nThe next test should generate only a message about NAU being a bad codon:\n");

$Bio::Frescobi::Sequtil::debug = 0;
ok(join("", na_2_prot("ATGNNN", 1)) eq 'M?', "na_2_prot 2");
$Bio::Frescobi::Sequtil::debug = 1;
ok(join("", na_2_prot("ATGNAT", 1)) eq 'M?', "na_2_prot 3");

diag("\nThe next test should generate only one message about non-divisibility by 3:\n");
$Bio::Frescobi::Sequtil::debug = 0;
ok(join("", na_2_prot("ATGGATA", 1)) eq 'MD', "na_2_prot 4");
$Bio::Frescobi::Sequtil::debug = 1;
ok(join("", na_2_prot("ATGGATAC", 1)) eq 'MD', "na_2_prot 5");

my $colors = colorize_seq("WEREGREENWITHENVY");
ok(file2string("t/color.test") eq $colors . "\n", "colorize_seq");
ok(degenerate_index("ATCGATCG", "CGA") == 2, "degenerate_index 1");
ok(degenerate_index("ATCGNTCG", "GCT") == 3, "degenerate_index 2");
my $dna_seq = 'NNNATCGNNNNNNATATNNN';
my @ns = find_ns($dna_seq);
ok(scalar(@ns) == 3, "Number of N's");
ok($ns[0]->[0] == 1, "N 1,1");
ok($ns[0]->[1] == 3, "N 1,2");
ok($ns[1]->[0] == 8, "N 2,1");
ok($ns[1]->[1] == 13, "N 2,2");
ok($ns[2]->[0] == 18, "N 3,1");
ok($ns[2]->[1] == 20, "N 3,2");
$dna_seq = 'gcgnatat';
@ns = find_ns($dna_seq);
ok(scalar(@ns) == 1, "Number of N's (2)");
ok($ns[0]->[0] == 4, "N 1,1");
ok($ns[0]->[1] == 4, "N 1,2");

