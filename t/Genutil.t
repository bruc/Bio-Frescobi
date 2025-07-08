# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Genutil.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More (tests => 66);
# use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::Genutil');
};
use Data::Dumper;

#########################

ok(set_line_buffering(), "set_line_buffering");
system_with_check('date'); ok(1, "system_with_check");
diag("\nThe next test generates an expected error message:\n");
eval('system_with_check("foo_fah_no_command_here", 1);');
if ($@) {
    ok(1, "system_with_check: Failure detected");
}
else {
    ok(0, "system_with_check: Failure missed");
}
ok(commify("123456789") eq "123,456,789", "commify");

open (IN, "<t/Genutil.inp") ||
    fail("Unable to open t/Genutil.inp: $!");

{
    my $total_file = "";
    my @lines = file2array("t/Genutil.inp");
    foreach my $line (@lines) {
	my $file_line = <IN>;
	$total_file .= $file_line;
	chomp $file_line;
	ok($line eq $file_line, "Genutil file2array line: $line");
    }
    my $file_line = <IN>;
    if (defined $file_line) {
	ok(0, "Genutil file2array registration 1");
    }
    else {
	ok(1, "Genutil file2array registration 2");
    }
    ok ($total_file eq file2string("t/Genutil.inp"), "file2string");
}
my $cursor = 0;
my @scan_test = (3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5);
scanarray(\@scan_test, $cursor, '1', 0);
ok($cursor == 1, "scanarray 1");
incrarray(\@scan_test, $cursor, 2);
scanarray(\@scan_test, $cursor, '1', 0);
ok($cursor == 3, "scanarray 2");
incrarray(\@scan_test, $cursor, 2);
ok($cursor == 5, "incrarray 1");
scanarray(\@scan_test, $cursor, '5', 0);
ok($cursor == 8, "scanarray 3");

diag("\nThe next test should print 'At <date>: TeSTed\n");
print STDERR "\n";
dated_mesg("TeSTed");

ok(max(2,3,4) == 4, "max");
ok(min(-4, -5, -7) == -7, "min");
ok(abs(2-4) == 2, "abs");
ok(sum(2,2) == 4, "sum");
my @strings = ('this', 'is', 'hard', 'for', 'computers');
# The following logic is used to avoid locale issues.
ok(minst(@strings) eq (sort {$a cmp $b} @strings)[0], "minst");
ok(maxst(@strings) eq (sort {$a cmp $b} @strings)[$#strings], "maxst");
my @piper_output = piper("grep r", \@strings);
ok(shift(@piper_output) == 0, "piper status");
my $piper_results = join("", @piper_output);
ok($piper_results eq "hard\nfor\ncomputers\n", "piper results");
my $trim = " middle word \t ";
ok(trim($trim) eq "middle word", "trim");
my $squeeze = " ah one and ah      two  ";
squeeze_blanks($squeeze);
ok($squeeze eq "ah one and ah two", "squeeze_blanks");
ok(join(" ", permute("on", 3)) eq 'ooo oon ono onn noo non nno nnn', 'permute');
my $where = "";
add_to_where($where, "base = 'A'", "and");
add_to_where($where, "position between 1000 and 1100", "and");
ok($where eq "base = 'A' and position between 1000 and 1100", "add_to_where");
ok(sprint1f_undef('%4d', 2718) . " " . sprint1f_undef('%3d', undef) eq '2718    ',
   "sprint1f_undef");
ok(string_truth("on") == 1, "string_truth 1");
ok(string_truth("yes") == 1, "string_truth 2");
ok(string_truth("1") == 1, "string_truth 3");
ok(string_truth("no") == 0, "string_truth 4");
ok(string_truth("off") == 0, "string_truth 5");
ok(string_truth("0") == 0, "string_truth 6");
ok(string_truth(undef, 1) == 1, "string_truth 7");
my %test = ("A" => 29, "B" => 42);
my %lctest = lc_hash_keys(%test);
ok($lctest{"b"} == 42, "lc_hash_keys");
ok(bool2logical(undef) == 0, "bool2logical 1");
ok(bool2logical('f') == 0, "bool2logical 2");
ok(bool2logical('t') == 1, "bool2logical 3");
ok(logical2bool(1) eq 't', "logical2bool 1");
ok(logical2bool(0) eq 'f', "logical2bool 2");
ok(quotify('string') eq "'string'", "quotify easy");
ok(quotify('\a " ' . "'" . ' b') eq q#'\\\\a " '' b'#, "quotify hard");
ok(find_smallest(3,2,4,1,5) == 3, "find_smallest");
ok(commify(123456789012) eq "123,456,789,012", "commify");
ok(join(",", clean_undef(1, undef, 2)) eq '1,,2', "clean_undef");
ok(join(",", mark_undef(1, undef, 2)) eq '1,undef,2', "mark_undef");
ok(find_num_in_array(5, [1,2,5,99]) == 2, "find_num_in_array 1");
ok(find_num_in_array(9, [1,2,5,99]) == -1, "find_num_in_array 2");
ok(find_string_in_array('fox', ['the', 'quick', 'brown', 'fox', 'jumped', 'over', 'the', 'lazy', 'dog']) == 3, "find_string_in_array 1");
ok(find_string_in_array('nevermore', ['Edgar', 'Allan', 'Poe']) == -1, "find_string_array 2");
my $line = "line\r\n";
dos_chomp($line);
ok($line eq 'line', "dos_chomp");
ok(md5sum("t/Genutil.inp") eq '131eec460e4391a552e7027616131858', "md5sum");
my @parse = parse_ms_config("t/Genutil.ini");
ok(join(',', @{$parse[0]}) eq 'header,Hamlet,', "parse_ms_config 1");
ok(join(',', @{$parse[1]}) eq 'body,To be or not to be.,', "parse_ms_config 2");
ok(join(',', @{$parse[2]}) eq 'trailer,by Shakespeare', "parse_ms_config 3");
ok(scalar(@parse) == 3, "parse_ms_config 4");
ok(range_overlap(1,3,1,3) == 1, "range_overlap 1");
ok(range_overlap(1,3,5,6) == 0, "range_overlap 2");
ok(range_overlap(1,3,-3,-1) == 0, "range_overlap 3");
ok(range_overlap(1,3,-3,2) == 1, "range_overlap 4");
ok(range_overlap(1,3,-3,4) == 1, "range_overlap 5");
ok(range_overlap(1,5,2,4) == 1, "range_overlap 6");
ok(range_overlap(1,5,2,6) == 1, "range_overlap 7");
ok(range_overlap(1,5,6,8) == 0, "range_overlap 8");
ok(range_overlap(1,5,5,8) == 1, "range_overlap 9");
ok(range_overlap(10,20,5,10) == 1, "range_overlap 10");

# No tests for display_environment or display_variables.
