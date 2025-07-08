# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl CgPg.t'

# Tests to add:
# check the driver.
# run all tests using both Pg and SQLite

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 106);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::CgPg');
};
use Bio::Frescobi::Genutil;
use Data::Dumper;
use File::Basename;

#########################
foreach my $driver (("Pg", "SQLite")) {
    my $dbname;
    if ($driver eq "Pg") {
	$dbname = $ENV{"PG_TEST_DB"} || "regression";
    }
    else {
	# Put the temporary SQLite database in the "t" directory or
	# wherever this script is running.
	my $dirname = dirname($0);
	$dbname = "${dirname}/CgPg.sqlite";
    }
    my $pg = new Bio::Frescobi::CgPg(driver => $driver,
				     dbname => $dbname,
				     cgi => 0,
				     default_echo => 1);
    ok($pg, "open connection");
    BAIL_OUT("Unable to open connection to $dbname. Aborting test.\n") if not $pg;
    ok($pg->driver eq $driver, "Driver function");
    ok($pg->command("drop table if exists cgpg_test_table"), "Drop table");
    ok($pg->command("create table cgpg_test_table (" .
		    "   a     text , " .
		    "   i     integer, " .
		    "   f     float)"),
       "Create table");
    print STDERR Dumper(\$pg);
    printf STDERR '$0 is %s\n', $0;
    my $cmd;
    $cmd = "create index cgpg_test_table_a on cgpg_test_table";
    if ($driver eq 'Pg') {
	$cmd .= " using btree(a)";
    }
    else {
	$cmd .= " (a)";
    }
    ok($pg->command($cmd), "Create index");
    my @data = ( [ 'no', 3, 2e0 ],
		 [ 'man', 1, 7e-1 ],
		 [ 'can', 4, 1e-2 ],
		 [ 'live', 1, 8e-3 ], 
		 [ 'on', 5, 2e-4 ],
		 [ 'bread', 9, 8e-5 ],
		 [ 'alone', 2, 1e-6 ] );
    
    for (my $i = 1; $i <= scalar(@data); $i++) {
	my $row = $data[$i - 1];
	ok($pg->command("insert into cgpg_test_table (a, i, f) values ( " .
			join(", ", quotify($row->[0]), $row->[1], $row->[2]) . ")"),
	   "Insert row $i");
    }
    
    ok($pg->get_single_value("select a from cgpg_test_table where i = 3") eq 'no',
       "get_single_value");
    my @row = $pg->get_single_row("select a, i, f from cgpg_test_table where a = 'on'");
    ok(scalar(@row) == 3, "get_single_row properties.");
    ok(is_deeply(\@row, $data[4]), "get_single_row values.");
    my $rows = $pg->query("select a, i, f from cgpg_test_table where i = 1 order by a");
    ok($rows, "Query OK");
    ok($rows->[0]->[0] eq 'live', "Query data 1");
    ok($rows->[1]->[2] == 0.7, "Query data 1");
    ok(is_deeply(\@data,
		 $pg->get_all_rows("select a, i, f from cgpg_test_table order by f desc")),
       "get_all_rows");
    @row = $pg->get_array_for_field("select i from cgpg_test_table order by f desc");
    ok(is_deeply([3, 1, 4, 1, 5, 9, 2],
		 \@row),
       "get_array_for_field");
    my %h = $pg->get_hash_for_field("select a, i from cgpg_test_table");
    my %testh;
    foreach my $r (@data) {
	$testh{$r->[0]} = $r->[1];
    }
    ok(is_deeply(\%testh, \%h), "get_hash_for_field");
    ok($pg->table_exists("cgpg_test_table"), "table_exists 1");
    ok(! $pg->table_exists("foo_test_table"), "table_exists 2");
    if ($driver eq 'Pg') {
	@row = $pg->get_fields_for_table("cgpg_test_table");
	ok(&check_table_list($pg->get_tables_with_field('^a$')),
	   "get_tables_with_field");
	ok(is_deeply(\@row, [ "a",  "i",  "f" ]),
	   "get_fields_for_table");
    }
    ok(&check_table_list($pg->get_all_tables),
       "get_all_tables");
    @row = $pg->get_indexes('cgpg_test_table');
    print STDERR Dumper(@row);
    if ($driver eq 'Pg') {
	ok(uc($row[0]) =~ m/CREATE INDEX CGPG_TEST_TABLE_A ON .*CGPG_TEST_TABLE USING BTREE \(A\)/,
	   "get_indexes");
    }
    else {
	ok(uc($row[0]) =~ m/CREATE INDEX CGPG_TEST_TABLE_A ON CGPG_TEST_TABLE \(A\)/,
	   "get_indexes");
    }
    ok($pg->drop_indexes('cgpg_test_table'), "drop indexes");
    
    my @more_data = ( [ 'but',6, 8e-7 ],
		      [ 'dogs',5,  2e-8 ],
		      [ 'are',3,  8e-9 ],
		      [ 'fine',5,  4e-10 ],
		      [ 'with',8,  5e-11 ],
		      [ 'meat', 9,  9e-12 ]);
    my @lines = map { join("\t", @{$_}) . "\n"; } @more_data;
    ok($pg->copy_into_postgres("cgpg_test_table", \@lines),
       "copy_into_postgres 1");
    ok(abs($pg->get_single_value("select sum(f) from cgpg_test_table") - 2.718281828459) < 1.0e-13 ,
       "copy_into_postgres 2");
    print Dumper($pg->get_single_value("select sum(f) from cgpg_test_table"));
    if ($driver eq 'Pg') {
	ok($pg->command("drop sequence if exists testseq"),
	   "Drop sequence if needed");
	ok($pg->command("create sequence testseq start 1"),
	   "Create sequence");
	ok($pg->get_single_value("select nextval('testseq')") == 1,
	   "Prime the sequence");
    }
    else {
	ok($pg->command("drop table if exists testseq"),
	   "Drop sequence table");
	ok($pg->command("create table testseq (i integer)"),
	   "create sequence table");
	ok($pg->command("insert into testseq (i) values(1)"),
	   "initialize sequence");
    }
    ok($pg->nextval("testseq") == 2,
       "nextval");
    ok($pg->dbname eq $dbname, "dbname");
    ok($pg->default_echo(0) == 0, "default_echo");
    ok($pg->cgi(1) == 1, "cgi");
    ok($pg->die_on_error(1) == 1, "die_on_error");
}
ok(nullify("") eq '\N', "nullify 1");
ok(nullify("foo") eq 'foo', "nullify 2");
ok(bool_quotify(0) eq "'f'", "bool_quotify 1");
ok(bool_quotify(undef) eq "'f'", "bool_quotify 2");
ok(bool_quotify(1) eq "'t'", "bool_quotify 3");
ok(bool_quotify("") eq "'f'", "bool_quotify 4");
ok(bool_quotify("foo") eq "'t'", "bool_quotify 5");
ok(bool2perl("t"), "bool2perl 1");
ok(! bool2perl("f"), "bool2perl 2");
ok(! bool2perl(undef), "bool2perl 3");
ok(! defined(bool2perl('foo')), "bool2perl 4");
ok(perl2bool(1) eq 't', "perl2bool 1");
ok(perl2bool(0) eq 'f', "perl2bool 2");
print STDERR pg_array_join('a', 'b', '"c', '\d'), "\n";
ok(pg_array_join('a', 'b', '"c', '\d') eq '{a,b,"\"c","\\\\d"}', "pg_array_join");
&byte_compare(pg_array_join('a', 'b', '"c', '\d'),
	      '{a,b,"\"c","\\\\d"}');
my @row = pg_array_split('{e,f,g,1}');
ok(is_deeply(\@row, [ '{e', 'f', 'g', '1}']), "pg_array_split");

sub check_table_list {
    foreach my $e (@_) {
	if (ref($e) eq 'ARRAY') {
	    if ($e->[1] eq 'cgpg_test_table') {
		return 1;
	    }
	}
	elsif ($e eq 'cgpg_test_table') {
	    return 1;
	}
    }
    return 0;
}

sub byte_compare {
    my ($a, $b) = @_;
    printf STDERR "'%s'\n'%s'\n", $a, $b;
    
    if (length($a) != length($b)) {
	print STDERR "Length mismatch.\n";
	return 0;
    }
    else {
	for (my $i = 0; $i < length($a); $i++) {
	    my $cha = substr($a, $i, 1);
	    my $chb = substr($b, $i, 1);
	    if ($cha ne $chb) {
		printf STDERR "Mismatch at position $i: cha = (%s: %d)  chb = (%s: %d)\n",
		$cha, ord($cha), $chb, ord($chb);
	    }
	}
    }
}
__END__
  

__END__
  
