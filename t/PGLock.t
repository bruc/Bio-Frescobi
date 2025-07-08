# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl PGLock.t'

use warnings;
use strict;

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

# use Test::More (tests => 56);
use Test::More qw(no_plan);
BEGIN {
    use_ok('Bio::Frescobi::PGLock');
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
ok(Bio::Frescobi::PGLock::lock_wait_time(11) == 11,
   "lock_wait_time");
ok(Bio::Frescobi::PGLock::lock_wait_cycles(1) == 1,
   "lock_wait_cycles");
Bio::Frescobi::PGLock::grab_lock($pg, "test_table");
my $pid = $$;
ok($pg->get_single_value("select count(*) from locks where table_name = 'test_table'") == 1,
   "Confirm we have a lock");
ok($pg->get_single_value("select pid from locks where table_name = 'test_table'") == $pid,
   "Confirm PID");
ok($pg->command("update locks set pid = 0 where table_name = 'test_table'"),
   "Fake another process has the lock.");
diag("\nThe next test generates an expected error message:\n");
eval('Bio::Frescobi::PGLock::grab_lock($pg, "test_table");');
if ($@) {
    ok(1, "Lock worked.");
}
else {
    ok(0, "Lock failed.");
}

Bio::Frescobi::PGLock::free_lock($pg, "test_table", 0);
ok($pg->get_single_value("select count(*) from locks where table_name = 'test_table'") == 0,
   "Confirm lock is gone");
__END__
