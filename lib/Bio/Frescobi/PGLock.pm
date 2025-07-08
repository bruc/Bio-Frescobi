# This package provides a locking mechanism within a PostgreSQL database.

# The following table needs to be in the database for this mechanism to work:
# create table locks (
#        table_name   text,
#        pid	    text,
#        lock_time    timestamp with time zone,
#        host       text);

package Bio::Frescobi::PGLock;
use strict;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use Exporter;
use Sys::Hostname;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::CgPg;
use Carp;
@ISA = ('Exporter');
@EXPORT = qw();

my $lock_wait_time = 10;
my $lock_wait_cycles = 20000;	# About 3 days.

$VERSION = '0.01';

1;

=head1 NAME

Bio::Frescobi::PGLock -- Module containing database locking functions.

=head1 SYNOPSIS

 use Bio::Frescobi::PGLock;

=head2 Functions

 Bio::Frescobi::PGLock::lock_wait_time
 Bio::Frescobi::PGLock::lock_wait_cycles
 Bio::Frescobi::PGLock::grab_lock
 Bio::Frescobi::PGLock::free_lock

=head1 DESCRIPTION

PGLock provides a simple locking mechanism to control execution of
critical sections in parallel Perl code. It uses the exclusive locks
in the database engine. Each lock is identified by a table name, but
the table name is purely for identification purposes -- it is not
necessary for the table to exist in the database.

PGLock uses the C<locks> table in the database as well as the
transactional mechanism in the database engine to serialize this code
across multiple processes and hosts. The C<locks> table holds the
locks in use.

If a process using C<PGLock> attempts to grab a lock that's held by
another process, C<PGLock> enters a waiting loop defined by two
variables: the C<$lock_wait_time> which is how long it sleeps between
testing the lock, and the C<$lock_wait_cycles> which is how many
iterations the loop executes before giving up and dying with an
error. These two variables are accessed or changed with function
calls.

=head1 METHODS

=over 4

=item C<lock_wait_time([$new_wait_time])>

Returns the wait time between attempts to grab a lock. If a parameter
is passed to this function, then it is used to update the
C<$lock_wait_time> in units of seconds.

=cut

sub lock_wait_time {
    if (defined($_[0])) {
	$lock_wait_time = shift;
	if ($lock_wait_time !~ m/^\d+/) {
	    carp "Lock_wait_time value ($lock_wait_time) is not numeric.\n" .
		"It will be set to 10\n";
	    $lock_wait_time = 10;
	}
    }
    return $lock_wait_time;
}

=item C<lock_wait_cycles([$new_wait_cycles])>

Returns the number of iterations allowed for attempts to grab a
lock. If a parameter is passed to this function, then it is used to
update the C<$lock_wait_cycles>.

=cut

sub lock_wait_cycles {
    if (defined($_[0])) {
	$lock_wait_cycles = shift;
	if ($lock_wait_cycles !~ m/^\d+/) {
	    carp "Lock_wait_cycles value ($lock_wait_cycles) is not numeric.\n" .
		"It will be set to 20000.\n";
	    $lock_wait_cycles = 20000;
	}
	elsif ($lock_wait_cycles < 1) {
	    carp "Lock_wait_cycles value ($lock_wait_cycles) is too small.\n" .
		"It will be set to 1.\n";
	}
    }
    return $lock_wait_cycles;
}

=item C<grab_lock($pg, $table)>

Attempt to grab the lock on C<$table> that stored in the database
referred to by the database interface variable, C<$pg>. If the lock is
held by another process, then this function will enter a loop waiting
for the lock to be released. The loop will be executed
C<$lock_wait_cycles> times and it will sleep for C<$lock_wait_time>
seconds between each iteration.

=cut

sub grab_lock {
    my $pg = shift;
    my $table = shift;
    my ($count, $lock_table, $lock_pid, $lock_host);

    if (not $table) {
	croak "No table name specified for grab_lock.\n";
    }
    my $qtable = quotify($table);
    my $hostname = hostname();
    my $qhostname = quotify($hostname);
    for ($count = 1; $count <= $lock_wait_cycles; $count++) {
	if ($pg->driver eq 'SQLite') {
	    $pg->command("begin exclusive");
	}
	else {
	    $pg->command("begin");
	    $pg->command("lock table locks");
	}
	$lock_table = $pg->get_single_value("select table_name from locks " .
					    " where table_name = $qtable");
	if (not defined $lock_table) {
	    $pg->command("insert into locks (table_name, pid, lock_time, host) " .
			 "       values ($qtable, " .
			 "               '$$', " .
			 quotify($pg->curtime) .
			 ", $qhostname)");
	    $pg->command("end");
	    return;
	}
	else {
	    ($lock_pid, $lock_host) = $pg->get_single_row("select pid, host from locks " .
							  " where table_name = $qtable");
	    if ($lock_pid eq $$ and $lock_host eq $hostname) {
		croak "Table already locked by us!\n";
	    }
	    $pg->command("end");
	    sleep($lock_wait_time);
	}
    }
    die "Unable to grab $table lock after $lock_wait_cycles tries.\n";
}

=item C<free_lock($pg, $table[, $check_pid])>

Free mpt to grab the lock on C<$table> that stored in the database
referenced to by the database interface variable, C<$pg>. The option,
C<$check_pid>, specifies that this function will confirm that the lock
was held by this process on this system. The use of C<$check_pid> is
useful for confirming proper execution of lock management.

=cut

sub free_lock {
    my $pg = shift;
    my $table = shift;
    my $check_pid = shift;
    $check_pid = 1 if not defined $check_pid;
    my ($lock_pid, $lock_host);
     
    my $hostname = hostname();
    my $qhostname = quotify($hostname);
    if (not $table) {
	croak "No table name specified for free_lock.\n";
    }
    my $qtable = quotify($table);
    if ($pg->driver eq 'SQLite') {
	$pg->command("begin exclusive");
    }
    else {
	$pg->command("begin");
	$pg->command("lock table locks");
    }
    ($lock_pid, $lock_host) = $pg->get_single_row("select pid, host from locks " .
						  " where table_name = $qtable");
    if ($check_pid) {
	if ($lock_pid ne $$ or $hostname ne $lock_host) {
	    die "Attempt to free lock we do not own. Owner is $lock_pid and host is $lock_host. hostname() = $hostname\n";
	}
    }
    $pg->command("delete from locks " .
		 " where pid = '$lock_pid' " .
		 "   and host = $qhostname " .
		 "   and table_name = $qtable");
    $pg->command("end");
}

=back

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
