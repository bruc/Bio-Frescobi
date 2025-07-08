# This package contains simplified calls to the PostgreSQL database interface using
# the DBI and DBD::Pg driver.
#
# The big change is that the old query function is now split into
# query (for Select's that return data) and command (for everything
# else).

package Bio::Frescobi::CgPg;
use strict;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use DBI;
use DBD::Pg qw(:async :pg_limits :pg_types);
use DBD::SQLite;
use Bio::Frescobi::Genutil;
use Exporter;
use Carp;

# The following is unfortunate -- we need to specify these names ourselves.
# These were taken from the Postgres installation include/libpq-fe.h

use constant {
	      PGRES_EMPTY_QUERY => 0,	   # empty query string was executed 
	      PGRES_COMMAND_OK => 1,	   # a query command that doesn't return
					   # anything was executed properly by the
					   # backend 
	      PGRES_TUPLES_OK => 2,	   # a query command that returns tuples was
				           # executed properly by the backend, PGresult
					   # contains the result tuples 
	      PGRES_COPY_OUT => 3,	   # Copy Out data transfer in progress 
	      PGRES_COPY_IN => 4,	   # Copy In data transfer in progress 
	      PGRES_BAD_RESPONSE => 5,	   # an unexpected response was recv'd from the
					   # backend 
	      PGRES_NONFATAL_ERROR => 6,   # notice or warning message 
	      PGRES_FATAL_ERROR => 7,	   # query failed 
	      PGRES_COPY_BOTH => 8,	   # Copy In/Out data transfer in progress 
	      PGRES_SINGLE_TUPLE => 9,     # single tuple from larger resultset
	      PGRES_PIPELINE_SYNC => 10,   # pipeline synchronization point 
	      PGRES_PIPELINE_ABORTED => 11 # Command didn't run because of an abort
	      };

@ISA = ('Exporter');
@EXPORT = qw(PGRES_EMPTY_QUERY
	     PGRES_COMMAND_OK 
	     PGRES_TUPLES_OK 
	     PGRES_COPY_OUT
	     PGRES_COPY_IN 
	     PGRES_BAD_RESPONSE 
	     PGRES_NONFATAL_ERROR
	     PGRES_FATAL_ERROR 
	     PGRES_COPY_BOTH 
	     PGRES_SINGLE_TUPLE 
	     PGRES_PIPELINE_SYNC
	     PGRES_PIPELINE_ABORTED
	     &nullify
	     &bool_quotify
	     &bool2perl
	     &perl2bool
	     &pg_array_join
	     &pg_array_split);

$VERSION = '0.02';

my %flag_components = (DEFAULT_ECHO => 1,
		       CGI => 1,
		       DIE_ON_ERROR => 1);

1;


=head1 NAME

Bio::Frescobi::CgPg -- class for convenient access to PostgreSQL or SQLite

=head1 SYNOPSIS

 use Bio::Frescobi::CgPg;

=head2 Constructor
    
 $pg = Bio::Frescobi::CgPg->new(option => value, option => value, ...);

=head2 Object Methods
    
 $pg->close_connection;
 $pg->query($expected_status, $sql_statement[, $echo]);
 $pg->command($sql[, $echo]);
 $pg->reset;
 $pg->get_single_value($sql_statement[, $echo]);
 $pg->get_single_row($sql_statement[, $echo]);
 $pg->get_all_rows($sql_statement[, $echo]);
 $pg->get_array_for_field($sql_statement[, $echo]);
 $pg->get_hash_for_field($sql_statement[, $check_uniqueness[, $echo]]);
 $pg->table_exists($table_name[, $echo]);
 $pg->get_fields_for_table($table_name);
 $pg->get_tables_with_field($field_pattern);
 $pg->get_all_tables;
 $pg->get_indexes($table_name);
 $pg->drop_indexes($table_name);
 $pg->startcopy($table);
 $pg->putline(@lines);
 $pg->endcopy;
 $pg->copy_into_postgres($table, $datap);
 $pg->default_echo[($echo_option)];
 $pg->die_on_error[($value)];
 $pg->cgi[($value)];
 $pg->autocommit[($value)];

=head2 Convenience Functions
    
 &nullify($string)
 &bool_quotify($val)
 &bool2perl($val)
 &perl2bool($val)
 &pg_array_join(@list)
 &pg_array_split($string)

=cut

sub new {
    my ($pkg, $class, $self);

=head1 CONSTRUCTOR

The constructor for this object is C<new> and
it expects a hash argument list, i.e. option value pairs.

 new({option => value, option => value})

The constructor opens a new connection to a PostgreSQL database. If the connection
is successful, then a new object is returned. If it fails and C<die_on_error> is off,
an undef is returned. Options are
specified using name, value pairs as in a hash. In addition
some of the options can be set from the environment, although
anything in the arguments to C<new> will override environmental
settings.

=head2 Options:

=over 4

=item C<driver> = DBI driver to use.

Only two choices are provided currently: Pg or SQLite. The default is Pg (Postgres).

=item C<dbname> = Database name

Specify the PostgreSQL database name or the SQLite database file. Can also be specified
by the PGDATABASE environment variable when Postgres is being used.j

=item C<cgi> = C<yes | no | on | off | 0 | 1>

Specify whether the script is run in a CGI environment. If not specified,
the module looks for the existence of the C<HTTP_USER_AGENT> environment
variable, and if found, then CGI is presumed.

=item C<default_echo> = C<yes | no | on | off | 0 | 1>

Specify whether queries are echoed to STDERR by default. Off by default.

=item C<die_on_error> = C<yes | no | on | off | 0 | 1>

Specify whether the module should terminate execution if
the connection fails. On by default.

=item C<host> = Server host

Specify name or address of PostgreSQL server. This value can also
be read from the C<PGHOST> environment variable.

=item C<port> = Server port

Specify the port for the PostgreSQL server. Connections to different
versions of PostgreSQL or different server environments are specified
through the port number. You must check the PostgreSQL setup script
to identify the specific port. 
This value can be read from the C<PGPORT> environment variable.

=item C<user> = User name

Specify the username for the connection. Can be specified by the C<PGUSER>
environment variable.

=item C<password> = password

Specify the password for the connection. Can be specified by the C<PGPASSWORD>
environment variable.

=item C<options> = options

Specify addition options for the connection. Can be specified by the C<PGOPTIONS>
environment variable.

=back

=cut

    $class = shift;
    eval {($pkg) = caller(0);};
    if ($class ne $pkg) {
	unshift @_, $class;
    }
    $self = {};
    bless $self;
    my %args = lc_hash_keys(@_);
    my $msgs = "";
    if (not exists($args{dbname})) {
	if (exists($ENV{PGDATABASE})) {
	    $args{dbname} = $ENV{PGDATABASE};
	}
	else {
	    $msgs .= "No database name specified.\n";
	}
    }
    $self->{DBNAME} = $args{dbname};
    if (not exists($args{driver})) {
	$args{driver} = "Pg";
    }
    $self->{DRIVER} = $args{driver};
    if ($self->{DRIVER} !~ m/^(Pg|SQLite)$/) {
	croak sprintf("Unsupported driver, '%s', specified.\n", $self->{DRIVER});
    }
    $self->{CGI} = string_truth($args{cgi},
				exists($ENV{"HTTP_USER_AGENT"}),
				"cgi");
    $self->{DEFAULT_ECHO} = string_truth($args{default_echo},
					 0,
					 "default_echo");
    $self->{DIE_ON_ERROR} = string_truth($args{die_on_error},
					 1,
					 "die_on_error");
    if ($msgs) {
	if ($self->{DIE_ON_ERROR}) {
	    croak $msgs;
	}
	else {
	    carp $msgs;
	    return undef;
	}
    }
    my $username = "";
    my $password = "";
    if ($self->{DRIVER} eq 'Pg') {
	my $connect_string = "dbname=$args{dbname}";
	foreach my $option (("host", "port", "options", "user", "password")) {
	    if (exists($args{$option})) {
		if ($option eq 'user') {
		    $username = $args{$option};
		} elsif ($option eq 'password') {
		    $password = $args{$option};
		} else {
		    $connect_string .= sprintf(";%s=%s", $option, $args{$option});
		}
	    } elsif (exists($ENV{"PG" . uc($option)})) {
		if ($option eq 'user') {
		    $username = $ENV{"PG" . uc($option)}
		} elsif ($option eq 'password') {
		    $password = $ENV{"PG" . uc($option)}
		} else {
		    $connect_string .= sprintf(";%s=%s", $option, $ENV{"PG" . uc($option)});
		}
	    }
	}
	$self->{CONNECTION_STRING} = "dbi:Pg:$connect_string";
    }
    else {
	$self->{CONNECTION_STRING} = sprintf("dbi:SQLite:dbname=%s", $self->{DBNAME});
    }
    $self->{USERNAME} = $username;
    $self->{PASSWORD} = $password;
    $self->connect;
    my $version;
    if ($self->{DRIVER} eq 'Pg') {
	my $version_line = $self->get_single_value("select version()");
	($version) = ($version_line =~ m/^PostgreSQL (\S+)\s/);
	if ($version =~ m/^7\.[012]/) {
	    $self->command("set enable_seqscan = 'off'");
	}
    }
    else {
	$version = $DBD::SQLite::sqlite_version;
    }
    $self->{VERSION} = $version;
    $self->{INCOPY} = 0;
    $self->{AUTOCOMMIT} = 1;
    return $self;
}

=head1 METHODS

The following methods are provided:

=over 4

=cut


sub close_connection {
    my $self = shift;

=item C<close_connection>

Closes the connection. There is no way to reopen the connection except
for recreating the object.

=cut
    my $dbh = $self->{DBH};
    if ($dbh->disconnect) {
	$self->{DBH} = undef;
    }
    else {
	if ($self->{DIE_ON_ERROR}) {
	    croak "Error detected in CgPg: " . $dbh->errstr . "\n";
	}
	else {
	    carp "Warning detected in CgPg: " . $dbh->errstr . "\n";
	}
    }
}

sub query {
    my $self = shift;

=item C<query($sql_statement[, $echo])>

Execute the C<$sql_statement> and return the resulting array reference to the data.
The C<$echo> variable overrides the default echo setting, and
controls whether the C<$sql_statement> is echoed to C<STDERR> prior to
execution.

=cut
    my ($sql_statement, $echo) = @_;
    my ($result, $result2, $msg, $break);
    local $/ = "\n";

    $self->_handle_echo($sql_statement, $echo);
    my $dbh = $self->{DBH};
    my $sth = $dbh->prepare($sql_statement);
    $self->{STH} = $sth;
    if (not $sth) {
	$msg = sprintf("Command %s returned status %d (%s)\n",
		       $sql_statement,
		       $dbh->err,
		       $dbh->errstr);
	$dbh->do("rollback transaction", { RaiseError => $self->{DIE_ON_ERROR},
					   PrintError => 1});
	return undef;
    }
    my $rv = $sth->execute;
    if (not $rv) {
	$msg = sprintf("Command %s returned status %d (%s)\n",
		       $sql_statement,
		       $sth->err,
		       $sth->errstr);
	$dbh->do("rollback transaction", { RaiseError => $self->{DIE_ON_ERROR},
					   PrintError => 1});
	return undef;
    }
    return $sth->fetchall_arrayref;
}

sub _handle_echo {
    my $self = shift;
    my $sql = shift;
    my $echo = shift;

    if (not defined($echo)) {
	$echo = $self->{DEFAULT_ECHO};
    }
    my $break = "";
    if ($self->{CGI}) {
	$break = "<br>";
    }
    if ($echo) {
	my $date = &date;
	print STDERR "At $date: executing ", $sql, "$break\n";
    }
}

sub command {
    my $self = shift;
    my $sql = shift;
    my $echo = shift;

=item C<command($sql[, $echo])>

Execute an SQL command and checks for success. The C<die_on_error>
setting will determine the outcome.

=cut

    $self->_handle_echo($sql, $echo);
    my $dbh = $self->{DBH};
    my $rv = $dbh->do($sql);
    if (not defined $rv) {
	my $msg = sprintf("Command $sql returned status %d (%s)\n",
			  $dbh->err,
			  $dbh->errstr);
	croak $msg if $self->{DIE_ON_ERROR};
    }
    return $rv;
}
sub DESTROY {
    my $self = shift;
    if ($self->{INCOPY}) {
	$self->endcopy;
    }
    $self->close_connection;
}

sub get_single_value {
    my $self = shift;

=item C<get_single_value($sql_statement[, $echo])>

Execute the SQL statement and return the first element
of the first record. If there is more than one field or more than one column,
the method will complain, but still return the first element of the first record.
If nothing is returned, then this method will return undef.

=cut

    my ($sql_statement, $echo) = @_;

    my $rows = $self->query($sql_statement,
			    $echo);
    if (scalar(@{$rows}) == 0) {
	return undef;
    }
    if (scalar(@{$rows}) != 1) {
	carp "get_single_value returned multiple rows.";
    }
    my @row = @{$rows->[0]};
    if (scalar(@row) != 1) {
	carp "get_single_value returned multiple columns.";
    }
    return $row[0];
}

sub get_single_row {
    my $self = shift;

=item C<get_single_row($sql_statement[, $echo])>

Execute an SQL statement and return an
array containing the first record. If multiple rows are
returned, the method will complain, but otherwise return the array.
If there is no data, the empty array will be returned.

=cut

    my ($sql_statement, $echo) = @_;

    my $rows = $self->query($sql_statement,
			    $echo);
    if (scalar(@{$rows}) == 0) {
	my @ret = ();
	return @ret;;
    }
    if (scalar(@{$rows}) != 1) {
	carp "get_single_row returned multiple rows.";
    }
    return @{$rows->[0]};
}

sub get_all_rows {
    my $self = shift;

=item C<get_all_rows($sql_statement[, $echo])>

Execute an SQL statement and return a reference
to an array containing references to arrays for each of the records.

=cut

    my ($sql_statement, $echo) = @_;
    my ($result, @ret, $i);

    my $rows = $self->query($sql_statement,
			   $echo);
    return $rows;
}

sub get_array_for_field {
    my $self = shift;

=item C<get_array_for_field($sql_statement[, $echo])>

Execute an SQL statement and return an array
containing all the elements from the first field of all the records.

=cut

    my ($sql_statement, $echo) = @_;
    my $rows = $self->query($sql_statement,
			    $echo);
    my @ret;

    foreach my $row (@{$rows}) {
	push (@ret, $row->[0]);
    }
    return @ret;
}

sub get_hash_for_field {
    my $self = shift;

=item C<get_hash_for_field($sql_statement[, $check_uniqueness[, $echo]])>

Execute an SQL statement and return a hash whose keys are all the
elements from the first field of all the records. If the SQL returns
just one field, then the value of the hash will be the number 1.  If
the SQL returns two fields, then the value of the hash will be the
second field.  If the SQL returns more than two fields, then the value
of the hash will be a pointer to an array containing all the fields
returned after the first.

The parameter, C<$check_uniqueness>, is a flag causing this method to
display warnings whenever a duplicate key is found.

=cut

    my ($sql_statement, $check_uniqueness, $echo) = @_;
    my %ret = ();

    my $rows = $self->query($sql_statement,
			    $echo);

    foreach my $row (@{$rows}) {
	my $key = $row->[0];
	if ($check_uniqueness) {
	    if (exists($ret{$key})) {
		warn "get_hash_for_field: Duplicate key: $key found\n";
	    }
	}
	if (scalar(@{$row}) == 1) {
	    $ret{$key} = 1;
	}
	elsif (scalar(@{$row}) == 2) {
	    $ret{$key} = $row->[1];
	}
	else {
	    for (my $j = 1; $j < scalar(@{$row}); $j++) {
		$ret{$key}->[$j-1] = $row->[$j];
	    }
	}
    }
    return %ret;
}

sub table_exists {
    my $self = shift;

=item C<table_exists($table_name[, $echo])>

Checks for the existence of $table_name. 

=cut

    my ($table_name, $echo) = @_;

    my $sth = $self->{DBH}->table_info();
    while (my $hash = $sth->fetchrow_hashref) {
	if ($hash->{TABLE_TYPE} eq 'TABLE' and
	    $hash->{TABLE_NAME} eq $table_name) {
	    return 1;
	}
    }
    return 0;
}

sub get_fields_for_table {
    my $self = shift;

=item C<get_fields_for_table($table_name)>

Get all the fields for a table in their proper order. This
function depends upon the specific implementation of the catalog
in PostgreSQL and could break in future releases.

=cut

    my ($table) = @_;
    my @ret = ();
    my $i = 1;

    croak "Table $table does not exist.\n" unless $self->table_exists($table);
    my $sth = $self->{DBH}->column_info(undef, undef, $table, undef);
    while (my $hash = $sth->fetchrow_hashref) {
	push (@ret, $hash->{COLUMN_NAME});
	if ($i != $hash->{ORDINAL_POSITION}) {
	    croak sprintf("Table $table did not return rows in ordinal order. Last column is %s\n",
			  $hash->{COLUMN_NAME});
	}
	$i++;
    }
    return @ret;
}

sub get_tables_with_field {
    my $self = shift;

=item C<get_tables_with_field($field_pattern)>

Return a list of tables which contain fields that match the given
pattern. The list contains two element arrays where the first element
is the field name, and the second is the table name. This
function depends upon the specific implementation of the catalog
in PostgreSQL and could break in future releases.

=cut

    my ($fieldpat) = @_;
    my @ret = ();
    my ($i, $field, $table);

    if ($self->{DRIVER} ne 'Pg') {
	carp "get_tables_with_field is only implemented for Postgres.\n";
	return @ret;
    }
    my $rows = $self->get_all_rows("select a.attname, c.relname " .
				   "  from pg_class c, pg_attribute a " .
				   " where a.attnum > 0 and " .
				   "       a.attrelid = c.oid and " .
				   "       a.attname ~ " . quotify($fieldpat) . " and " .
				   "       c.relkind = 'r'");
    foreach my $row (@{$rows}) {
	my ($field, $table) = @{$row};
	push (@ret, [ $field, $table ]);
    }
    return @ret;
}

sub get_all_tables {
    my $self = shift;

=item C<get_all_tables>

Return a list of all user table names. This
function depends upon the specific implementation of the catalog
in PostgreSQL and could break in future releases.

=cut

    my ($i, $table);
    my @ret;
    if ($self->{DRIVER} eq 'Pg') {
	@ret = $self->get_array_for_field("select tablename " .
					  "  from pg_tables " .
					  " where schemaname not in ('information_schema', 'pg_catalog')");
    }
    else {
	my $sth = $self->{DBH}->table_info();
	while (my $hash = $sth->fetchrow_hashref) {
	    if ($hash->{TABLE_TYPE} eq 'TABLE') {
		push(@ret, $hash->{TABLE_NAME});
	    }
	}
    }
    return @ret;
}

sub get_all_table_info {
    my $self = shift;

=item C<get_all_table_info>

Return all the data from the DBI function C<table_info> as an array of hashrefs.

=cut

    my @ret;

    my $sth = $self->{DBH}->table_info();
    while (my $hash = $sth->fetchrow_hashref) {
	push(@ret, $hash);
    }
    return @ret;
}

sub get_indexes {
    my $self = shift;

=item C<get_indexes($table[, $echo])>

Return a list of all index creation commands for a table whose name
is C<$table>. Echoing is controlled by the $echo variable.
The result is returned as a list of 
commands necessary to create the indexes. This
function depends upon the specific implementation of the catalog in
PostgreSQL and could break in future releases.

=cut

    my $table = shift;
    my $echo = shift;
    
    my @ret;
    if ($self->{DRIVER} eq 'Pg') {
	@ret = $self->get_array_for_field("select indexdef " .
					  "  from pg_indexes " .
					  " where tablename = " . quotify($table),
					  $echo);
    }
    else {
	my @info = $self->get_all_table_info;
	foreach my $hash (@info) {
	    if ($hash->{"TABLE_TYPE"} eq "INDEX") {
		push(@ret, $hash->{"sqlite_sql"});
	    }
	}
    }
    return @ret;
}

sub drop_indexes {
    my $self = shift;
    
=item C<drop_indexes($table_name[, $echo])>

Delete all the indexes associated with table whose name is
C<$table_name>. Echoing is controlled by the C<$echo> variable.
This function depends upon the specific implementation
of the catalog in PostgreSQL and could break in future releases.

=cut

    my $table = shift;
    my $echo = shift;

    # Force the database to die if we get an error in here. Restore state afterwards.
    
    my $die_on_error = $self->die_on_error;
    $self->die_on_error(1);
    my @indices;
    if ($self->{DRIVER} eq 'Pg') {
	@indices = $self->get_array_for_field("select indexname " .
					      "  from pg_indexes " .
					      " where tablename = " . quotify($table),
					      $echo);
    }
    else {
	my @info = $self->get_all_table_info;
	foreach my $hash (@info) {
	    if ($hash->{"TABLE_TYPE"} eq "INDEX") {
		my $index_name = (split(/\s+/, $hash->{"sqlite_sql"}))[2];
		push(@indices, $index_name);
	    }
	}
    }
    foreach my $indexname (@indices) {
	$self->command("DROP INDEX $indexname", $echo);
    }
    $self->die_on_error($die_on_error);
    return 1;
}

sub startcopy {
    my $self = shift;

=item C<startcopy($table)>

Initiate a PostgreSQL COPY command for $table.

=cut

    my ($table) = @_;

    croak "No table specified\n" if $table eq "";
    croak "COPY operation already in progress\n" if $self->{INCOPY};
    $self->{INCOPY} = 1;
    $self->{INCOPY_TABLE} = $table;
    $self->{DBH}->{AutoCommit} = 0; # This change is intended to be
                                    # temporary, so the object setting
                                    # is left alone.
    if ($self->{DRIVER} eq "Pg") {
	$self->{DBH}->do("copy $table from stdin");
    }
    else {
	$self->{FIELDS} = [ $self->get_fields_for_table($table) ];
    }
}

sub putline {
    my $self = shift;

=item C<putline(@lines)>

Load a list of lines into the database as part of a COPY command.

=cut

    my $line;
    croak "COPY not in progress\n" unless $self->{INCOPY};
    my $dbh = $self->{DBH};
    if ($self->{DRIVER} eq 'Pg') {
	foreach $line (@_) {
	    $dbh->pg_putcopydata($line);
	}
    }
    else {
	my $cmd = sprintf("INSERT INTO %s (%s) VALUES ",
			  $self->{INCOPY_TABLE},
			  join(", ", @{$self->{FIELDS}}));
	for (my $i = 1; $i <= scalar(@_); $i++) {
	    my $line = $_[$i-1];
	    chomp($line);
	    my @data;
	    foreach my $datum (split(/\t/, $line)) {
		if ($datum eq '\N') {
		    $datum = 'NULL';
		}
		else {
		    $datum = $dbh->quote($datum);
		}
		push(@data, $datum);
	    }
	    $cmd .= sprintf("( %s )", join(", ", @data));
	    if ($i < scalar(@_)) {
		$cmd .= ", ";
	    }
	}
	$self->command($cmd, 0);
    }
    return 1;
}

sub endcopy {
    my $self = shift;

=item C<endcopy>

Ends a COPY operation.

=cut

    my ($status);

    if (not $self->{INCOPY}) {
	carp "No COPY operation in progress\n" if not $self->{INCOPY};
	return;
    }
    if ($self->{DRIVER} eq "Pg") {
	$self->{DBH}->pg_putcopydata("\\.\n");
	if (not $self->{DBH}->pg_putcopyend) {
	    croak sprintf("putcopyend failed. Code = %d\n", $status);
	}
    }
    $self->{INCOPY} = 0;
    $self->{DBH}->{AutoCommit} = $self->{AUTOCOMMIT}; # Reset it.
    return 1;
}

sub copy_into_postgres {
    my $self = shift;

=item C<copy_into_postgres($table, $datap)>

Copies an array of input line for the COPY command referred to
by the $datap argument into $table.

=cut

    my ($table, $datap) = @_;

    $self->startcopy($table);
    $self->putline(@$datap);
    $self->endcopy();
}

sub dbname {
    my $self = shift;

=item C<dbname>

Return the database name for the CgPg object.

=cut

    return $self->{DBNAME};
}

sub driver {
    my $self = shift;

=item C<driver>

Return the driver for the CgPg object.

=cut

    return $self->{DRIVER};
}

sub autocommit {
    my $self = shift;

=item C<autocommit>

Returns the current autocommit mode and
optionally sets it to the argument you specify.

=cut

    my $arg = shift;

    if (defined $arg) {
	$self->{DBH}->{AutoCommit} = $arg;
	$self->{AUTOCOMMIT} = $arg;
    }
    return $self->{AUTOCOMMIT};
}

sub nextval {
    my $self = shift;
    my $table = shift;

=item C<nextval($table)>

Returns the next value of a sequence. In Postgres, this is a database function,
but in SQLite, we need to update the table.

=cut
    my $value;
    
    if ($self->{DRIVER} eq 'Pg') {
	$value = $self->get_single_value("select nextval('$table')");
    }
    else {
	$self->{DBH}->{AutoCommit} = 0;
	$self->command("begin exclusive");
	$self->command("update $table set i = i + 1");
	$value = $self->get_single_value("select i from $table");
	$self->command("end");
	$self->{DBH}->{AutoCommit} = $self->{AUTOCOMMIT};
    }
    return $value;
}

sub curtime {
    my $self = shift;

=item C<curtime()>

Returns the current time appropriately formatted for the database driver. 

=cut
    my $curtime;

    if ($self->driver eq 'Pg') {
	$curtime = $self->get_single_value("select 'now'::timestamp with time zone");
    }
    else {
	$curtime= $self->get_single_value("select datetime('now', 'localtime')");
    }
    return $curtime;
}

# documentation for autoloaded methods goes here.

=item C<default_echo[($echo_option)]>

Returns the current default_echo mode and
optionally sets it to the argument you specify.

=item C<die_on_error[($value)]>

Returns the setting for die_on_error mode,
and optionally sets it to the argument you specify.

=item C<cgi[($value)]>

Returns the setting for cgi mode,
and optionally sets it to the argument you specify.

=cut

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self);
    my $pkg = __PACKAGE__;
    my ($pos, $new_val);

    croak "$self is not an object\n" if not $type;
    croak "$pkg AUTOLOAD function fails on $type\n"
	if $type !~ m/^$pkg$/;
    my $name = uc($AUTOLOAD);
    $name =~ s/^.*:://;
    if (exists($flag_components{$name})) {
	if (not exists($self->{$name})) {
	    croak "Logic error: No data found for $name\n";
	}
	$new_val = shift;
	if (defined $new_val) {
	    $self->{$name} = string_truth($new_val,
					  $self->{$name},
					  $name);
	}
	return $self->{$name};
    }
    else {
	if (not exists($self->{$name})) {
	    croak "$name is not a valid method for $type\n";
	}
	elsif (defined $_[0]) {
	    croak "No permission to change value of $name\n";
	}
	else {
	    return $self->{$name};
	}
    }
}

=back

=head1 CONVENIENCE FUNCTIONS

=over

=cut

sub nullify {
    if ($_[0] eq "") {
	return "\\N";
    }
    else {
	return $_[0];
    }

=item C<nullify($string)>

Converts a blank argument into PostgreSQL NULL strings, i.e. \N, otherwise
just returns the argument.

=cut

}

sub bool_quotify {
    if ($_[0]) {
	return "'t'";
    }
    else {
	return "'f'";
    }

=item C<bool_quotify($val)>

Converts a Perl boolean value into
a quoted PostgreSQL boolean value ("'t'" or "'f'").
The result is three characters in length.

=cut


}

sub bool2perl {
    my $val = shift;
    if (not defined($val)) {
	return 0;
    }
    elsif ($val eq 't') {
	return 1;
    }
    elsif ($val eq 'f') {
	return 0;
    }
    else {
	carp "Bad bool value ($val)\n";
	return undef;
    }

=item C<bool2perl($val)>

Converts a PostgreSQL boolean value ('t' or 'f') into
a Perl boolean value (1 or 0). If the argument is not either 't' or 'f',
then C<undef> is returned.

=cut

}

sub perl2bool {
    my $val = shift;
    if ($val) {
	return 't';
    }
    else {
	return 'f';
    }

=item C<perl2bool($val)>

Converts a Perl boolean value into a PostgreSQL boolean value ('t' or 'f').
The value is not quoted on return as would be the case for bool_quotify.

=cut

}

sub pg_array_join {
    my @args = @_;
    my @pieces;
    foreach my $arg (@args) {
	if ($arg =~ m/[ {},"\\]/) { # " ] / ) {
	    $arg =~ s/(["\\])/\\$1/g; # " ] ) ;
	    $arg = '"' . $arg . '"';
	}
	push (@pieces, $arg);
    }
    return "{" . join(",", @pieces) . "}";
}

=item C<pg_array_join(@list)>

Converts a list of values into a string suitable for use in any array
input string.

=cut

sub pg_array_split {
    my $arg = shift;

    if ($arg =~ m/^\{(.*)\}$/) {
	my $core = $1;
	if ($core =~ m/"/) { # " Keep Emacs happy.
	    carp "String $arg will not be split by pg_array_split correctly.\n";
	}
	return split(/,/, $arg);
    }
    else {
	croak "Malformed array string: $arg\n";
    }
}

=item C<pg_array_split($string)>

Converts a string representing an array value into a list of values.
N.B. This function is not smart about handling nested braces.    

=cut

sub reset {
    my $self = shift;

=item C<reset>

Reset the connection.

=cut
    if ($self->{INCOPY}) {
	$self->endcopy;
    }
    $self->close_connection;
    $self->connect;
    
}

sub connect {
    # Make the connection using the object's parameters and check results.
    my $self = shift;

    my $dbh = DBI->connect($self->{CONNECTION_STRING},
			   $self->{USERNAME},
			   $self->{PASSWORD},
			   { AutoCommit => 1,
			     RaiseError => 0,
			     PrintError => 1});
    if (not $dbh) {
	my $msg = sprintf("Bad connection for %s: %s\n",
			  $self->{CONNECTION_STRING},
			  $dbh->errstr);
	croak($msg);
    }
    $self->{DBH} = $dbh;
    if ($self->{DRIVER} eq 'SQLite') {
	$self->command("PRAGMA synchronous = OFF");
    }
}



=back

=head1 SEE ALSO

DBI(3), DBD::Pg(3), DBD::SQLite(3)

=head1 AUTHOR

Robert E. Bruccoleri,
C<bruc@acm.org>,
Congenomics, LLC

=cut

__END__
