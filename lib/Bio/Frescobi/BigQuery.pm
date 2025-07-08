package Bio::Frescobi::BigQuery;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Carp;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Genutil;

@ISA = ('Exporter');
@EXPORT = qw();

$VERSION = '0.01';

my $cursor_count = 0;

=head1 NAME

Bio::Frescobi::BigQuery -- class for getting results from a big PostgreSQL query

=head1 SYNOPSIS

 use Bio::Frescobi::BigQuery;

=head2 Constructor

 $query = Bio::Frescobi::BigQuery->new($pg, $sql[, $fetch_size]);

=head2 Object Methods

 $query->next;
 $query->verbose_mode([$new_mode]);

=head1 DESCRIPTION

The C<Bio::Frescobi::BigQuery> class implements a simple method for handling results
from a large query in a memory efficient manner. The need for this class arose from
the fact that PostgreSQL will slurp into memory all the results from
a SELECT statement. This class uses an SQL cursor to define a scanning mechanism for
getting the results. PostgreSQL will put the query results into temporary disk files.
    
=head1 CONSTRUCTOR

The constructor for this object is C<new>, and it expects two or three arguments,
a connection to PostgreSQL (from the L<Bio::Frescobi::CgPg> method), an SQL statement,
and an optional fetch size. The fetch size defaults to 10000 if not specified.
a PostgreSQL object created by C<Bio::Frescobi::CgPg>, i.e.

 $query = Bio::Frescobi::BigQuery->new($pg, $sql[, $fetch_size]);

The constructor can take a long time to execute because the method executes
the SQL statement and then issues the first FETCH. If the query is complex,
most of the work is done at this point.

=head1 METHODS

The following methods are provided:

=over 4

=cut


sub new {
    my ($pkg, $class, $self);
    my ($i);

    $class = shift;
    eval {($pkg) = caller(0);};
    if ($class ne $pkg) {
	unshift @_, $class;
    }
    $self = {};
    bless $self;
    
    my ($pg, $sql, $fetch_size) = @_;
    
    if (not defined $sql) {
	croak "SQL command must be defined when creating a new " . ref($self) . "\n";
    }
    $fetch_size = 10000 if not defined $fetch_size;
    $self->{PG} = $pg;
    $self->{SQL} = $sql;
    $self->{FETCH_SIZE} = $fetch_size;
    $cursor_count += 1;
    $self->{CURSOR} = "bigquery_${cursor_count}_$$";
    # $pg->command("BEGIN");
    # $pg->command("DECLARE $self->{CURSOR} CURSOR for $sql");
    my $sth = $pg->{DBH}->prepare($sql);
    $self->{STH} = $sth;
    my $rv = $self->{STH}->execute;
    if (not $rv) {
	my $msg = sprintf("Command %s returned status %d (%s)\n",
			  $sql,
			  $sth->err,
			  $sth->errstr);
	$pg->{DBH}->do("rollback transaction", { RaiseError => $self->{DIE_ON_ERROR},
						 PrintError => 1});
	croak($msg);
    }
    $self->{DONE} = 0;
    $self->{IREC} = -1;
    $self->{VERBOSE} = 0;
    $self->{DATA} = [ $self->{STH}->fetchrow_arrayref ];
    $self->{NRECORDS} = scalar(@{$self->{DATA}});
    return $self;
}

sub verbose_mode {

=item C<verbose_mode[($mode)]>

Return and optionally set the verbose mode to C<$mode>.
Verbose mode will generate more output during command execution.
Note, the debugging display of SQL commands is controlled by settings
on the PostgreSQL connection object originally passed to the constructor
of this object.

=cut

    my $self = shift;
    if (defined($_[0])) {
	$self->{VERBOSE} = $_[0];
    }
    return $self->{VERBOSE};
}

sub next {

=item C<next>

Return a reference to an array containing the next row of data
from the SQL command. If there is no more data available, then this method
will return the undefined value.

=cut

    my $self = shift;
    my $pg = $self->{PG};
    my (@rows);

    my $verbose = $self->{VERBOSE};
    if ($self->{DONE}) {
	dated_mesg("BigQuery called already done.") if $verbose;
	return undef;
    }
    $self->{IREC} += 1;
    if ($self->{IREC} >= $self->{NRECORDS}) {
#	dated_mesg(sprintf("IREC (%d) >= ntuples (%d)", $self->{IREC}, $self->{NRECORDS})) if $verbose;
	if ($self->{NRECORDS} == 0) {
	    $self->{DONE} = 1;
	    # $pg->command("END");
	    return undef;
	}
	$self->{DATA} = [ $self->{STH}->fetchrow_arrayref ];
	$self->{NRECORDS} = scalar(@{$self->{DATA}});
#	dated_mesg(sprintf("Fetch yielded %d tuples", $self->{NRECORDS})) if $verbose;
	if ($self->{NRECORDS} == 0) {
	    $self->{DONE} = 1;
	    # $pg->command("END");
	    return undef;
	}
	$self->{IREC} = 0;
    }
#    dated_mesg(sprintf("Record %d: %s", $self->{IREC}, join("\t", @{$self->{DATA}->[$self->{IREC}]}))) if $verbose;
    return $self->{DATA}->[$self->{IREC}];
}

sub DESTROY {
    my $self = shift;
    my $pg = $self->{PG};

    dated_mesg("DESTROY called") if $self->{VERBOSE};
    # $pg->command("END") if not $self->{DONE};
    $self->{DATA} = undef;
}
    
1;
