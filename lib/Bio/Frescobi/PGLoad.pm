# This package provides a simple interface to the PostgreSQL copy command

package Bio::Frescobi::PGLoad;
use strict;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use Exporter;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::CgPg;
use Carp;
@ISA = ('Exporter');
@EXPORT = qw();

$VERSION = '0.01';

1;

=head1 NAME

Bio::Frescobi::PGLoad -- class for loading data using the PostgreSQL C<copy> command. 

=head1 SYNOPSIS

 use Bio::Frescobi::PGLoad;

=head2 Constructor

 $loader = Bio::Frescobi::PGLoad->new($pg, $table);

=head2 Object Methods

 $loader->swap_table_mode[($mode)]
 $loader->delimiter[($delimiter)]
 $loader->verbose[($mode)]
 $loader->blocksize[($size)]
 $loader->truncate[($flag)]
 $loader->die_on_unknown_field[($flag)]
 $loader->drop_indexes[($flag)]
 $loader->length_limit[($limit)]
 $loader->output_function[($function_reference)]
 $loader->add_data({ "field" => value, ... }) or
 $loader->add_data(value, ...);
 $loader->record_count
 $loader->finish
 $loader->write_data

=head1 DESCRIPTION


=cut


sub new {
    my ($pkg, $class, $self, $pg);

    $class = shift;
    eval {($pkg) = caller(0);};
    if ($class ne $pkg) {
	unshift @_, $class;
    }
    $self = {};
    bless $self;
    $pg = $self->{PG} = shift;
    $self->{TABLE} = shift;
    croak "A table name must be specified when creating a $pkg object.\n"
	unless $self->{TABLE};
    $self->{FIELDS} = {};
    my $ind = 0;
    foreach my $field ($pg->get_fields_for_table($self->{TABLE})) {
	$self->{FIELDS}->{$field} = $ind++;
    }
    croak "A valid table name must be specified when creating a $pkg object.\n"
	unless scalar(keys %{$self->{FIELDS}}) > 0;
    @{$self->{INDEXING_COMMANDS}} = $pg->get_indexes($self->{TABLE});
    $self->{BLOCKSIZE} = 10000;
    $self->{SWAP_TABLE_MODE} = 0;
    $self->{DELIMITER} = "\t";
    $self->{OUTPUT_FUNCTION} = undef;
    $self->{DATA} = [];
    $self->{VERBOSE} = 1;
    $self->{STARTED} = 0;
    $self->{ENDED} = 0;
    $self->{TRUNCATE} = 0;
    $self->{DROP_INDEXES} = 0;
    $self->{STARTED} = 0;
    $self->{DIE_ON_UNKNOWN_FIELD} = 0;
    $self->{RECORD_COUNT} = 0;
    if ($pg->version ge "7.3") {
	$self->{LENGTH_LIMIT} = 500000000;
    }
    elsif ($pg->version ge "7.1") {
	$self->{LENGTH_LIMIT} = 100000000;
    }
    else {
	$self->{LENGTH_LIMIT} = 8000;
    }
    
    return $self;
}

sub swap_table_mode {
    my $self = shift;
    if (defined $_[0]) {
	if ($self->{STARTED}) {
	    croak "Swap table mode cannot be changed once loading has begun\n";
	}
	else {
	    $self->{SWAP_TABLE_MODE} = $_[0];
	}
    }
    return $self->{SWAP_TABLE_MODE};
}

sub delimiter {
    my $self = shift;
    if (defined $_[0]) {
	if ($self->{STARTED}) {
	    carp "The delimiter cannot be changed once loading has begun\n";
	}
	else {
	    $self->{DELIMITER} = $_[0];
	}
    }
    croak "delimiter can only be a tab character in this release.\n"
	if $self->{DELIMITER} ne "\t";
    return $self->{DELIMITER};
}

sub verbose {
    my $self = shift;
    if (defined $_[0]) {
	$self->{VERBOSE} = $_[0];
    }
    return $self->{VERBOSE};
}

sub blocksize {
    my $self = shift;
    if (defined $_[0]) {
	$self->{BLOCKSIZE} = $_[0];
    }
    return $self->{BLOCKSIZE};
}

sub truncate {
    my $self = shift;
    if (defined $_[0]) {
	$self->{TRUNCATE} = $_[0];
    }
    return $self->{TRUNCATE};
}

sub die_on_unknown_field {
    my $self = shift;
    if (defined $_[0]) {
	$self->{DIE_ON_UNKNOWN_FIELD} = $_[0];
    }
    return $self->{DIE_ON_UNKNOWN_FIELD};
}

sub drop_indexes {
    my $self = shift;
    if (defined $_[0]) {
	$self->{DROP_INDEXES} = $_[0];
    }
    return $self->{DROP_INDEXES};
}

sub length_limit {
    my $self = shift;
    if (defined $_[0]) {
	$self->{LENGTH_LIMIT} = $_[0];
    }
    return $self->{LENGTH_LIMIT};
}

sub output_function {
    my $self = shift;
    if (defined $_[0]) {
	$self->{OUTPUT_FUNCTION} = $_[0];
    }
    return $self->{OUTPUT_FUNCTION};
}

sub record_count {
    my $self = shift;
    return $self->{RECORD_COUNT};
}

sub add_data {
    my $self = shift;
    my $pg = $self->{PG};
    my $delim = $self->{DELIMITER};
    
    if (not $self->{STARTED}) {
	$self->_initialize_loading;
    }
    if ($self->{ENDED}) {
	croak "Data added after data load completed.\n";
    }
    my $arg1 = shift;
    if (ref($arg1) eq "HASH") {
	$self->_add_data_by_hash($arg1);
    }
    else {
	unshift(@_, $arg1);
	my @data;
	my $val;
	foreach my $field (@_) {
	    if (not defined $field) {
		push(@data, '\N');
	    }
	    else {
		$val = $field;
		$val =~ s/\\/\\\\/g;
		$val =~ s/\n/\\n/g;
		$val =~ s/\r/\\r/g;
		$val =~ s/$delim/\\$delim/g;
 		$val =~ s/\x00/\\x00/g;
		push(@data, $val);
	    }
	}
	my $line = join($delim, @data) . "\n";
	$self->_add_line($line);
    }
}

sub _initialize_loading {
    my $self = shift;
    if ($self->{SWAP_TABLE_MODE}) {
	$self->{NEW_TABLE} = "new_" . $self->{TABLE};
	my $pg = $self->{PG};
	$pg->command("select * into " . $self->{NEW_TABLE} .
		     "  from " . $self->{TABLE} .
		     " limit 0");
	$self->{LOAD_TABLE} = $self->{NEW_TABLE};
    }
    else {
	if ($self->{DROP_INDEXES}) {
	    $self->_drop_indexes;
	}
	if ($self->{TRUNCATE}) {
	    $self->{PG}->command("truncate " . $self->{TABLE});
	}
	$self->{LOAD_TABLE} = $self->{TABLE};
    }
    $self->{STARTED} = 1;
}

sub _add_data_by_hash {
    my $self = shift;
    my $pg = $self->{PG};
    my $hash = shift;
    my ($field, $line, $delim, $msg, @data, $i);

    croak "add_data_by_hash requires a hash reference.\n"
	unless ref($hash) eq "HASH";
    $delim = $self->{DELIMITER};
    $msg = "";
    for $field (keys %{$hash}) {
	if (not exists($self->{FIELDS}->{$field})) {
	    $msg .= "$field ";
	}
    }
    if ($msg ne "") {
	if ($self->{DIE_ON_UNKNOWN_FIELD}) {
	    croak "$msg were not fields in " . $self->{TABLE} . "\n";
	}
	else {
	    carp "$msg were not fields in " . $self->{TABLE} . "\n";
	    return;
	}
    }
    $line = "";
    @data = ();
    for $field (keys %{$self->{FIELDS}}) {
	my $val;
	my $pos = $self->{FIELDS}->{$field};
	if (not exists($hash->{$field}) or
	    not defined($hash->{$field})) {
	    $val = '\N';
	}
	else {
	    $val = $hash->{$field};
	    $val =~ s/\\/\\\\/g;
	    $val =~ s/\n/\\n/g;
	    $val =~ s/\r/\\r/g;
	    $val =~ s/$delim/\\$delim/g;
	}
	$data[$pos] = $val;
    }
    $line = join($self->{DELIMITER}, @data) . "\n";
    $self->_add_line($line);
}

sub _add_line {
    my $self = shift;
    my $line = shift;
    
    if (length($line) > $self->{LENGTH_LIMIT}) {
	$line =~ s/\t/\n/g;
	carp sprintf("add_data has an oversize record (%d). Fields:\n%s",
		     length($line),
		     $line);
    }
    else {
	push (@{$self->{DATA}}, $line);
	$self->{RECORD_COUNT} += 1;
    }
    if (scalar(@{$self->{DATA}}) > $self->{BLOCKSIZE}) {
	$self->write_data;
    }
}

sub finish {
    my $self = shift;
    my $pg = $self->{PG};

    return if $self->{ENDED};
    if (scalar(@{$self->{DATA}}) > 0) {
	$self->write_data;
    }
    if ($self->{SWAP_TABLE_MODE}) {
	my @index_names;
	my $table = $self->{TABLE};
	my $new_table = $self->{NEW_TABLE};
	my $old_table = "old_${table}";
	foreach my $index_def (@{$self->{INDEXING_COMMANDS}}) {
	    my $cmd = lc($index_def);
	    push (@index_names, (split(/\s+/, $cmd))[2]);
	    $cmd =~ s/create index /create index new_/;
	    $cmd =~ s/on $table/on $new_table/;
	    $pg->command($cmd);
	}
	$pg->command("alter table $table rename to $old_table");

	foreach my $index_name (@index_names) {
	    $pg->command("alter index ${index_name} rename to old_${index_name}");
	}
	$pg->command("alter table $new_table rename to $table");

	foreach my $index_name (@index_names) {
	    $pg->command("alter index new_${index_name} rename to ${index_name}");
	}
	$pg->command("analyze $table");
	$pg->command("drop table $old_table");
    }
    else {
	if ($self->{DROP_INDEXES}) {
	    $self->_create_indexes;
	}
    }
    $self->{ENDED} = 1;
}

sub DESTROY {
    my $self = shift;
    $self->finish;
}

sub write_data {
    my $self = shift;
    my $pg = $self->{PG};
    my $default_echo = $pg->default_echo;

    if ($self->{VERBOSE}) {
	$pg->default_echo(1);
    }
    else {
	$pg->default_echo(0);
    }
    if ($self->{OUTPUT_FUNCTION}) {
	&{$self->{OUTPUT_FUNCTION}}($self->{LOAD_TABLE}, $self->{DATA});
    }
    else {
	$pg->copy_into_postgres($self->{LOAD_TABLE},
				$self->{DATA});
    }
    $self->{DATA} = [];
    $pg->default_echo($default_echo);
}

sub _drop_indexes {
    my $self = shift;
    my $pg = $self->{PG};

    dated_mesg("Indexes are being deleted. Creation commands are:");
    foreach my $cmd (@{$self->{INDEXING_COMMANDS}}) {
	print STDERR "  $cmd\n";
    }
    $pg->drop_indexes($self->{TABLE});
}

sub _create_indexes {
    my $self = shift;
    my $pg = $self->{PG};

    foreach my $command (@{$self->{INDEXING_COMMANDS}}) {
	$pg->command($command);
    }
}

    
