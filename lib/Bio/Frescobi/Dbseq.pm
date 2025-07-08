# This package provides an interface to the relational database containing
# segmented sequences. It is designed for operations where we wish to
# insert a large number of records into a table where each record
# represents a sequence with associated data.

# There is an inherent bug in this module: if concatenated id's are used
# with uniquification, the actual name field will not be updated properly.

package Bio::Frescobi::Dbseq;
use strict;
use warnings;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use SHA;
use Exporter;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::PGLock;
use Bio::Frescobi::BigQuery;
use Bio::Frescobi::Config;
use Carp qw(verbose croak carp confess);
use Data::Dumper;

@ISA = ('Exporter');
@EXPORT = qw(find_sequence
	     retrieve_sequence_data);

our $AUTOLOAD;

my %component =
    map { ($_ => 1) }
    qw(REPLACE ATTRIBUTE_FOR_NAME DEFAULT_SEQ_TYPE
       WORKDIR SORTMEM DROP_INDEXES UNIQUIFY_NAMES UNIQUIFY_SEPARATOR
       VERBOSE CHECKSEQ BLOCK_SIZE USEDISK DEBUG UNIQUIFY_RANDOM
       UNIQUIFY_RANDOM_LENGTH UNIQUIFY_MAXITER);

my $max_piece_length = 3900;

Bio::Frescobi::PGLock::lock_wait_time(10);
Bio::Frescobi::PGLock::lock_wait_cycles(20000);	# About 3 days.

$VERSION = '0.02';

1;

    
sub new {
    my ($pkg, $class, $self);

    $class = shift;
    eval {($pkg) = caller(0);};
    if ($class ne $pkg) {
	unshift @_, $class;
    }
    $self = {};
    bless $self;
    $self->{PG} = shift;
    &Bio::Frescobi::PGLock::grab_lock($self->{PG},
				      "seq_data");
    $self->{SEQID_HEAD} = "";
    $self->{SEQCNT} = 0;
    $self->{DEFAULT_SEQ_TYPE} = 'nucleic';
    $self->{FIELDS} = undef;
    $self->{DB_NAMES_SEEN} = undef;
    $self->{ATTRIBUTE_FOR_NAME} = "name";
    $self->{DATA_NAMES_SEEN} = {};
    $self->{UNIQUIFY_NAMES} = 0;
    $self->{UNIQUIFY_SEPARATOR} = "_";
    $self->{UNIQUIFY_RANDOM} = 0;
    $self->{UNIQUIFY_RANDOM_LENGTH} = 6;
    $self->{UNIQUIFY_MAXITER} = 1000;
    $self->{SEQ_PIECES} = [];
    $self->{SEQ_DATA} = [];
    $self->{TABLE_DATA} = [];
    $self->{NAMES_TO_DELETE} = [];
    $self->{TABLE_NAME} = undef;
    $self->{NAME_SEEN} = undef;
    $self->{REPLACE} = 0;
    if ($self->{PG}->driver eq 'Pg') {
	$self->{CURTIME} = $self->{PG}->get_single_value("select 'now'::timestamp with time zone");
    }
    else {
	$self->{CURTIME} = $self->{PG}->get_single_value("select datetime('now', 'localtime')");
    }
    $self->{BLOCK_SIZE} = 10000;
    $self->{VERBOSE} = 1;
    $self->{CHECKSEQ} = 1;
    $self->{WORKDIR} = $Bio::Frescobi::Config::tmpdir . "/Dbseq.$$";
    $self->{SORTMEM} = '3G';
    $self->{CHECKSUM_SIZE} = 8;
    $self->{DROP_INDEXES} = 0;
    $self->{STARTED} = 0;
    $self->{FINISHED} = 0;
    $self->{USEDISK} = 1;
    $self->{DEBUG} = 1;
    $self->{DUPLICATED_NAMES} = {};
    $self->{ORIGIN_FOR_NEW_NAMES} = {};
    return $self;
}

sub set_seqid_head {
    my $self = shift;
    my $head = shift;

    croak "seqid_head ($head) must be specified"
	if $head eq "";
    $self->{SEQID_HEAD} = $head;
    $self->{SEQCNT} = 0;
}

sub setup_table {
    my $self = shift;

    my $pg = $self->{PG};
    my ($field, $count, $seqid_seen);

    croak "setup_table can only be called once per session"
	if defined $self->{FIELDS} or defined $self->{TABLE_NAME};
    $self->{TABLE_NAME} = shift;
    $self->{FIELDS} = {};
    $seqid_seen = 0;
    $self->{NAME_SEEN} = 0;
    $count = 0;
    foreach $field ($pg->get_fields_for_table($self->{TABLE_NAME})) {
	$self->{FIELDS}->{$field} = $count++;
	if ($field eq $self->{ATTRIBUTE_FOR_NAME}) {
	    $self->{NAME_SEEN} = 1;
	}
	elsif ($field eq "seqid") {
	    $seqid_seen = 1;
	}
    }
    croak "Seqid field must be specified." if not $seqid_seen;
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self);
    my $pkg = __PACKAGE__;
    my ($pos, $new_val);

    confess "$self is not an object. AUTOLOAD = $AUTOLOAD\n" if not $type;
    confess "$pkg AUTOLOAD function fails on $type\n"
	if $type !~ m/^$pkg$/;
    my $name = uc($AUTOLOAD);
    $name =~ s/^.*:://;
    if (not exists($component{$name})) {
	confess "Error 1: $name is not a valid method for $type.\n";
    }
    elsif (not exists($self->{$name})) {
	confess "Error 2: $name is not a valid method for $type\n";
    }
    elsif (defined $_[0]) {
	$self->{$name} = shift;
    }
    if ($name =~ m/^(REPLACE|UNIQUIFY_NAMES)$/) {
 	&_check_replace_uniquify_conflicts($self);
    }
    if ($name eq 'ATTRIBUTE_FOR_NAME') {
	if (defined $self->{TABLE_NAME}) {
	    $self->{NAME_SEEN} = 0;
	    foreach my $field ($self->{PG}->get_fields_for_table($self->{TABLE_NAME})) {
		if ($field eq $self->{ATTRIBUTE_FOR_NAME}) {
		    $self->{NAME_SEEN} = 1;
		}
	    }
	}
    }
    if ($name eq 'UNIQUIFY_MAXITER') {
	if ($self->{UNIQUIFY_MAXITER} < 10) {
	    confess sprintf("uniquify_maxiter (%s) must be at least 10.\n",
			    $self->{UNIQUIFY_MAXITER});
	}
    }
    if ($name eq 'UNIQUIFY_SEPARATOR') {
	if (length($self->{UNIQUIFY_SEPARATOR}) != 1) {
	    confess sprintf("uniquify_separator (%s) must be a single character.\n",
			    $self->{UNIQUIFY_SEPARATOR});
	}
    }
    return $self->{$name};
}

sub _check_replace_uniquify_conflicts {
    my $self = shift;

    if ($self->{REPLACE} and $self->{UNIQUIFY_NAMES}) {
	carp <<EOF;
The replace and uniquify_names options conflict with each other.
The uniquify_names option will take precedence
EOF
    }
}

sub setup_index {
    my $self = shift;
    print STDERR "setup_index is now obsolete.\n";
}

sub set_index_reconstruction_threshold {
    my $self = shift;

    print STDERR "set_index_reconstruction_threshold is now obsolete.\n";
}

sub begin {
    my $self = shift;

    if (not $self->{USEDISK}) {
	croak "Code is temporarily broken for in memory use.\n";
    }
    
    # Get the sort order right.
    $ENV{LC_ALL} = "C";
    
    croak "begin may only be called once.\n" if $self->{STARTED};
    my $test = `(echo 0old; echo 1new) | sort | head -1`;
    chomp $test;
    croak "Sort order violation: 1new after 0old\n" if $test ne '0old';
    croak "Table fields not specified prior to begin call\n"
	if not defined($self->{FIELDS});
    my $pg = $self->{PG};
    $self->{WORKDIR_CREATED} = 0;
    my $workdir = $self->{WORKDIR};
    if (not -d $workdir) {
	system_with_check("mkdir -p $workdir",
			  $self->{VERBOSE});
	$self->{WORKDIR_CREATED} = 1;
    }
    if (not -d $workdir or
	not -w $workdir) {
	croak "Unable to write working directory $workdir\n";
    }
    set_checksum_size($self->{CHECKSUM_SIZE});
    my $sortmem = $self->{SORTMEM};
    open ($self->{SEQS_FH},
    	  "| sort --stable -T $workdir -S $sortmem >$workdir/seqs.1") ||
    	croak "Unable to open pipe for $workdir/seqs.1: $!\n";
    dated_mesg("Opened $workdir/seqs.1 for write");
    open ($self->{NAMES_FH},
    	  "| sort --stable -T $workdir -S $sortmem >$workdir/names.1") ||
    	croak "Unable to open pipe for $workdir/names.1: $!\n";
    dated_mesg("Opened $workdir/names.1 for write");
    # open ($self->{SEQS_FH}, ">$workdir/seqs.1") ||
    # 	croak "Unable to open pipe for $workdir/seqs.1: $!\n";
    # dated_mesg("Opened $workdir/seqs.1 for write");
    # open ($self->{NAMES_FH}, ">$workdir/names.1") ||
    # 	croak "Unable to open pipe for $workdir/names.1: $!\n";
    # dated_mesg("Opened $workdir/names.1 for write");
    my $query = Bio::Frescobi::BigQuery->new($pg,
					     "select d.checksum, " .
					     "       coalesce(p.seq_piece, ''), " .
					     "       coalesce(p.qual_piece, ''), " .
					     "       d.seqid, " . 
					     "       d.seq_type, " .
					     "       p.segnum " .
					     "  from seq_pieces p," .
					     "       seq_data d " .
					     " where d.seqid = p.seqid " . 
					     " order by p.seqid, p.segnum");
    my $row = $query->next;
    while ($row) {
	my (@seq_pieces, @qual_pieces);
	my ($checksum, $seq_piece, $qual_piece,
	    $seqid, $seq_type, $segnum) = @{$row};
	push (@seq_pieces, $seq_piece);
	push (@qual_pieces, $qual_piece);
	my $prev_segnum = $segnum;
	while ($row = $query->next and
	       $row->[3] eq $seqid) {
	    ($checksum, $seq_piece, $qual_piece, $seqid, $seq_type, $segnum) = @{$row};
	    push (@seq_pieces, $seq_piece);
	    push (@qual_pieces, $qual_piece);
	    if ($prev_segnum + 1 != $segnum) {
		dated_mesg("For $seqid; segnum's out of order at segnum $prev_segnum");
	    }
	    $prev_segnum = $segnum;
	}
	if (not print {$self->{SEQS_FH}} join("\t",
					      $checksum,
					      join("", @seq_pieces),
					      join("", @qual_pieces),
					      '0old',
					      $seq_type,
					      $seqid), "\n") {
	    croak "Print of seqs to sort command failed: $!\n";
	}
    }
    if ($self->{NAME_SEEN}) {
	$query = Bio::Frescobi::BigQuery->new($pg,
					      "select " . $self->{ATTRIBUTE_FOR_NAME} .
					      "  from " . $self->{TABLE_NAME});
	while ($row = $query->next) {
	    my $name = $row->[0];
	    if (not print {$self->{NAMES_FH}} join("\t",
						   $name,
						   '0old'), "\n") {
		croak "Print of names to sort command failed: $!\n";
	    }
	}
    }
    $self->{STARTED} = 1;
}

sub add_data {

    my $self = shift;
    my $pg = $self->{PG};
    my ($key, $val, %fields_seen, @data, $seq_type);
    my ($default_echo, $seq, $qual, $pos, $checksum, $line);
    my ($seqid, $must_insert, $attr_for_name, $name);
    my ($db_seen, $data_seen, $seen, $ok, $new_name, $seen_count);

    if ($self->{VERBOSE}) {
	$default_echo = $pg->default_echo;
	$pg->default_echo(1);
    }

    croak "add_data requires a prior call to set_seqid_head\n"
	if not $self->{SEQID_HEAD};

    croak "add_data required a prior call to begin\n"
	if not $self->{STARTED};

    croak "add_data called after a prior call to finish\n"
	if $self->{FINISHED};

    croak "add_data requires a hash containing table data"
	if scalar(@_) != 1 and ref($_[0]) ne "HASH";
    %fields_seen = ();
    $seq = "";
    $qual = "";
    $seq_type = "";
    $must_insert = 0;
    if ($self->{DROP_INDEXES}) {
        $must_insert = 1;
    }
    $attr_for_name = $self->{ATTRIBUTE_FOR_NAME};
    $name = "";
    foreach $key (keys %{$_[0]}) {
	$val = $_[0]->{$key};
	if ($key eq "seq") {
	    if ($self->{CHECKSEQ}) {
		$seq = _checkseq($val, $_[0]->{$attr_for_name});
	    }
	    else {
		$seq = $val;
	    }
	}
	elsif ($key eq "qual") {
	    $qual = $val;
	}
	elsif ($key eq "seqid") {
	    croak "seqid may not specified to add_data";
	}
	elsif ($key eq "seq_type") {
	    $seq_type = $val;
	}
	else {
	    croak "$key not specified in table."
		if not exists($self->{FIELDS}->{$key});
	    croak "$key specified twice." if exists($fields_seen{$key});
	    $pos = $self->{FIELDS}->{$key};
	    $data[$pos] = $val;
	    if ($key eq $attr_for_name) {
		$name = $val;
		if ($name =~ m/%/) {
		    croak "Name cannot contain a % sign.\n";
		}
		if (not print {$self->{NAMES_FH}} join("\t",
						       $name,
						       '1new'), "\n") {
		    croak "Print of names (2) to sort command failed: $!\n";
		}
	    }
	}
    }
    foreach $val (values %{$self->{FIELDS}}) {
	if (not defined($data[$val])) {
	    $data[$val] = "\\N";
	}
    }
    my $seq_qual = $seq . $qual;
    $checksum = &seq_checksum($seq_qual);
    if ($seq_type eq "") {
	if ($self->{DEFAULT_SEQ_TYPE} eq "") {
	    carp "Input sequence named $name had no sequence type specified.\n";
	    return;
	}
	else {
	    $seq_type = $self->{DEFAULT_SEQ_TYPE};
	}
    }
    if ($self->{USEDISK}) {
	unshift (@data, $checksum, $seq, $qual, '1new', $seq_type);
        if (not print {$self->{SEQS_FH}} join("\t", @data), "\n") {
	    croak "Print of seqs (2) to sort command failed: $!\n";
	}
    }
    else {
	croak "Temporary programming error.\n";
    }
}

sub finish {

    my $self = shift;
    return unless $self->{STARTED};
    return if $self->{FINISHED};
    my $pg = $self->{PG};
    my ($record, $qual, $seq, $seqid, @lines, $result, $line, $values);
    my ($pos, $val, $name, $checksum, $field, $seqid_pos, $seq_type, $status);
    my $attr_for_name = $self->{ATTRIBUTE_FOR_NAME};
    my $name_pos = $self->{FIELDS}->{$attr_for_name};
    $seqid_pos = $self->{FIELDS}->{"seqid"};
    local (*IN, $_);

    $self->_write_data;
    $self->_handle_names;
    $self->{INDEXES} = [ ];
    if ($self->{DROP_INDEXES}) {
	foreach my $table (('seq_pieces', 'seq_data', $self->{TABLE_NAME})) {
	    push (@{$self->{INDEXES}}, $pg->get_indexes($table));
	    $pg->drop_indexes($table);
	}
	print STDERR "In case of error, please execute the following index creation commands:\n";
	print STDERR join("\n", @{$self->{INDEXES}}), "\n";
    }
    my $unique_st;
    if ($self->{USEDISK}) {
        close ($self->{SEQS_FH});
	my $workdir = $self->{WORKDIR};
	open (IN, "<$workdir/seqs.1") ||
	    croak "Unable to open $workdir/seqs.1 for read: $!\n";
	if ($seqid_pos < 0) {
	    croak "finish logic error: pos ($pos) too small.";
	}
	if ($self->{UNIQUIFY_NAMES}) {
	    my $done = 0;
	    my $iter = 0;
	    until ($done) {
		$unique_st = random_string($self->{UNIQUIFY_RANDOM_LENGTH});
		my $name_check =
		    $pg->get_single_value("select $attr_for_name " . 
					  "  from " . $self->{TABLE_NAME} .
					  " where $attr_for_name like '%$unique_st%' " .
					  " limit 1");
		if (not $name_check) {
		    $done = 1;
		}
		else {
		    if (++$iter > $self->{UNIQUIFY_MAXITER}) {
			croak "Unable to find unique string after 1000 tries.\n";
		    }
		}
	    }
	}
	my $prev_checksum = "";
	my $prev_seq_qual = "";
	my $prev_seqid = "";
	while (<IN>) {
	    chomp;
	    my @record = split(/\t/);
	    $checksum = shift(@record);
	    $seq = shift(@record);
	    $qual = shift(@record);
	    $status = shift(@record);
	    $seq_type = shift(@record);
	    my $seq_qual = $seq . $qual;
	    my $ok_to_add = 1;
	    if ($status eq '1new') {
		if ($self->{NAME_SEEN}) {
		    my $name = $record[$name_pos];
		    if (not $self->{UNIQUIFY_NAMES} and
			not $self->{REPLACE} and
			exists($self->{DUPLICATED_NAMES}->{$name})) {
			$ok_to_add = 0;
			print STDERR "Duplicate sequence for $name skipped in new data.\n";
		    }
		    elsif ($self->{REPLACE} and
			   exists($self->{DUPLICATED_NAMES}->{$name})) {
			if ($self->{DUPLICATED_NAMES}->{$name} == 1) {
			    $ok_to_add = 0;
			    print STDERR "Sequence $name skipped because of replace flags and multiple copies in input data.\n";
			}
			else {
			    $self->{DUPLICATED_NAMES}->{$name} = 1;
			}
		    }
		}
	    }
	    if ($checksum eq $prev_checksum and
		$prev_seq_qual eq $seq_qual) {
		$seqid = $prev_seqid;
	    }
	    else {
		$seqid = undef;
		if ($status eq '0old') {
		    $seqid = $record[0];
		}
		elsif (not defined $seqid and $status eq '1new' and $ok_to_add) {
		    $seqid = $self->{SEQID_HEAD} . ++$self->{SEQCNT};
		    push(@{$self->{SEQ_PIECES}},
			 &create_seqpieces($seqid, $seq, $qual));
		    my $line = sprintf("%s\t%s\t%d\t%s\t%s\n",
				       $seqid,
				       $checksum,
				       length($seq),
				       $self->{CURTIME},
				       $seq_type);
		    push(@{$self->{SEQ_DATA}}, $line);
		}
		$prev_checksum = $checksum;
		$prev_seq_qual = $seq_qual;
		$prev_seqid = $seqid;
	    }
	    if ($status eq '1new' and $ok_to_add) {
		if ($seqid_pos >= scalar(@record)) {
		    croak "finish logic error: pos ($pos) too big.";
		}
		croak "seqid ($seqid) specified for record."
		    if $record[$seqid_pos] ne "\\N";
		$record[$seqid_pos] = $seqid;
		if ($self->{NAME_SEEN}) {
		    my $name = $record[$name_pos];
		    if ($self->{UNIQUIFY_NAMES} and
			exists($self->{DUPLICATED_NAMES}->{$name})) {
			my $new_name = $name .
			    $self->{UNIQUIFY_SEPARATOR} .
				$unique_st .
				    $self->{DUPLICATED_NAMES}->{$name}++;
			$record[$name_pos] = $new_name;
			if (not $self->{UNIQUIFY_RANDOM}) {
			    $self->{ORIGIN_FOR_NEW_NAMES}->{$new_name} = $name;
			}
		    }
		}
		push(@{$self->{TABLE_DATA}},
		     join("\t", @record) . "\n");
		if (scalar(@{$self->{TABLE_DATA}}) > $self->{BLOCK_SIZE}) {
		    $self->_write_data;
		}
	    }
	}
	close (IN);
	$self->_write_data;
    }
    foreach my $index_sql (@{$self->{INDEXES}}) {
	$pg->command($index_sql);
    }
    if (not $self->{UNIQUIFY_RANDOM}) {
	$self->_fix_random_names;
    }
    $self->{FINISHED} = 1;
#    print STDERR Dumper($self);
}

sub _handle_names {
    my $self = shift;
    my $pg = $self->{PG};
    my $workdir = $self->{WORKDIR};
    my $attr_for_name = $self->{ATTRIBUTE_FOR_NAME};
    local (*IN, $_);

    return if not $self->{NAME_SEEN};
    close ($self->{NAMES_FH});
    open (IN, "<$workdir/names.1") ||
	croak "Unable to open $workdir/names.1: $!\n";
    my $prev_name = "";
    my $prev_state = "";
    my @names_to_delete;
    while (<IN>) {
	chomp;
	my ($name, $state) = split(/\t/);
	if ($state eq '1new' and $prev_state eq "") {
	    # print STDERR "new $name is OK\n";
	}
	elsif ($state eq '1new' and
	       $prev_state eq '0old' and
	       $prev_name eq $name) {
#	    print STDERR "new $name is new and replicated in old database.\n";
	    if ($self->{REPLACE}) {
		push (@names_to_delete, quotify($name));
		if (scalar(@names_to_delete) > 1000) {
		    $pg->command("delete from " . $self->{TABLE_NAME} .
				 " where $attr_for_name in (" .
				 join(",", @names_to_delete) . ")");
		    @names_to_delete = ();
		}
	    }
	    else {
		$self->{DUPLICATED_NAMES}->{$name} = 1;
	    }
	}
	elsif ($state eq '1new' and
	       $prev_state eq '1new' and
	       $prev_name eq $name) {
#	    print STDERR "new $name is new and replicated in new data.\n";
	    if ($self->{REPLACE}) {
		print STDERR "Only one replacement of $name will be used.\n";
		$self->{DUPLICATED_NAMES}->{$name} = -1;
	    }
	    else {
		$self->{DUPLICATED_NAMES}->{$name} = 1;
	    }
	}
	elsif ($state eq '1new' and $prev_name ne $name) {
	    # print STDERR "new $name is OK also.\n";
	}
	elsif ($state eq '0old' and $prev_name eq "") {
	    # print STDERR "old $name is ok.\n";
	}
	elsif ($state eq '0old' and $prev_name eq $name) {
	    print STDERR "old $name is not unique.\n";
	}
	$prev_state = $state;
	$prev_name = $name;
    }
    if (scalar(@names_to_delete) > 0) {
	$pg->command("delete from " . $self->{TABLE_NAME} .
		     " where $attr_for_name in (" .
		     join(",", @names_to_delete) . ")");
	@names_to_delete = ();
    }
}

sub _fix_random_names {
    my $self = shift;
    my $pg = $self->{PG};
    my $name_attribute = $self->{ATTRIBUTE_FOR_NAME};
    my $table = $self->{TABLE_NAME};
    my $delim = $self->{UNIQUIFY_SEPARATOR};
    my $backslash = "\\";
    foreach my $name (keys %{$self->{ORIGIN_FOR_NEW_NAMES}}) {
	my $original_name = $self->{ORIGIN_FOR_NEW_NAMES}->{$name};
	my @name_candidates;
	if ($pg->driver eq 'Pg') {
	    @name_candidates =
		$pg->get_array_for_field("select $name_attribute " .
					 "  from $table " .
					 " where $name_attribute ~ " . quotify("^" . ${original_name} . ${delim}));
	}
	else {    
	    my $pattern = $original_name . $delim;
	    $pattern =~ s/$delim/${backslash}${delim}/g;
	    @name_candidates =
		$pg->get_array_for_field("select $name_attribute " .
					 "  from $table " .
					 " where like('" . $pattern . "', $name_attribute, '" .  $backslash . "')");
	}
	my $max_n = 1;
	foreach my $cand (@name_candidates) {
	    my ($n) = ($cand =~ m/(\d+)$/);
	    if (defined $n) {
		$max_n = max($n, $max_n);
	    }
	}
	my $new_name = "";
	for (my $i = 1; $i <= $self->{UNIQUIFY_MAXITER}; $i++) {
	    my $test_name = sprintf("%s%s%d", $original_name, $delim, $max_n + $i);
	    my $count = $pg->get_single_value("select count(*) from $table " .
					      " where $name_attribute = " . quotify($test_name));
	    if ($count == 0) {
		$new_name = $test_name;
		last;
	    }
	}
	if ($new_name eq "") {
	    confess "Unable to find unique name for $original_name.\n";
	}
	$pg->command("update $table set $name_attribute = " . quotify($new_name) .
		     " where $name_attribute = " . quotify($name));
    }
}


sub DESTROY {
    my $self = shift;
    my $pg = $self->{PG};
    $self->finish;
    if ($self->{WORKDIR_CREATED}) {
	my $modifier = $self->{DEBUG} ? "echo " : "";
	system_with_check("$modifier rm -rf $self->{WORKDIR}",
			  $self->{VERBOSE});
    }
    &Bio::Frescobi::PGLock::free_lock($pg, "seq_data");
}

sub create_seqpieces {
    # Returns an array of lines that can be passed to the PostgreSQL copy
    # command that represent the sequence in $seq whose quality scores are
    # in $qual and whose seqid is $seqid. A zero length sequence returns a
    # null list.

    my ($seqid, $seq, $qual) = @_;
    my ($segcount, $l, $seq_segment, $qual_segment, @ret);
    my ($len_seq, $len_qual, $seq_cursor, $qual_cursor);

    if ($seq =~ /[\\\t]/ or $qual =~ /[\\\t]/ or $seqid =~ /[\\\t]/) {
        print STDERR "create_seqpieces cannot correctly handle backslashes or tabs in its data. Please reprogram.\n";
        exit 1;
    }
    @ret = ();
    $len_seq = length($seq);
    $len_qual = length($qual);
    if ($len_seq > 0) {
	$segcount = 1;
	$seq_cursor = 0;
	$qual_cursor = 0;
	for (;; $segcount++) {
	    $l = &min($len_seq - $seq_cursor, $max_piece_length);
	    $seq_segment = substr($seq, $seq_cursor, $l);
	    $seq_cursor += $l;
	    if ($seq_segment eq "") {
		$seq_segment = "\\N";
	    }
	    $l = &min($len_qual - $qual_cursor, $max_piece_length);
	    $qual_segment = substr($qual, $qual_cursor, $l);
	    $qual_cursor += $l;
	    if ($qual_segment eq "") {
		$qual_segment = "\\N";
	    }
	    push(@ret,
		 sprintf "%s\t%d\t%s\t%s\n", $seqid, $segcount, $seq_segment, $qual_segment);
	    last if $len_seq == $seq_cursor  and $len_qual == $qual_cursor;
	}
    }
    return @ret;
}

sub retrieve_sequence_data {
    my ($pg, $seqid) = @_;
    my ($data, $seq, $qual);
    my (@ret);

    $seq = "";
    $qual = "";
    $data = $pg->query("select seq_piece, qual_piece, seqid, segnum from seq_pieces" .
		       "       where seqid = '$seqid' order by segnum",
		       0);
    my (@seq_pieces, @qual_pieces);
    foreach my $rowp (@{$data}) {
	push (@seq_pieces, $rowp->[0]);
	if (defined $rowp->[1]) {
	    push (@qual_pieces, $rowp->[1]);
	}
    }
    $seq = join("", @seq_pieces);
    $qual = join("", @qual_pieces);
    @ret = ($seq, $qual);
    return @ret;
}

sub find_sequence {
    my ($pg, $seq, $qual, $checksum, $verbose) = @_;
    my ($saved_seq, $saved_qual, $seqid, $i);

    if (not defined $checksum or $checksum eq "") {
	$checksum = &seq_checksum($seq . $qual);
    }
    if (not defined $verbose or $verbose eq "") {
	$verbose = 0;
    }

    my $data_ref = $pg->query("select seqid from seq_data where checksum = '$checksum'");
    foreach my $row (@{$data_ref}) {
	$seqid = $row->[0];
	($saved_seq, $saved_qual) = retrieve_sequence_data($pg, $seqid);
	if (uc($saved_seq) eq uc($seq) and $saved_qual eq $qual) {
	    return $seqid;
	}
    }
    print STDERR "Sequence $seq and quality $qual do not match stored sequences with checksum $checksum.\n";
    return undef;
}


sub _checkseq {
    my ($seq, $name) = @_;
    my ($msg, $st);

    if ($seq =~ m/[^A-Za-z0-9*]/) {
	$msg = "Sequence named $name has non alphanumeric characters: ";
	($st = $seq) =~ s/[A-Za-z0-9*]+//g;
	$msg .= $st;
	carp "$msg\n";
	$seq =~ s/[^A-Za-z0-9*]+//g;
    }
    return $seq;
}

sub _look_up_sequence {
    my $self = shift;
    my $pg = $self->{PG};
    my ($seq, $qual, $checksum) = @_;
    my $seq_qual = $seq . $qual;

    my $seqid = &find_sequence($pg, $seq, $qual, $checksum, $self->{VERBOSE});
    return $seqid;
}

sub _write_data {
    my $self = shift;
    my $pg = $self->{PG};
    my ($name, $result);

    if (scalar(@{$self->{SEQ_PIECES}}) > 0) {
	$pg->copy_into_postgres("seq_pieces",
				$self->{SEQ_PIECES});
	$self->{SEQ_PIECES} = [];
    }
    if (scalar(@{$self->{SEQ_DATA}}) > 0) {
	$pg->copy_into_postgres("seq_data",
				$self->{SEQ_DATA});
	$self->{SEQ_DATA} = [];
    }
    if (scalar(@{$self->{TABLE_DATA}}) > 0) {
	$pg->copy_into_postgres($self->{TABLE_NAME},
				$self->{TABLE_DATA});
	$self->{TABLE_DATA} = [];
    }
}

__END__

=head1 NAME

Bio::Frescobi::Dbseq -- class for managing sequence loading into FRESCOBI

=head1 SYNOPSIS

 use Bio::Frescobi::Dbseq;

=head2 Constructor

 $seqobj = Bio::Frescobi::Dbseq->new($pg);

=head2 Object Methods

 $seqobj->set_seqid_head($head);
 $seqobj->setup_table($table);
 $seqobj->replace_mode[($mode)];
 $seqobj->attribute_for_name[($attribute)];
 $seqobj->default_seq_type[($seq_type)];
 $seqobj->uniquify_names_mode[($mode]);
 $seqobj->uniquify_separator[($separator)];
 $seqobj->verbose_mode[($verbosity)];
 $seqobj->checkseq_mode[($mode)];
 $seqobj->block_size[($size)];
 $seqobj->add_data({FIELD => $value, ...});
 $seqobj->finish;

=head2 Convenience Functions

 $seqid = &find_sequence($pg, $sequence, $quality[, $checksum[, $verbose]]);

 ($sequence, $quality) = &retrieve_sequence_data($pg, $seqid);

=head1 DESCRIPTION

The C<Bio::Frescobi::Dbseq> class implements the sequence loading process for
a FRESCOBI database. FRESCOBI provides a sequence repository, where
each unique sequence is stored only once, and the system provides a surrogate
field, the C<SEQID>, to stand in the place of the sequence in other tables
in the database. This class, C<Bio::Frescobi::Dbseq>, is used for implementing
the sequence repository loading process as well as some simple retrieval functions.

The schema also handles quality scores using the PHRED encoding,
namely each quality score is stored as a string of numbers in
text form. Each number can be in the range of 0 to 99 and is
separated from the others by single spaces. If a quality score is
provided with the sequence, then the C<SEQID> is set based on the
concatenation of the sequence and quality score strings. Thus,
two identical sequences with different quality scores are treated
as different sequences and have different C<SEQID>'s.

A checksum is used on the sequence and quality scores to
optimize loading and access. Unfortunately, the checksum is
not sufficient long to assure that difference sequences will
always have different checksums, but it does guarantee that
different checksums imply different sequences.

This class works with three tables to load sequences. Two of the
tables, C<seq_data> and C<seq_pieces>, are essential parts of the
FRESCOBI schema. The third table is specified by the programmer, and
must have a C<SEQID> field within it which will be automatically set
by this class. This third table will be referred to as the high level
table.  The C<seq_data> table contains the C<SEQID>, checksum, length,
creation date, and sequence type (nucleic or protein). The
C<seq_pieces> table contains the segments of the sequences since long
sequences are broken into 3900 character pieces.

When this class is used, the programmer must first create the object.
Then, the various modes must be set to the programmer's need.
Next, the high level table to hold the new C<SEQID>'s must be
specified.  Then, sequence records can be added using the C<add_data>
method.  Each C<add_data> call must include a hash hey (parameter)
named "C<seq>" which holds the actual sequence. The quality is
optionally specified using a hash key named "C<qual>".  Finally, the
programmer must destroy the object by setting the variable which holds
it to C<undef> so that the data is all written out.

This class also provides mechanisms for handling duplicate keys in the
high level table. By default, C<Bio::Frescobi::Dbseq> sets the field name for
the high level table to be C<name>. Using the C<attribute_for_name>
method, the programmer can set the key to any other field in the high
level table.  If a duplicate key is passed in new records, there are
three possible outcomes. By default, duplicate records will be
ignored. If C<replace_mode> is set, then new records will replace
existing records, except that the previous sequences will be still be
left in the database in the C<seq_data> and C<seq_pieces> tables. If
C<uniquify_names_mode> is set, then a new key will be automatically
generated by appending the C<uniquify_separator> (which defaults to
"_") followed by a number starting with 2 and increasing until a
unique name is found.

When data is loaded using this class, substantial data from the
database must loaded into RAM which is used throughout the loading
process. In particular, the complete list of checksums and high level
table keys must be retrieved. In addition, the method accumulates records
with checksum conflicts in RAM and loads them at the end of the
process. Therefore, only one process can load
data at any one time. Therefore, this class will set a lock using the
C<Bio::Frescobi::PGLock> class on the C<seq_data> table.

In order to make the loading process fast, the class uses three
techniques to improve performance.  The first technique is the
creation of unique sequence identifiers. The sequence identifiers have
two serial numbers in them.  The first serial number is generated from
incrementing a sequence within the PostgreSQL database and incremented
once per class invocation.  The second serial number is generated from
a counter within the class.  The advantage of this approach is that
only one database transaction is needed per class invocation to get a
large number of unique identifiers. The <Bio::Frescobi::Dbseq> class does not handle
getting the first part of the serial number -- the programmer must do this through
a call to C<set_seqid_head>.

The second technique is the separation of the easy loads from the hard
ones. During each call of C<add_data>, which is the primary data
loading method, the checksum of each input sequence is compared
against all existing sequences. If not found, the data associated with
the sequence is copied into an array of new records and periodically
loaded into the database.  The periodic load is done with the
PostgreSQL C<copy> command, which is very fast.  If the checksum of a
sequence has already been seen, the record is appended to an array of
records to be added at the end of the loading process.

When the first phase of loading is complete, the class then handles
those sequences which may already be present in the sequence
repository. First, all old records which are being replaced are deleted now.
Then,
the sequence associated with each record is looked up in
the sequence repository. If found, then the high level record gets the
existing C<SEQID>. If not, a new sequence identifier is created.

The third optimization comes from the final loading process.
Every additional sequence requires the lookup of sequences in the database.

=cut

# ' For emacs

=head1 CONSTRUCTOR

The constructor for this object is C<new>, and it expects one argument,
a PostgreSQL object created by C<Bio::Frescobi::CgPg>, i.e.

 $seqobj = Bio::Frescobi::Dbseq->new($pg);

=head1 METHODS

The following methods are provided:

=over 4

=cut

=item C<set_seqid_head($head)>

Specifies the leading string for new C<SEQID>'s created by this
object.  This string must be set to a unique value with respect to the
prefixes of other C<SEQID>'s in the database. Typically, a sequence is
defined in the database, and the C<SEQID> head will be a letter
followed by the next value in the sequence followed by an underscore.
Currently, no check is made for this leading string, so make sure you
provide one.

=cut

=item C<setup_table($table)>

Setup the field information for the high level table. The table must contain
a seqid field, and can contain a name field which can be automatically
checked for duplications. The name field is normally named "C<name>", but
can be renamed by the C<attribute_for_name> method.

=cut

=item C<replace_mode[($mode)]>

Retrieve and optionaly set the replace mode. Replace mode means that
records that have matching names will be replaced, except that
duplicates within the data stream will be ignored. If an argument is
specified, it will be used to set the replace mode. In all cases, the
current value of the replace mode is returned. Replace mode and
uniquify mode cannot be set both at the same time.

=cut

=item C<attribute_for_name[($field)]>

Retrieve and optionally set the field name which is the key for the
high level table.

=cut

=item C<default_seq_type[($field)]>

Retrieve and optionally set the default sequence type for the sequences to come.
Currently, only two values can be used, "C<nucleic>" and "C<protein>".

=cut

=item C<workdir[($size)]>

Retrieve and optionally set the working directory for temporary storage.

=cut

=item C<sortmem[($size)]>

Retrieve and optionally set the Unix sort command memory use option..

=cut

=item C<drop_indexes[($mode)]>

Retrieve and optionally set the drop indexes mode.

=cut

=item C<uniquify_names_mode[($mode)]>

Retrieve and optionally set the uniquify names mode. If uniquify names mode
is set, then any duplicate names entries will have new names set by appending
successive integers to the name. This option conflicts with replace mode.

=cut

=item C<uniquify_separator[($separator)]>

Retrieve and optionally set the separator for uniquifying names
in the high level table.

=cut

=item C<verbose_mode[($mode)]>

Retrieve and optionally set the verbose mode, which is either true or false (non-zero or zero).
This controls how verbiage is output when the class operates.

=cut

=item C<checkseq_mode[($mode)]>

Retrieve and optionally set the checkseq mode, which is either true or
false (non-zero or zero).  This mode controls whether the class checks
for non-alphanumeric characters in sequences and whether it will
delete them when found, along with a error message.

=cut

=item C<block_size[($size)]>

Retrieve and optionally set the number of records loaded per batch by
C<add_data>,

=cut

=item C<< add_data({ field_name => value, ...}) >>

 The C<add_data> method is the primary method for adding new records to
 the database. The input parameter is a single hash reference
 containing the various fields in a single record of the high level
 table and a sequence and optional quality score.  The keys in the hash
 reference correspond to the field names in the high level table, and
 the values in the hash reference are the new values for each field.
 Two keys in the input hash references are special, C<seq> and C<qual>,
 which contain the sequence and quality score for the record.

 This method implements the loading optimizations described in the
 DESCRIPTION section above.

=cut

=item C<< finish >>

Complete the loading of data for this class.

=cut

=back

=head1 SEE ALSO

Pg(3).

=head1 AUTHOR

Robert E. Bruccoleri,
C<bruc@acm.org>,
Congenomics, LLC

=cut
