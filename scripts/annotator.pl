#!/usr/bin/perl
# Perl script to annotate a set of contigs
#
# Copyright (c) 2009 Robert E. Bruccoleri
# Licensed under the GNU Public License version 2

use warnings;
use strict;

my $usage = <<EOF;
annotator.pl -dbname=name
             [-driver=[Pg|SQLite] ]
             [-db=db]         
             [-method=string]
             [-extraopts=<string>]
             [-run=integer]
             [-library=string]
             [-tmp=dir]
             [-maxproc=integer]
             [-maxhits=integer]
             [-cutoff=float]
             [-idcutoff=float]
             [-replace]
             [-timeout=integer]
             [-genome_size=integer]
             [-sql=string]
             [-[no]filter]
             [-minscore=integer]
             [-gapopen=integer]
             [-block=integer]  Default: 1
             [-[no]debug]
             [-minlength=integer]
             [-wordsize=integer]
             [-[no]selfself]
             [-nice=string]
             [-fasty=string]
             [-ggsearch=string]
             [-exonerate=string]

    Default driver is "Pg", i.e. Postgres.
    Default db is nr.
    Default method is blastn
    Either run, library, sql can be specified, but not more than one.
    Default run is blank which means all sequences.
    Default tmp is /tmp/annotator.<pid>
    Default maxproc = 1
    Default cutoff is 0.1
    Default idcutoff is 0.0
    dbname must be specified.
    Default maxhits is 1000
    Default timeout is 31415926 seconds (about a year)
    -replace means replacing any existing annotations.
    Default genome size is 2000000000
    Default filter value is blastall default.
    Default minscore is undefined.
    Default gapopen is undefined.
    Default block is 1
    Default debug is false.
    Default minlength is 6. Lower will not work with current version of blastall.
    Default wordsize is up to blastall
    Default nice (command) is "nice -19"
    Default fasty (program) is "fasty35"
    Default ggsearch (program) is "ggsearch35"
    Default exonerate (program) is "exonerate";


EOF
    ;

use Carp;
use Getopt::Long;
use POSIX ();
use File::Basename;
use Bio::SearchIO;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Dbseq;
use Bio::Frescobi::PGLoad;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::Frescobi;
use Bio::Frescobi::Config;
use Bio::Frescobi::PGLock;

my ($run_loader, $hsps_loader, $hits_loader);
my ($internal);
my $hit_count = 0;
my @seqids_to_annotate;

my $dbname = "";
my $common_db_name = "nr";
my $method = 'blastn';
my $extra_options = "";
my $run_number = -1;
my $library = "";
my $workdir = "";
my $maxproc = 1;
my $cutoff = 0.1;
my $idcutoff = 0.0;
my $maxhits = 1000;
my $replace = 0;
my $sql_select = "";
my $timeout = 31415926;
my $genome_size = 2000000000;
my $filter = undef;
my $minscore = undef;
my $gapopen = undef;
my $block = 1;
my $debug = 0;
my $minlength = 6;
my $selfself = 0;
my $wordsize = -1;
my $driver = "Pg";
my $nice_command = "nice -19";
my $fasty_program = "fasty35";
my $ggsearch_program = "ggsearch35";
my $exonerate_program = "exonerate";

print STDERR join(" ", @ARGV), "\n";
GetOptions("dbname=s" => \$dbname,
	   "driver=s" => \$driver,
	   "db=s" => \$common_db_name,
	   "method=s" => \$method,
	   "extraopts=s" => \$extra_options,
	   "run=i" => \$run_number,
	   "library=s" => \$library,
	   "tmp=s" => \$workdir,
	   "maxproc=i" => \$maxproc,
	   "maxhits=i" => \$maxhits,
	   "cutoff=f" => \$cutoff,
	   "idcutoff=f" => \$idcutoff,
	   "replace!" => \$replace,
	   "timeout=i" => \$timeout,
	   "genome_size=i" => \$genome_size,
	   "sql=s" => \$sql_select,
	   "filter!" => \$filter,
	   "minscore=i" => \$minscore,
	   "gapopen=i" => \$gapopen,
	   "block=i" => \$block,
	   "debug!" => \$debug,
	   "minlength=i" => \$minlength,
	   "wordsize=i" => \$wordsize,
	   "selfself!" => \$selfself,
	   "nice=s" => \$nice_command,
	   "fasty=s" => \$fasty_program,
	   "ggsearch=s" => \$ggsearch_program,
	   "exonerate=s" => \$exonerate_program);

$method = lc($method);

if ($idcutoff > 1.0 or $idcutoff < 0.0) {
    die "idcutoff switch ($idcutoff) must be between 0 and 1\n";
}

if ($dbname eq "") {
    die "The database name (-dbname option) must be specified.\n";
}

if ($minlength < 6) {
    print STDERR "minlength raised to 6.\n";
    $minlength = 6;
}

if ($block > 1 and $method =~ m/exonerate/) {
    print STDERR "Blocking of sequences for exonerate is not supported. Blocking set to 1.\n";
    $block = 1;
}


if ($workdir eq "") {
    $workdir = $tmpdir . "/annotator.$$";
}

if ($workdir !~ /^\//) {
    die "Temporary directory ($workdir) must begin with a /\n";
}

my $delete_workdir = 0;
if (! -d $workdir) {
    system_with_check("mkdir -p $workdir");
    system_with_check("chmod 755 $workdir");
#    $delete_workdir = 1;
}

if (scalar(@ARGV) != 0) {
    print STDERR $usage;
    exit 1;
}
if ($method =~ /^t?blast[xpn]$/) {
    if (-r "$workdir/blastall") {
	system_with_check("rm -f $workdir/blastall", 1);
    }
    my $blastall_program = `which blastall`;
    chomp $blastall_program;
    system_with_check("ln -s $blastall_program $workdir/blastall", 1);
}
elsif ($method =~ /^blast\+t?blast[xpn]$/) {
    if (-r "$workdir/blast+bin") {
	system_with_check("rm -f $workdir/blast+bin", 1);
    }
    my $blastplus_makeblastdb = `which makeblastdb`;
    chomp $blastplus_makeblastdb;
    my $blastplus_bindir = dirname($blastplus_makeblastdb);
    system_with_check("ln -s $blastplus_bindir $workdir/blast+bin", 1);
}
    
$ENV{WORKDIR} = $workdir;

$SIG{ALRM} = sub { die "timeout"; };
set_line_buffering;
my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => 0);
$pg->default_echo(1);
if ($driver eq 'Pg') {
    $pg->command("set enable_seqscan = 'off'");
}

my $spec_count = 0;
$spec_count++ if $run_number >= 0;
$spec_count++ if $library ne "";
$spec_count++ if $sql_select ne "";

if ($spec_count == 0) {
    print STDERR "All sequences will be annotated.\n";
}
elsif ($spec_count > 1) {
    die "You may specify only one value for -run, -library or -sql\n";
}
elsif ($library ne "") {
    if ($driver ne 'Pg') {
	die "Programming limitation: code only works for Postgres.\n";
    }
    $run_number = $pg->get_single_value("select number from contig_run " .
					" where library = " . quotify($library) .
					"   and create_time = (select max(create_time) " .
					"                        from contig_run " .
					"                       where library = " . quotify($library) . ")");
    if (not defined $run_number) {
	die "Unable find run number for $library\n";
    }
}

if ($run_number >= 0) {
    &check_run_number;
}

my $hit_id_value = $pg->nextval('hit_id_serial');
die "Unable to get access hit_id_serial\n" if not defined $hit_id_value;
my $hit_id = "H" . $hit_id_value . "_";

my $curtime;
if ($pg->driver eq 'Pg') {
    $curtime = $pg->get_single_value("select 'now'::timestamp with time zone");
}
else {
    $curtime = $pg->get_single_value("select datetime('now', 'localtime')");
}

$pg->command("insert into hit_id_record (hit_id_val, create_time) " .
	     "       values($hit_id_value, " . quotify($curtime) . ")");

my $dbfile = &set_dbfile;

my %seq_lengths;

Bio::Frescobi::PGLock::lock_wait_time(60); 
Bio::Frescobi::PGLock::lock_wait_cycles(20000);	# About two weeks.

my $run_count = $pg->nextval('annotator_run_count');

&select_data;
 
&run_annotator;

$pg->reset;

&load_results;

&update_annotator_history($common_db_name, $method, $hit_id_value);

system_with_check("rm -rf $workdir") if $delete_workdir;

print STDERR "Done at ", `date`;

sub check_run_number {
    my ($err);
	
    if (not $pg->table_exists("contig_$run_number")) {
	print STDERR "No table named contig_$run_number found\n";
	$err = 1;
    }
    if (not $pg->table_exists("singlets_$run_number")) {
	print STDERR "No table named singlets_$run_number found\n";
	$err = 1;
    }
    die if $err;
}

sub set_dbfile {
    # Retrieve the datafile from the database, or create one if needed.
    #
    
    my ($data, $i, $db_name, $file_name, $update_time, $dbfile, $dbfile2);
    my ($gencmd, $prepcmd);

    $data = $pg->query("select common_db_name, file_name, last_update_time, " .
		       "       internal, generation_command, prepare_command " .
		       "  from dbs " .
		       " where common_db_name = " . quotify($common_db_name));
    $dbfile = "";
    if (scalar(@{$data}) == 0) {
	die "No database found for $common_db_name\n";
    }
    elsif (scalar(@{$data}) > 1) {
	die "Multiple databases found for $common_db_name\n";
    }
    ($db_name, $file_name, $update_time, $internal, $gencmd, $prepcmd) = @{$data->[0]};
    if (defined($internal) and $internal eq 't') {
	$dbfile = "$workdir/$db_name";
	system_with_check("$gencmd -dbname=$dbname -driver=$driver $dbfile");
	system_with_check("$prepcmd $dbfile");
    }
    else {
	$internal = 0;
	$dbfile = $file_name;
	if ($dbfile eq "") {
	    die "Unable to find datafile for $common_db_name in the database.\n";
	}
	if (not -r $dbfile) {
	    die "The file, $dbfile, is not readable.\n";
	}
    }
    return $dbfile;
}

sub select_data {
    my ($data, $i, $seqid, %seqids_to_annotate, $sql);
    
    &Bio::Frescobi::PGLock::grab_lock($pg, "annotation_runs");
    
    if ($run_number >= 0) {
	$sql = "select seqid from contig_$run_number " .
	       " union " .
   	       "select seqid from singlets_$run_number ";
    }
    elsif ($sql_select ne "") {
	$sql = $sql_select;
    }
    else {
	$sql = "select seqid from seq_data where seqid !~ '^GP'";
    }
    $pg->command("create table annotator_selected_${run_count} " .
		 " as " . $sql);
    if ($replace) {
	$pg->command("delete from annotator_selected_${run_count} " .
		     " where seqid in (select seqid " .
		     "                   from annotation_runs " .
		     "                  where status = 'in progress' " .
		     "                    and common_db_name = " . quotify($common_db_name) .
		     "                    and method = " . quotify($method) . ")");
    }
    else {
	# If we're not replacing, then skip any annotations, either completed
	# or not in progress.
	$pg->command("delete from annotator_selected_${run_count} " .
		     " where seqid in (select seqid " .
		     "                   from annotation_runs " .
		     "                  where common_db_name = " . quotify($common_db_name) .
		     "                    and method = " . quotify($method) . ")");
    }	
    $data = $pg->query("select seqid from annotator_selected_${run_count} ");
    foreach my $rowp (@{$data}) {
	$seqid = $rowp->[0];
	$seqids_to_annotate{$seqid} = 1;
    }
    # If there are no seqids to annotate, the query below is
    # syntactically invalid, so we need to skip it in that case.

    if (scalar(keys %seqids_to_annotate) > 0) {
	$data = $pg->query("select seqid, length " .
			   "  from seq_data " .
			   " where seqid in (select seqid from annotator_selected_${run_count})");
	foreach my $rowp (@{$data}) {
	    $seqid = $rowp->[0];
	    my $length = $rowp->[1];
	    $seq_lengths{$seqid} = $length;
	}
    }
    if ($replace) {
	$data = $pg->query("select seqid from annotation_runs " .
			   " where common_db_name = " . &quotify($common_db_name) .
			   "   and method = " . &quotify($method) .
			   "   and seqid in (select seqid " .
			   "                   from annotator_selected_${run_count})" );
	foreach my $rowp (@{$data}) {
	    $seqid = $rowp->[0];
	    if (exists($seqids_to_annotate{$seqid})) {
		&delete_annotation($seqid, $common_db_name, $method);
	    }
	}
    }
    @seqids_to_annotate = sort { $seq_lengths{$b} <=> $seq_lengths{$a}} keys(%seqids_to_annotate);
    my $curtime;
    if ($pg->driver eq 'Pg') {
	$curtime = $pg->get_single_value("select 'now'::timestamp with time zone");
    }
    else {
	$curtime = $pg->get_single_value("select datetime('now', 'localtime')");
    }
    $pg->command("insert into annotation_runs " .
		 " (seqid, common_db_name, method, start_time, status) " .
		 "select seqid, " .
		 "   " . quotify($common_db_name) . " as common_db_name, " .
		 "   " . quotify($method) . " as method, " .
		 "   " . quotify($curtime) . " as start_time, " .
		 "       'in progress' as status " .
		 "  from annotator_selected_${run_count} ");
    &Bio::Frescobi::PGLock::free_lock($pg, "annotation_runs");
}

sub print_query {
    my $sql = shift;

    foreach my $row (@{$pg->get_all_rows($sql)}) {
	print STDERR join("\n", @{$row}), "\n";
    }
}

sub delete_annotation {
    my ($seqid, $common_db_name, $method) = @_;
    my ($data, $result2, $i, $j, $hit_id);

    my $seqid_select = "hits.seqid = '$seqid'";
    if ($internal) {
	$seqid_select = "($seqid_select or hits.key = '$seqid')";
    }
    $data = $pg->query("select hsps.hit_id from hits, hsps " .
		       " where hsps.hit_id = hits.hit_id and " .
		       "       $seqid_select and " .
		       "       common_db_name = " . &quotify($common_db_name) . " and " .
		       "       method = " . &quotify($method) .
		       " order by hsps.hit_id");
    my %hit_id_seen = ();
    my ($proc_id, $hit_serial, @hit_count_by_length, $hit_count);
    my $optimization_failed = 1;
    # my $prev_hit_serial = "";
    # my $prev_proc_id = "";
    # my $prev_hit_id = "";
    # for ($i = 0; $i < $result->ntuples; $i++) {
    # 	$hit_id = $result->getvalue($i, 0);
    # 	next if exists($hit_id_seen{$hit_id});
    # 	$hit_id_seen{$hit_id} = 1;
    # 	if ($hit_id !~ m/^H(\d+)_(\d+)_(\d+)$/) {
    # 	    print STDERR "Unable to parse hit_id: $hit_id\n";
    # 	    $optimization_failed = 1;
    # 	    last;
    # 	}
    # 	$hit_serial = $1;
    # 	$proc_id = $2;
    # 	$hit_count = $3;
    # 	if (($prev_hit_serial ne "" and $prev_hit_serial ne $hit_serial) or
    # 	    ($prev_proc_id ne "" and $prev_proc_id ne $proc_id)) {
    # 	    print STDERR "Unable to apply range find rules for hit_id ($hit_id) and previous hit_id ($prev_hit_id)\n";
    # 	    $optimization_failed = 1;
    # 	    last;
    # 	}
    # 	$prev_hit_id = $hit_id;
    # 	$prev_hit_serial = $hit_serial;
    # 	$prev_proc_id = $proc_id;
    # 	push (@{$hit_count_by_length[length($hit_id)]}, $hit_count);
    # }
    # unless ($optimization_failed) {
    #   GAP_CHECK:
    # 	for ($l = 1; $l < scalar(@hit_count_by_length); $l++) {
    # 	    next if not defined($hit_count_by_length[$l]);
    # 	    next if scalar(@{$hit_count_by_length[$l]}) == 0;
    # 	    @{$hit_count_by_length[$l]} = sort @{$hit_count_by_length[$l]};
    # 	    my $first_count = $hit_count_by_length[$l]->[0];
    # 	    my $prev_count = $first_count;
    # 	    for ($j = 1; $j < scalar(@{$hit_count_by_length[$l]}); $j++) {
    # 		if ($prev_count + 1 != $hit_count_by_length[$l]->[$j]) {
    # 		    printf STDERR "Successive hits are gapped. Hit_count = %s\n and previous count = %s\n",
    # 		    $hit_count_by_length[$l]->[$j], $prev_count;
    # 		    $optimization_failed = 1;
    # 		    last GAP_CHECK;
    # 		}
    # 		$prev_count = $hit_count_by_length[$l]->[$j];
    # 	    }
    # 	    my $first_hit_id = sprintf("H%s_%s_%s",
    # 				       $prev_hit_serial,
    # 				       $prev_proc_id,
    # 				       $first_count);
    # 	    my $last_hit_id = sprintf("H%s_%s_%s",
    # 				      $prev_hit_serial,
    # 				      $prev_proc_id,
    # 				      $prev_count);
    # 	    $pg->command("delete from hsps " .
    # 			 " where hit_id between '$first_hit_id' and '$last_hit_id' " .
    # 			 "   and length(hit_id) = " . length($first_hit_id));
    # 	    $pg->command("delete from hits " .
    # 			 "where hit_id between '$first_hit_id' and '$last_hit_id' " .
    # 			 "   and length(hit_id) = " . length($first_hit_id));
    # 	}
    # }
    if ($optimization_failed) {
	my $deletes = 0;
	%hit_id_seen = ();
	foreach my $rowp (@{$data}) {
	    $hit_id = $rowp->[0];
	    next if exists($hit_id_seen{$hit_id});
	    $hit_id_seen{$hit_id} = 1;
	    $pg->command("delete from hsps where hit_id = '$hit_id'", $debug);
	    $pg->command("delete from hits where hit_id = '$hit_id'", $debug);
	    $deletes += 2;
	}
	print STDERR "$deletes deletion commands executed.\n";
    }
    $pg->command("delete from annotation_runs " .
		 " where common_db_name = " . &quotify($common_db_name) . " and " .
		 "       method = " . &quotify($method) . " and " .
		 "       seqid = '$seqid'");
}

sub run_annotator {
    my ($pid, @pids, $file, $iproc, $seq, $seqid);
    my ($limit, $i, $err);

    $pg->close_connection;
    for ($iproc = 0; $iproc < $maxproc; $iproc++) {
 	if ($pid = fork) {
 	    print STDERR "Process $pid forked. iproc = $iproc.\n";
 	    push @pids, $pid;
 	}
 	elsif (defined $pid) {
	    sleep $iproc;
	    $pg->connect;
	    $hit_id .= $iproc;
	    $hit_count = 0;
	    open (NEWHITS, ">$workdir/load.hits.$iproc") || die "Unable to create $workdir/load.hits.$iproc: $!\n";
	    select((select(NEWHITS), $| = 1)[$[]);
	    open (NEWHSPS, ">$workdir/load.hsps.$iproc") || die "Unable to create $workdir/load.hsps.$iproc: $!\n";
	    select((select(NEWHSPS), $| = 1)[$[]);
	    open (NEWRUNS, ">$workdir/load.runs.$iproc") || die "Unable to create $workdir/load.runs.$iproc: $!\n";
	    select((select(NEWRUNS), $| = 1)[$[]);
	    $limit = scalar(@seqids_to_annotate);
	    # We want to distribute sequences so that each block has
	    # roughly the same amount of sequence.
	    my @blocks;
	    my $nblocks = POSIX::ceil(scalar(@seqids_to_annotate) / $block);
		
	    my $iblock = 0;
	    for (my $i = 1; $i <= $limit; $i++) {
		my $seqid = $seqids_to_annotate[$i-1];
		push (@{$blocks[$iblock++]}, $seqid);
		if ($iblock >= $nblocks) {
		    $iblock = 0;
		}
	    }
	    for (my $ib = 0; $ib < scalar(@blocks); $ib++) {
		if (($ib % $maxproc) == $iproc and
		    scalar(@{$blocks[$ib]}) > 0) {
		    &dated_mesg("search_one_block: iproc = $iproc  maxproc = $maxproc ib = $ib");
		    &search_one_block(@{$blocks[$ib]});
		}
	    }
	    print NEWHITS "\\.\n";
	    close NEWHITS;
	    print NEWHSPS "\\.\n";
	    close NEWHSPS;
	    print NEWRUNS "\\.\n";
	    close NEWRUNS;
	    exit 0;
	}
	else {
	    die "Can't fork: $!\n";
	}
    }
    $err = 0;
    foreach $pid (@pids) {
	waitpid($pid, 0);
	if ($? >> 8 != 0) {
	    printf STDERR "Process $pid return status was %x\n", $?;
	    $err = 1;
	}
    }
    carp "Subprocess failed unexpectedly. \n" if $err;
    $pg->connect;
}

sub write_table {
    my ($data, $fh) = @_;
    foreach my $line (@{$data}) {
	print $fh $line;
    }
}

sub search_one_block {
    my @seqids = @_;

    local (*SEQ);
    my ($file);

    $file = &make_tmp_seq_file_name($seqids[0]);
    open (SEQ, ">$file") || die "Unable to open for write $file: $!\n";
    my $good_count = 0;
    foreach my $seqid (@seqids) {
	my ($seq) = retrieve_sequence_data($pg, $seqid);
	if (&is_totally_masked($seq) or $seq_lengths{$seqid} < $minlength) {
	    my $runline = sprintf("%s\t%s\t%s\t%s\t%s\n",
				  $seqid, 
				  $common_db_name,
				  $method,
				  $curtime,
				  'completed');
	    print NEWRUNS $runline;
	}
	else {
	    print SEQ ">$seqid\n";
	    print SEQ &format_seq($seq);
	    $good_count += 1;
	}
    }
    close SEQ;
    if ($good_count > 0) {
	if ($method eq "fasty") {
	    search_with_fasta($method, $file, @_);
	}
	elsif ($method eq "ggsearch") {
	    search_with_fasta($method, $file, @_)
	}
	elsif ($method =~ /^t?blast[xpn]$/ and $method ne "tblastp") {
	    search_with_blast($method, $file, @_);
	}
	elsif ($method =~ /^blast\+t?blast[xpn]$/ and $method ne "blast\+tblastp") {
	    search_with_blast($method, $file, @_);
	}
	elsif ($method =~ m/^exonerate_(protein2dna|dna2dna|cdna2genome)$/) {
	    my $exonerate_method = substr($method, length("exonerate_"));
	    search_with_exonerate($exonerate_method, $minscore, $file, @_);
	}
	elsif ($method eq "hmm") {
	    search_with_hmm($file, @_);
	}
	else {
	    die "Search method $method not yet implemented.\n";
	}
    }
    else {
	printf STDERR "For block containing %s, no sequences were annotated because of size or total masking.\n",
	join(" ", @seqids);
    }
#    system("rm -f $file"); # Zero length sequences will not make a file so
                           # don't check for failure.
}

sub make_tmp_seq_file_name {
    my ($seqid) = @_;
    my $file = $workdir . "/seq.$seqid.$$";
    return $file;
}

sub search_with_fasta {
    my $method = shift;
    my $file = shift;
    my @seqids = @_;

    my ($cmd, $date, $count, @words, $key, $score, $description, $found);
    local ($_);
    my ($bo, $out_file, $empty, $hits, $runline);
    my ($num_identical, $num_conserved, $length, $frac_identical, $frac_conserved);
    $out_file = $workdir . "/fasta.out.$$";
    my $args = $extra_options;
    my $program;
    if ($method eq "fasty") {
	$program = $fasty_program;
    }
    elsif ($method eq "ggsearch") {
	$program = $ggsearch_program;
    }
    else {
	die "Unknown method: $method specified for fasta analysis.\n";
    }
    if (defined $gapopen) {
	$args .= " -f $gapopen";
    }
    my @q_seqid = map { quotify($_) } @seqids;
    my @seq_type = $pg->get_array_for_field("select distinct seq_type" .
					    "  from seq_data " .
					    " where seqid in (" . join(",", @q_seqid) . ")");
    if (scalar(@seq_type) == 0) { 
	print STDERR "Programming error in search_with_fasta: seq_type array is empty.\n";
	exit(1);
    }
    elsif (scalar(@seq_type) > 1) {
	print STDERR "Error in search_with_fasta: seq_types must all be the same.\n";
	exit(1);
    }
    if ($seq_type[0] eq 'protein') {
	$args .= " -p";
    }
    elsif ($seq_type[0] eq "nucleic") {
	$args .= " -n";
    }
    $cmd = "$nice_command $program $args -E $cutoff -w 200 -b $maxhits -d $maxhits -L -H -q -Z $genome_size $file $dbfile >$out_file";
    $date = `date`;
    chomp $date;
    print STDERR "At $date: start $cmd\n";
    $found = 0;
    eval {
	alarm $timeout;
	system_with_check($cmd);
	alarm 0;
    };
    alarm 0;
    if ($@) {
	if ($@ !~ /timeout/) {
	    die $@;
	}
	print STDERR "At $date: timeout occurred on $cmd\n";
	$out_file = "/dev/null";
    }
    my $parser = new Bio::SearchIO(-format => "fasta",
				   -file => $out_file);
    if (not $parser) {
	print STDERR "Unable to open $out_file: $!\n";
	print NEWRUNS $runline;
	return;
    }
    while (my $result = $parser->next_result) {
	$count = 0;
	my $seqid = $result->query_name;
	$runline = sprintf("%s\t%s\t%s\t%s\t%s\n",
			   $seqid, 
			   $common_db_name,
			   $method,
			   $curtime,
			   'completed');
	while (my $hit = $result->next_hit) {
	    while (my $hsp = $hit->next_hsp) {
		$key = $hit->name;
		next if (not $selfself) and $seqid eq $key;
		my $desc = $hit->description;
		my $q_strand = $hsp->strand('query') eq '1' ? 't' : 'f';
		my $s_strand = 't';
		my $expect = $hsp->evalue;
		my $score = $hsp->sw_score;
		my $bits = $hsp->bits;
		my $pct_id = $hsp->frac_identical * 100;
		$length = $hsp->length;
		my ($q_start, $q_end) = $hsp->range('query');
		my ($s_start, $s_end) = $hsp->range('hit');
		next if $expect > $cutoff and $expect != 999;
		next if $pct_id / 100.0 < $idcutoff;
		$hit_count += 1;
		my $new_hit_id = $hit_id . "_" . $hit_count;
		$desc =~ s/\t/ /g;
		$desc =~ s/\s\s+/ /g;
		&check_length($desc, $seqid, $key);
		$count += 1;
		printf NEWHITS "%s\t%s\t%s\t%s\t%s\t%g\t%s\n",
		$seqid, $common_db_name, $method, $key, $desc, $expect, $new_hit_id;
		$num_identical = $hsp->num_identical;
		$num_conserved = $hsp->num_conserved;
		$frac_conserved = $hsp->frac_conserved;
		$frac_identical = $hsp->frac_identical;
		printf NEWHSPS "%s\t%d\t%g\t%g\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%s\t%d\t%d\t%d\t%d\n",
		$new_hit_id, $score, $bits, $expect, $length,
		0, $num_identical, $num_conserved,
		$frac_identical, $frac_conserved, 
		$q_strand, $s_strand,
		$q_start, $q_end,
		$s_start, $s_end;
	    }
	}
	print STDERR "At $date: $count matches found for $seqid.\n";
	print NEWRUNS $runline;
    }
}

sub search_with_hmm {
    my $file = shift;
    my @seqids = @_;

    my $cmd = "hmmpfam $extra_options $dbfile $file >$file.out";
    my $date = `date`;
    chomp $date;
    print STDERR "At $date: $cmd\n";
    eval {
	alarm $timeout;
	system($cmd);
	alarm 0;
    };
    alarm 0;
    if ($@) {
	if ($@ !~ /timeout/) {
	    die $@;
	}
	print STDERR "At $date: timeout occurred on $cmd\n";
	return;
    }
    my @queries = &parse_hmm("$file.out");
    if (scalar(@queries) == 0) {
	return;
    }
    my $count = 0;
    foreach my $query (@queries) {
	my $seqid = $query->{QUERY_ID};
	my $runline = sprintf("%s\t%s\t%s\t%s\t%s\n",
			      $seqid, 
			      $common_db_name,
			      $method,
			      $curtime,
			      'completed');
	my $count = 0;
	foreach my $hit (@{$query->{DOMAINS}}) {
	    $count += 1;
	    $hit_count += 1;
	    my $new_hit_id = $hit_id . "_" . $hit_count;
	    $count += 1;
	    printf NEWHITS "%s\t%s\t%s\t%s\t%s\t%g\t%s\n",
	    $seqid, $common_db_name, $method,
	    $hit->{MODEL},
	    $hit->{MODEL},
	    $hit->{EVALUE},
	    $new_hit_id;
	    my $len_seq = $hit->{SEQT} - $hit->{SEQF} + 1;
	    my $len_hmm = $hit->{HMMT} - $hit->{HMMF} + 1;
	    my $length = max($len_seq);
	    my $num_conserved = 0;
	    my $num_identical = 0;
	    my $frac_conserved = $num_conserved / $length;
	    my $frac_identical = $num_identical / $length;
	    printf NEWHSPS "%s\t%d\t%g\t%g\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%s\t%d\t%d\t%d\t%d\n",
	    $new_hit_id,
	    $hit->{SCORE},
	    $hit->{SCORE},
	    $hit->{EVALUE},
	    $length,
	    0,
	    $num_identical,
	    $num_conserved,
	    $frac_identical,
	    $frac_conserved,
	    't',
	    't',
	    $hit->{SEQF},
	    $hit->{SEQT},
	    $hit->{HMMF},
	    $hit->{HMMT};
	}
	print NEWRUNS $runline;
	print STDERR "At $date: $cmd $count domains found for $seqid.\n";
    }
#     &system_with_check("rm $file.out");
}

sub search_with_blast {
    my $method = shift;
    my $file = shift;
    my @seqids = @_;
    
    my ($cmd, $date, $count, @words, $key, $score, $description, $found);
    local (*IN, *OUT, $_);
    my ($parser, $out_file, $empty, $hits, $hsp_count, $runline);
    my ($num_identical, $num_conserved, $length, $frac_identical, $frac_conserved);
    my ($status, $name);
    my $filter_opt = "";
    
    # This logical variable specifies if we are using the old blastall
    # or a more recent Blast+
    my $is_blastall = $method =~ m/^blast\+/ ? 0 : 1;

    my $protein_search = $method =~ m/[xp]$/ ? 1 : 0;

    $out_file = $workdir . "/blast.out.$$.$seqids[0]";
    if ($is_blastall) {
	if (defined $filter) {
	    if ($filter) {
		$filter_opt = "-F T";
	    }
	    else {
		$filter_opt = "-F F";
	    }
	}
	else {
	    $filter_opt = $protein_search ? "-F F" : "";
	}
    }
    else {
	if (not defined $filter) {
	    $filter = 0;
	}
	if ($filter and $protein_search) {
	    $filter_opt = "-seg yes";
	}
	elsif ($filter and not $protein_search) {
	    $filter_opt = "-dust yes";
	}
    }
	    
    my $wordsize_opt = "";
    if ($wordsize > 0) {
	$wordsize_opt = $is_blastall ? "-W $wordsize" : "-word_size $wordsize";
    }
    if ($is_blastall) {
	$cmd = "$nice_command $workdir/blastall " .
	    "$wordsize_opt " .
	    "$extra_options " .
	    "-p $method " .
	    "-i $file " .
	    "-d $dbfile " .
	    "-v $maxhits " .
	    "-b $maxhits " .
	    "-z $genome_size " .
	    "-e $cutoff " .
	    "-a 1 " .
	    "$filter_opt " .
	    ">${out_file}.tmp";
    }
    else {
	my $method_binary;
	($method_binary = $method) =~ s/^blast\+//;
	if ($protein_search) {
	    $cmd = "";
	}
	else {
	    $cmd = "$nice_command $workdir/blast+bin/${method_binary} " .
		"$wordsize_opt " .
		"$extra_options " .
		"-query $file " .
		"-db $dbfile " .
		"-num_alignments $maxhits " .
		"-num_descriptions $maxhits " .
		"-dbsize $genome_size " .
		"-evalue $cutoff " .
		"-num_threads 1 " .
		"$filter_opt " .
		">${out_file}.tmp";
	}
    }
    eval {
	alarm $timeout;
	$status = system($cmd);
	alarm 0;
    };
    alarm 0;
    if ($@) {
	if ($@ !~ /timeout/) {
	    die $@;
	}
	print STDERR "At $date: timeout occurred on $cmd\n";
	$empty = 1;
    }
    else {
	if ($status != 0) {
	    print STDERR sprintf("$cmd returned %d error code.\n", $status);
	    $empty = 1;
	}
	else {
	    $empty = 1;
	    open (IN, "<${out_file}.tmp") ||
		die "Unable to open ${out_file}.tmp: $!\n";
	    open (OUT, ">$out_file") ||
		die "Unable to open $out_file: $!\n";
	    my @lines;
	    my $hits_found = 1;
	    while (<IN>) {
		push (@lines, $_);
		if (m/^S2/) {
		    if ($hits_found) {
			# There's a formatting problem with blastall where alignments with lots of
			# gaps put the gap percentage on the following line and that messes
			# up the BioPerl parser. The follow code rejoins the percentage.
			
			for (my $i = 1; $i < scalar(@lines); $i++) {
			    if ($lines[$i-1] =~ m/^ Identities =/ and
				$lines[$i] =~ m/^\(\d+%\)/) {
				chomp $lines[$i-1];
				$lines[$i-1] .= " " . $lines[$i];
				for (my $j = $i + 1; $j < scalar(@lines); $j++) {
				    $lines[$j-1] = $lines[$j];
				}
				$#lines = $#lines - 1;
			    }
			}
			print OUT @lines;
			$empty = 0;
		    }
		    else {
			$hits_found = 1;
		    }
		    @lines = ();
		}
		if (m/\*\*\*\*\* No hits found \*\*\*\*\*\*/) {
		    $hits_found = 0;
		}
	    }
	    close IN;
	    print OUT @lines if scalar(@lines) > 0;
	    close OUT;
	}
    }
    $date = `date`;
    chomp $date;
    if (-z $out_file) {
	print STDERR "At $date: $cmd 0 sequences found\n";
    }
    else {
	$parser = new Bio::SearchIO(-format => "blast",
				    -file => $out_file);
	if (not $parser) {
	    print STDERR "Unable to open $out_file: $!\n";
	}
	else {
	    while (my $result = $parser->next_result) {
		my $seqid = $result->query_name;
		$count = 0;
		$hsp_count = 0;
		while (my $hit = $result->next_hit) {
		    next if $hit->expect > $cutoff;
		    if (scalar($hit->hsps) > 0 and
			($hit->hsps)[0]->frac_identical  < $idcutoff) {
			next;
		    }
		    my $name = $hit->name;
		    next if not $selfself and $name eq $seqid;
		    $name =~ s/^[^:]+://g;
		    next if not $selfself and $name eq $seqid;
		    $hit_count += 1;
		    my $new_hit_id = $hit_id . "_" . $hit_count;
		    my $description = $hit->description;
		    $description =~ s/\t/ /g;
		    $description =~ s/\s\s+/ /g;
		    $count += 1;
		    &check_length($description, $seqid, $name);
		    printf NEWHITS "%s\t%s\t%s\t%s\t%s\t%g\t%s\n",
		    $seqid, $common_db_name, $method, $name, $description, $hit->expect,
		    $new_hit_id;
		    foreach my $hsp ($hit->hsps) {
			next if $hsp->expect > $cutoff;
			$hsp_count += 1;
			my $length = $hsp->length;
			my $num_identical = $hsp->num_identical;
			my ($num_conserved, $frac_identical, $frac_conserved);
			if ($method eq "blastn") {
			    $num_conserved = $num_identical;
			}
			else {
			    $num_conserved = $hsp->num_conserved;
			}
			if ($length == 0) {
			    $frac_identical = 0.0;
			    $frac_conserved = 0.0;
			}
			else {
			    $frac_identical = $num_identical / $length;
			    $frac_conserved = $num_conserved / $length;
			}
			printf NEWHSPS "%s\t%d\t%g\t%g\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%s\t%d\t%d\t%d\t%d\n",
			$new_hit_id, $hsp->score, $hsp->bits, $hsp->expect, $length,
			$hsp->gaps, $num_identical, $num_conserved,
			$frac_identical, $frac_conserved, 
			($hsp->strand('query') eq "Plus" or $hsp->strand('query') == 1) ? 't' : 'f',
			($hsp->strand('sbjct') eq "Plus" or $hsp->strand('sbjct') == 1) ? 't' : 'f',
			$hsp->start('query'), $hsp->end('query'),
			$hsp->start('sbjct'), $hsp->end('sbjct');
		    }
		}
		print STDERR "At $date: $cmd $count descriptions and $hsp_count HSP's found for $seqid.\n";
	    }
	}
    }
#    &system_with_check("rm -f $out_file ${out_file}.tmp");
    foreach my $seqid (@seqids) {
	$runline = sprintf("%s\t%s\t%s\t%s\t%s\n",
			   $seqid, 
			   $common_db_name,
			   $method,
			   $curtime,
			   'completed');
	print NEWRUNS $runline;
    }
}

sub load_results {
    my ($iproc, @words, $seqobj);
    my ($n_seqid, $seq, $rmethod, $gene);
    local ($_, *IN);

    &Bio::Frescobi::PGLock::grab_lock($pg, "annotation_runs");
    my @seqids = $pg->get_array_for_field("select seqid from annotator_selected_${run_count}");
    for (my $is = 0; $is < scalar(@seqids); $is += 1000) {
	my @list = @seqids[$is..min($is+999, scalar(@seqids)-1)];
	my @quoted_list = map { quotify($_) }  @list;
	$pg->command("delete from annotation_runs " .
		     " where seqid in (" . join(",", @quoted_list) . ")" . 
		     "   and common_db_name = " . quotify($common_db_name) .
		     "   and method = " . quotify($method));
    }
    for ($iproc = 0; $iproc < $maxproc; $iproc++) {
	&load_file("$workdir/load.runs.$iproc", "annotation_runs");
	if (not $internal) {
	    &load_file("$workdir/load.hits.$iproc", "hits");
	    &load_file("$workdir/load.hsps.$iproc", "hsps");
	}
    }
    if ($internal) {
	&load_internal_results;
    }
    $pg->command("drop table annotator_selected_${run_count}") unless $debug;
    &Bio::Frescobi::PGLock::free_lock($pg, "annotation_runs");
}

sub load_file {
    my $file = shift;
    my $table = shift;
    local (*IN, $_);

    open (IN, "<$file") ||
	die "Unable to open $file: $!\n";
    $pg->startcopy($table);
    while (<IN>) {
	last if /^\\\.$/;
	$pg->putline($_);
    }
    $pg->endcopy();
    close IN || die "Unable to close $file: $!\n";
}    
    
sub search_with_exonerate {

    # This code is intended for searches of protein against dna. As such, the order
    # of query and database is reversed, so the retrieval of information from
    # the results is inverted.

    my $model = shift;
    my $minscore = shift;
    my $file = shift;
    my @seqids = @_;

    my ($cmd, $date);
    local (*IN, $_);
    my ($out_file, $runline);
    my ($count);

    $date = &date;
    $out_file = $workdir . "/exonerate.out.$$";
    my $extra = $extra_options;
    if (defined $minscore) {
	$extra .= " --score $minscore";
    }
    if ($model =~ m/protein2dna/) {
	$extra .= " --frameshift -5 --proteinwordlen 3 --dnawordlen 6";
    }
    elsif ($model =~ m/dna2dna/) {
	$model = 'affine:bestfit';
	$extra .= " --exhaustive";
    }
    $cmd = "sh -c \"$nice_command $exonerate_program --showvulgar 1 --showcigar 0 " .
        "--bestn $maxhits " . 
	"$extra --model $model $dbfile $file >$out_file \"";
    eval {
	alarm $timeout;
	system_with_check($cmd, "echo");
	alarm 0;
    };
    alarm 0;
    if ($@) {
	if ($@ !~ /timeout/) {
	    die $@;
	}
	print STDERR "At $date: timeout occurred on $cmd\n";
	return;
    }
    my $parser = new Bio::SearchIO(-format => "exonerate",
				   -file => $out_file);
    if (not defined $parser) {
	print STDERR "Unable to open exonerate output file $out_file: \n";
	return;
    }
    $count = 0;
    my $hsp_count = 0;
    while (my $result = $parser->next_result) {
	while (my $hit = $result->next_hit) {
	    $count += 1;
	    $hit_count += 1;
	    my $seqid = $hit->name;
	    my $new_hit_id = $hit_id . "_" . $hit_count;
	    my $description = $result->query_description;
	    $description =~ s/\t/ /g;
	    $description =~ s/\s\s+/ /g;
	    my $name = $result->query_name;
	    my $best_score = 0;
	    while (my $hsp = $hit->next_hsp) {
		$name =~ s/^[^:]+://g;
		&check_length($description, $seqid, $name);
		$hsp_count += 1;
		my $length = $hsp->length('query');
		my $num_identical = -1;
		my $num_conserved = -1;
		my $frac_identical = -1;
		my $frac_conserved = -1;
		$best_score = max($best_score, $hsp->score);
		printf NEWHSPS "%s\t%d\t%g\t%g\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%s\t%d\t%d\t%d\t%d\n",
		$new_hit_id, $hsp->score, -1, 0.0, $length,
		$hsp->gaps, $num_identical, $num_conserved,
		$frac_identical, $frac_conserved, 
		($hsp->strand('hit') eq "Plus" or $hsp->strand('hit') == 1) ? 't' : 'f',
		($hsp->strand('query') eq "Plus" or $hsp->strand('query') == 1) ? 't' : 'f',
		$hsp->start('hit'), $hsp->end('hit'),
		$hsp->start('query'), $hsp->end('query');
	    }
	    printf NEWHITS "%s\t%s\t%s\t%s\t%s\t%g\t%s\n",
	    $seqid, $common_db_name, $method, $name, $description, $best_score, $new_hit_id;
	}
    }
    foreach my $seqid (@seqids) {
	my $runline = sprintf("%s\t%s\t%s\t%s\t%s\n",
			      $seqid, 
			      $common_db_name,
			      $method,
			      $curtime,
			      'completed');
	print NEWRUNS $runline;
    }
#     &system_with_check("rm $out_file");
    print STDERR "At $date: $cmd $count descriptions and $hsp_count HSP's found.\n";
}

sub load_internal_results {
    my ($iproc);
    local (*PIPE, *IN, $_);
    my (@words);

    my $tmpfile = "$workdir/load.hits.tmp";
    open (PIPE, "| sort -T $workdir >$tmpfile");
    my %skip_hit_id = ();
    my $skip_identity = 0;
    for ($iproc = 0; $iproc < $maxproc; $iproc++) {
	my $file = "$workdir/load.hits.$iproc";
	open (IN, "<$file") || die "Unable to open $file: ". $!. "\n";
	dated_mesg("Now reading $file...");
	while (<IN>) {
	    chomp;
	    next if m/^\\\.$/;
	    @words = split(/\t/);
	    #
	    # Skip the self comparison.
	    #
	    if ($words[0] eq $words[3]) {
		$skip_hit_id{$words[6]} = 1;
		$skip_identity += 1;
		next;
	    }
	    # Remove the duplicates by keeping the hit ordered.
	    next if $words[0] gt $words[3];
	    print PIPE join("\t", @words), "\n";
	}
	close IN || die "Unable to close $file: $!\n";
    }
    dated_mesg("$skip_identity identity hits skipped");
    close PIPE || die "Unable to close pipe: $!\n";
    open (IN, "<$workdir/load.hits.tmp") ||
	die "Unable to open $workdir/load.hits.tmp: $!\n";
    my $prev_seqid = "";
    my $prev_key = "";
    my $dupl_count = 0;
    my $hit_count = 0;
    my $last_dupl = 0;
    $pg->startcopy("hits");
    while (<IN>) {
	chomp;
	@words = split(/\t/);
	if ($prev_seqid eq $words[0] and
	    $prev_key eq $words[3]) {
	    if ($last_dupl) {
		print STDERR "The following hit was found more than twice: $_\n";
	    }
	    $last_dupl = 1;
	    $dupl_count += 1;
	    $skip_hit_id{$words[6]} = 1;
	    next;
	}
	$last_dupl = 0;
	$prev_seqid = $words[0];
	$prev_key = $words[3];
	$pg->putline("$_\n");
	$hit_count += 1;
    }
    $pg->endcopy();
    close IN;
    dated_mesg("$hit_count hits loaded and $dupl_count duplicate hits skipped.");
    my $hit_id;
    my $hsps_skipped = 0;
    my $hsps_loaded = 0;
    for ($iproc = 0; $iproc < $maxproc; $iproc++) {
	my $file = "$workdir/load.hsps.$iproc";
	open (IN, "<$file") || die "Unable to open $file: ". $!. "\n";
	dated_mesg("Now reading $file...");
	$pg->startcopy("hsps");
	while (<IN>) {
	    last if /^\\\.$/;
	    $hit_id = (split(/\t/, $_, 2))[0];
	    if (exists($skip_hit_id{$hit_id})) {
		$hsps_skipped += 1;
		next;
	    }
	    $pg->putline($_);
	    $hsps_loaded += 1;
	}
	$pg->endcopy();
	close IN || die "Unable to close $file: $!\n";
    }
    dated_mesg("$hsps_loaded HSPS were loaded.");
    dated_mesg("$hsps_skipped HSPS were skipped because of duplication.");
}
    
sub check_length {
    my ($desc, $seqid, $key) = @_;
    
    if (length($desc) > 7000) {
	$_[0] = substr($desc, 0, 7000) . " Truncated...";
	print STDERR "For $seqid, description for $key has been truncated.\n";
    }
}

sub update_annotator_history {
    my ($common_db_name, $method, $hit_id_value) = @_;
    if ($run_number > -1) {
	$pg->autocommit(0);
	$pg->command("begin");
	$pg->command("delete from annotator_history " .
		     " where library_number = $run_number and " .
		     "       common_db_name = " . quotify($common_db_name) . " and " .
		     "       method = " . quotify($method));
	$pg->command("insert into annotator_history " .
		     "       (library_number, common_db_name, method, hit_id_value) " .
		     " values($run_number, " .
		     "       " . quotify($common_db_name) . ", " .
		     "       " . quotify($method) . ", " .
		     "       $hit_id_value)");
	$pg->command("end");
	$pg->autocommit(1);
    }
}

sub parse_hmm {
    my $file = shift;
    local (*IN, $_);
    my @ret;

    my %threshold = (A_DOMAIN => 50,
		     C_DOMAIN => 100,
		     E_DOMAIN => 200,
		     Cy_DOMAIN => 100,
		     R_DOMAIN => 300,
		     T_DOMAIN => 10,
		     M_DOMAIN => 200,
		     TE_DOMAIN => 8.3,
		     AT_DOMAIN => 100,
		     DH_DOMAIN => 0.5,
		     ER_DOMAIN => 200,
		     KR_DOMAIN => 8,
		     KS_DOMAIN => 180);

    open (IN, "<$file") ||
	die "Unable to open $file: $!\n";
    while (<IN>) {
	chomp;
	last if m/^Query sequence/;
    }
    my $iquery = 0;
    while (defined $_) {
	my $query_id = substr($_, length("Query sequence: "));
	trim($query_id);
	$ret[$iquery]->{QUERY_ID} = $query_id;
	$ret[$iquery]->{DOMAINS} = [ ];
	while (<IN>) {
	    last if m/^Parsed for domains/;
	    last if m#^//#;
	}
	last if not defined $_;
	unless (m#^//#) {
	    while (<IN>) {
		last if m/^-----/;
	    }
	    last if not defined $_;
	    while (<IN>) {
		my ($model, $domain, $seqf, $seqt, $seq_bound,
		    $hmmf, $hmmt, $hmm_bound, $score, $evalue) = split(/\s+/);
		last if not defined $evalue;
		my $ok;
		if (exists($threshold{$model})) {
		    $ok =  ($score > $threshold{$model});
		}
		else {
		    print STDERR "No threshold found for $model. File is $file.\n";
		    $ok = 1;
		}
		if ($ok) {
		    push (@{$ret[$iquery]->{DOMAINS}}, { MODEL => $model,
							 DOMAIN => $domain,
							 SEQF => $seqf,
							 SEQT => $seqt,
							 SEQ_BOUND => $seq_bound,
							 HMMF => $hmmf,
							 HMMT => $hmmt,
							 HMM_BOUND => $hmm_bound,
							 SCORE => $score,
							 EVALUE => $evalue });
		}
	    }
	}
	while (<IN>) {
	    chomp;
	    last if m/^Query sequence/;
	}
	$iquery++;
    }
    return @ret;
}

=head1 NAME

annotator.pl - Annotate sequences in the database.

=head1 SYNOPSIS

 annotator.pl -dbname=name
             [-driver=[Pg|SQLite] ]
             [-db=db]         
             [-method=executable]
             [-extraopts=<string>]
             [-run=integer]
             [-library=string]
             [-tmp=dir]
             [-maxproc=integer]
             [-maxhits=integer]
             [-cutoff=float]
             [-idcutoff=float]
             [-replace]
             [-timeout=integer]
             [-genome_size=integer]
             [-sql=string]
             [-[no]filter]
             [-minscore=integer]
             [-gapopen=integer]
             [-block=integer]  Default: 1
             [-[no]debug]
             [-minlength=integer]
             [-wordsize=integer]
             [-[no]selfself]
             [-nice=string]
             [-fasty=string]
             [-ggsearch=string]
             [-exonerate=string]

    Default driver is "Pg", i.e. Postgres.
    Default db is nr.
    Default method is blastn
    Either run, library, sql can be specified, but not more than one.
    Default sql is blank which means all sequences.
    Default tmp is /tmp/annotator.<pid>
    Default maxproc = 1
    Default cutoff is 0.1
    Default idcutoff is 0.0
    dbname must be specified.
    Default maxhits is 1000
    Default timeout is 31415926 seconds (about a year)
    -replace means replacing any existing annotations.
    Default genome size is 2000000000
    Default filter value is blastall default.
    Default minscore is undefined.
    Default gapopen is undefined.
    Default block is 1
    Default debug is false.
    Default minlength is 6. Lower will not work with current version of blastall.
    Default wordsize is up to blastall
    Default nice (command) is "nice -19"
    Default fasty (program) is "fasty35"
    Default ggsearch (program) is "ggsearch35"
    Default exonerate (program) is "exonerate";

=head1 DESCRIPTION

The C<annotator.pl> script is the workhorse script for annotating
sequences. It selects a set of sequences from the database, divides
the set into smaller sets for running the annotations in parallel,
runs the similarity algorithms in parallel, and then load the results
into the database. The script keeps track of what sequences have
already been annotated, and it can do incremental annotations, i.e.,
only annotate those sequences that have not been annotated earlier. It
is also possible (using the C<-replace> option) to replace existing
annotations.

=head1 ARGUMENTS AND OPTIONS

=over 4

=item C<< -dbname=name >>

Specifies the database name from where the sequences will be read and
where the results of the annotation will be written. There are
currently two database systems that Frescobi can use, PostgreSQL and
SQLite. In the case of PostgreSQL, you must specify the name of the
database in this option, and you must configure environment variables
so the script can login. In the case of SQLite, the C<-dbname> option
is used to specify the filename of the SQLite database.

This option is mandatory.

=item C<< -driver=[Pg|SQLite]  >>

Specifies the database engine to be used for data retrieval and
storage. See the C<-dbname> option above.

=item C<< -db=db >>

The C<-db> option specifies the sequence database for use as a
reference for the sequence annotation. This database is different than
the database for storing information. The database must be found in
the C<common_db_name> field of the C<dbs> table.

=item C<< -method=string >>

Specifies the sequence similarity method to be used for
annotation. The method string is not case sensitive.

Please consult
the source code for the details of the method and how the results are
stored. Here is the list of methods:

=over 4

=item C<< blastn, blastx, blastp, tblastn, tblastx >>

Use the old C<blastall> program from NCBI. The specific search methods are as follows:

=over 4

=item blastn

For a nucleotide query, search a nucleotide database for similar sequences.

=item blastp

For a protein query, search a protein database for similar sequences.

=item blastx

For a nucleotide query, search all three frame translation of the
query in a protein database for matching sequences.

=item tblastn

For a protein query, search over the three frame translations in a nucleotide database.

=item tblastx

For a nucleotide query, search all three translations into protein
over the three frame translations in a nucleotide database. This is
the most expensive search.

=back

=item C<< blast+blastn, blast+blastx, blast+blastp, blast+tblastn, blast+tblastx >>

Use the current Blast+ programs from NCBI. The search methods are the
same as with the C<blastall> program above -- the method names have a
prefix of "blast+" to distinguish using the new method versus the old
method.

=item C<< fasty, ggsearch >>

Use Pearson's C<fasty> or C<ggsearch> program for protein sequence
comparison. Reference is William R. Pearson, Todd Wood, Zheng Zhang,
Webb Miller, "Comparison of DNA Sequences with Protein Sequences",
Genomics, Volume 46, Issue 1, 1997, pages 24-36. The options,
C<-fasty> and C<-ggsearch>, below will specify the paths to get to
these programs.

=item C<< exonerate_protein2dna, exonerate_dna2dna, exonerate_cdna2genome  >>

Use the Exonerate program, (Slater, G.S.C., Birney, E. Automated
generation of heuristics for biological sequence comparison. BMC
Bioinformatics 6, 31 (2005)) for sequence comparison.

Please note that in this version of Frescobi, there is no testing for
this option. Such testing will be provided in a future release.

=item C<< hmm >>

Use the Hidden Markov Modeling program, C<hmmpfam>, to search a family
of models for protein sequence similarity. Here is one reference:
Eddy, S. R. (1998). Profile hidden Markov models. Bioinformatics, 14:755-763.

Frescobi was written to use HMMER version 2.3.2.

The database is a PFAM file.

=back

=item C<< -extraopts=<string> >>

The C<-extraopts> option is used to add extra options to the command
line for running the search program above.. This provides you great
flexibility in running the program.

=item C<< -run=integer >>

This option is obsolete -- it was used when Frescobi could handle assemblies of short reads.

=item C<< -library=string >>

This option is obsolete -- it was used when Frescobi could handle assemblies of short reads.

=item C<< -tmp=dir >>

This option specifies a directory for temporary files. For a large run, this directory needs to be a file system with plenty of space. It is good practice to use this option with a specific directory, and use it to diagnose failures in the annotation process. The default is C<< /tmp/annotator.<process id> >>.

=item C<< -maxproc=integer >>

Specifies the maximum number of processes to run annotation. Please note that the annotator uses a Linux fork system call to create sub processes for parallel execution, so the system needs adequate cores and memory to support all the running processes. The default setting is 1.

=item C<< -maxhits=integer >>

Specifies the maximum number of hits to be saved in the database for each query sequence. 

=item C<< -cutoff=float >>

Specifies the maximum value for the search method scoring function for keeping hits. For example, this parameter is the limit for the E value for Blast. The default value is 0.1, which is probably too big for Blast algorithms.

=item C<< -idcutoff=float >>

Specifies the maximum identical sequence ratio for keeping a hit. The default is 0, so all matches meeting the C<-cutoff> criteria above are kept. This parameter must be a floating point number between 0 and 1.

=item C<< -[no]replace >>

Specifies whether or not existing annotations for sequences should be
replaced. If C<-noreplace> is specified, then only those sequences
which have not been annotated using the specified sequence database
(C<-db> option) and method (C<-method> option) will be annotated. If
C<-replace> is specified, then all existing sequence annotations using
the sequence database and method will be replaced with new ones.

=item C<< -timeout=integer >>

This option gives the user control over long running search
algorithms. For example, Blast can take an inordinately long time when
searching reference sequences with high levels of repetitive
sequences.  The C<-timeout> option specifies the maximum number of
seconds that any search can execute. If the timeout occurs, execution
terminates.

=item C<< -genome_size=integer >>

Specifies the database size option for search methods that use it, such as Blast. The units are letters in the reference sequences.

=item C<< -sql=string >>

The C<-sql> option is used to specify which sequences will be
annotated. The SQL statement much have a C<seqid> column in the
selection list.

=item C<< -[no]filter >>

The C<-filter> option is used to control if low complexity filtering is applied to any Blast search. It is not used for other search methods.

=item C<< -minscore=integer >>

The C<-minscore> option specifies the minimum score for an Exonerate
hit. It does not apply to any other search method.

=item C<< -gapopen=integer >>

The C<-gapopen> parameter sets the gap opening penalty (C<-f> option)
for Pearson's fasty or ggsearch programs. Use the C<-extraopts> option
for other search methods.

=item C<< -block=integer] >>

Specifies the number of query sequences per search program
execution. This option is important when the searches for individual
queries are fast. By putting multiple query sequences into one block,
the execution performance is improved because search program setup and
shutdown costs per query. are reduced. The default is one query per
execution, i.e. no blocking.

=item C<< -[no]debug >>

Turns on additional debugging output.

=item C<< -minlength=integer >>

Specifies the minumum sequence length for a query to have. If a query
sequence is shorter than this parameter, the query is marked as
searched, but no results are generated. The minimum value for this
parameter is 6.

=item C<< -wordsize=integer >>

Specifies the word size for the Blast search methods. If not
specified, then the default value is used.

=item C<< -[no]selfself >>

When the annotator is used to compare a genome against itself, there's
usually no point in including the sequences that match against
themselves. This option, if specified, will result in self self
comparisons being included.

=item C<< -nice=string >>

The annotator will typically try to lower the priority of the
annotation processes it spawns. This option specifies the command to
lower the priority. The default setting is "nice -19".

=item C<< -fasty=string >>

Specifies the path to the C<fasty> program. Otherwise, the PATH environment variable is used to find "fasty35".

=item C<< -ggsearch=string >>


Specifies the path to the C<ggsearch> program. Otherwise, the PATH environment variable is used to find "ggsearch35".

=item C<< -exonerate=string >>

Specifies the path to the C<exonerate> program. Otherwise, the PATH environment variable is used to find "exonerate".

=back

=head1 AUTHOR

Robert Bruccoleri <bruc@congen.com>
Congenomics, LLC.

=cut
