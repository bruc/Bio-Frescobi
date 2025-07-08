# This package contains the version number for the Frescobi suite
# and common functions for Frescobi.

package Bio::Frescobi::Frescobi;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Carp;
use CGI qw(:all);
use Bio::Frescobi::Genutil;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Dbseq;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Config;

@ISA = ('Exporter');
@EXPORT = qw(&update_duplicates
	     );

$VERSION = '0.33';

1;

sub update_duplicates {
    my ($pg, $base_table) = @_;

    my $duplicate_table = "duplicate_" . $base_table;
    my ($i, $name, $seqid, $cmd);
    my (%count_seqids);
    my (%seqid_for_name, %names_for_seqid);
    my (%duplicate_seqid_for_name, %duplicate_names_for_seqid);

    foreach my $rowp (@{$pg->query("select name, seqid from $duplicate_table")}) {
	($name, $seqid) = @{$rowp};
	$duplicate_seqid_for_name{$name} = $seqid;
	$duplicate_names_for_seqid{$seqid}{$name} = 1;
    } 
    my @data = @{$pg->query("select name, seqid from $base_table")};
    foreach my $rowp (@data) {
	($name, $seqid) = @{$rowp};
	if (not defined $count_seqids{$seqid}) {
	    $count_seqids{$seqid} = 0;
	}
	$count_seqids{$seqid} += 1;
    }
    foreach my $rowp (@data) {
	($name, $seqid) = @{$rowp};
	if ($count_seqids{$seqid} > 1) {
	    $seqid_for_name{$name} = $seqid;
	    $names_for_seqid{$seqid}{$name} = 1;
	}
    }
    $pg->autocommit(0);
    foreach $seqid (keys %names_for_seqid) {
	if (scalar(keys %{$names_for_seqid{$seqid}}) == 1) {
	    if (exists($duplicate_names_for_seqid{$seqid})) {
		$cmd = "delete from $duplicate_table where seqid = '$seqid'";
		$pg->command($cmd);
	    }
	}
	elsif (scalar(keys %{$names_for_seqid{$seqid}}) > 1) {
	    if (not exists($duplicate_names_for_seqid{$seqid}) or
		scalar(keys %{$duplicate_names_for_seqid{$seqid}}) == 1) { }
	    else {
		foreach $name (keys %{$names_for_seqid{$seqid}}) {
		    if (exists($duplicate_names_for_seqid{$seqid}{$name})) {
			delete $duplicate_names_for_seqid{$seqid}{$name};
			delete $names_for_seqid{$seqid}{$name};
		    }
		}
		foreach $name (keys %{$duplicate_names_for_seqid{$seqid}}) {
		    $cmd = "delete from $duplicate_table where seqid = '$seqid' " .
			" and name = " . &quotify($name);
		    $pg->command($cmd);
		}
	    }
	    foreach $name (keys %{$names_for_seqid{$seqid}}) {
		$cmd = "insert into $duplicate_table (seqid, name) values ('$seqid'," .
		    &quotify($name) . ")";
		$pg->command($cmd);
	    }
	}
	else {
	    die 'Unexpected value for scalar(keys %{$names_for_seqid{' . $seqid . '} = ' .
		scalar(keys %{$names_for_seqid{$seqid}}) . "\n";
	}
    }
    foreach $seqid (keys %duplicate_names_for_seqid) {
	if (not exists($names_for_seqid{$seqid})) {
	    $cmd = "delete from $duplicate_table where seqid = '$seqid'";
	    $pg->command($cmd);
	}
    }
    $pg->autocommit(1);
}

