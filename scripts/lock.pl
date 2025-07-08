#!perl

use strict;
use warnings;


use Bio::Frescobi::PGLock;
use Getopt::Long;


my $dbname="";
my $table="";
my $wait=10;
my $cycles=20000;
my $verbose = 1;

my $usage = <<EOF;
lock.pl -dbname=<dbname> -table=<table> <lock|unlock>
        [-driver=Pg|SQLite ]
        [-wait=<integer>   ]     Default: $wait seconds
        [-cycles=<integer> ]     Default: $cycles cycles
        [-[no]verbose      ]     Default: $verbose
EOF
    ;

GetOptions("dbname=s" => \$dbname,
	   "driver=s" => \$driver,
           "table=s" => \$table,
           "wait=i" => \$wait,
           "cycles=i" => \$cycles,
           "verbose!" => \$verbose);

if (not $dbname or
    not $table or
    scalar(@ARGV) != 1 or
    $ARGV[0] !~ m/^(un)?lock$/ ) {
    print STDERR $usage;
    exit(1);
}

my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => 0,
				  default_echo => $verbose);

Bio::Frescobi::PGLock::lock_wait_time($wait);
Bio::Frescobi::PGLock::lock_wait_cycles($cycles);      # About 3 days.

if ($ARGV[0] eq 'lock') {
    &Bio::Frescobi::PGLock::grab_lock($pg, $table);
}
elsif ($ARGV[0] eq 'unlock') {
    &Bio::Frescobi::PGLock::free_lock($pg, $table, 0);
}
else {
    die sprintf("Programming error: Unexpected lock command: %s\n",
                $ARGV[0]);
}
