#!perl

# Filter to implement identity filtering.

use warnings;
use strict;

use Bio::Frescobi::Genutil;

if ($#ARGV != 3) {
    die "Usage: $0 input-seq-file input-qual-file output-seq-file output-qual-file\n";
}

my ($dir, $file) = ($ARGV[0] =~ m%^(.*)/([^\/])+$%);
chdir($dir) || die "Unable to chdir $dir\n";
system_with_check("cp $ARGV[0] $ARGV[2]");
system_with_check("cp $ARGV[1] $ARGV[3]");
