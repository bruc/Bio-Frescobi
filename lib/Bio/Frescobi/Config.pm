#
# This file contains global variables used for configuring Frescobi's tools.

package Bio::Frescobi::Config;
use strict;
use vars  qw(@ISA @EXPORT $VERSION
	     $tmpdir );
use Exporter;
use File::Basename;
use Carp;

@ISA = ('Exporter');
@EXPORT = qw($tmpdir);

$VERSION = '0.02';

# $tmpdir holds the temporary directory used by the tools

if (exists($ENV{"TMPDIR"}) and $ENV{"TMPDIR"}) {
    $tmpdir = $ENV{"TMPDIR"};
}
else {
    $tmpdir = "/tmp";
}

1;
