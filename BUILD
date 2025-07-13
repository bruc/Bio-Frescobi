#!/bin/sh -x

# Environment variables used here:

# BUILD_PREFIX:   top level installation location. It is assumed that there is
# a setup.sh script to define the environment. This is optional.

# PERL_FOR_BUILD: The executable for running Perl. This can be used if the system
# has a custom Perl executable.

# PERL_PREFIX: Top level for Perl module libraries. The following
# statement for adding additional Perl module libraries is an example
# usage for the PERL_PREFIX environment variable:

# export PERL5LIB=$PERL_PREFIX/share/perl5:$PERL_PREFIX/lib64/perl5:$PERL_PREFIX/lib/perl5${PERL5LIB:+:${PERL5LIB}}

prefix=${BUILD_PREFIX:-/usr/local}

perl=${PERL_FOR_BUILD:-/usr/bin/perl}

source $prefix/setup.sh

# The following is an example of setting up environment variables for
# running PostgreSQL used in the testing code. The user must have the
# roles for creating databases.

export PGUSER=bruc
export PGHOST=localhost
export PGPORT=6543

make clean
set -e
$perl Makefile.PL PREFIX=${PERL_PREFIX:-/BAD/PERL/PREFIX/NEEDS/FIXING}
make
export PATH=$prefix/frescobi/scripts:$PATH
make test
make install

