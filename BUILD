#!/bin/sh -x

prefix=${BUILD_PREFIX:-/usr/local}

perl=${PERL_FOR_BUILD:-/usr/bin/perl}

source $prefix/setup.sh

export PGUSER=bruc
export PGHOST=localhost
export PGPORT=6543

make clean
set -e
# dropdb frescobi_test
$perl Makefile.PL PREFIX=${PERL_PREFIX:-/BAD/PERL/PREFIX/NEEDS/FIXING}
make
make test
make install

