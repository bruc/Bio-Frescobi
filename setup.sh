# Partial script to setup environment for Frescobi Note that the
# BUILD_PREFIX and PERL_PREFIX environment variables need to be set
# prior to this script as well as the correct PERL5LIB.

prefix=${BUILD_PREFIX:-/you/need/to/set/BUILD_PREFIX}
perl_perfix=${PERL_PREFIX:-/you/need/to/set/PERL_PREFIX}

export PATH=$prefix/binaries/ncbi-050604/bin:$prefix/binaries/blast-2.14.0+/bin:$prefix/frescobi/scripts:$PATH

if [[ "$MANPATH" == "" ]]
then
    export MANPATH=$(man -w)
fi

export MANPATH=$perl_prefix/share/man${MANPATH:+:${MANPATH}}



