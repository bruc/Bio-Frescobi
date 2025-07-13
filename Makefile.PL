use 5.008;
use ExtUtils::MakeMaker;
use strict;

my $prefix = $ENV{BUILD_PREFIX} || "/usr/local";

WriteMakefile(
	      NAME              => 'Bio::Frescobi',
	      VERSION_FROM      => 'lib/Bio/Frescobi/Frescobi.pm', # finds $VERSION
	      PREREQ_PM         => { 'CGI' => 3.0,
				     'Carp' => 1.0,
				     'Digest::MD5' => 2.36,
				     'Text::LevenshteinXS' => 0.03,
				     'SHA' => 2.01,
				     'DBI' => 1.646,
				     'DBD::Pg' => 3.18.0,
				     'DBD::SQLite' => 0,
				     'BioPerl' => 1.7.8,
				     'Data::Dumper' => 0,
				     'Exporter' => 0,
				     'File::Basename' => 0,
				     'Getopt::Long' => 0,
				     'POSIX' => 0,
				     'SHA' => 0,
				     'Sys::Hostname' => 0,
				     'Text::Wrap' => 0,
	      }, # e.g., Module::Name => 1.1
	      EXE_FILES         => [ 'scripts/annotator.pl',
				     'scripts/clean_date.pl',
				     'scripts/count_aa.pl',
				     'scripts/count_fastq.pl',
				     'scripts/create.sql',
				     'scripts/create_sqlite.sql',
				     'scripts/extract_seq.pl',
				     'scripts/fastq_sample.pl',
				     'scripts/identity.pl',
				     'scripts/loadseq.pl',
				     'scripts/lock.pl',
				     'scripts/pematch.pl',
				     'scripts/revcomp_fastq.pl',
				     'scripts/split_fasta_by_size.pl',
				     'scripts/update_filtered_seqs.pl',
				     'scripts/update_ns.pl',
				     'scripts/update_libraries.pl' ],
              INST_SCRIPT       => $prefix . '/frescobi/scripts',
              MAN1PODS          => { 'scripts/pematch.pl' => 'blib/man1/pematch.pl.1',
                                     'scripts/annotator.pl' => 'blib/man1/annotator.pl.1',
				     'scripts/loadseq.pl' => 'blib/man1/loadseq.pl.1',
				     'scripts/count_fastq.pl' => 'blib/man1/count_fastq.pl.1',
				     'scripts/revcomp_fastq.pl' => 'blib/man1/revcomp_fastq.pl.1',
				     'scripts/extract_seq.pl' => 'blib/man1/extract_seq.pl.1'
                                    },
              MAN3PODS            => { 'lib/Bio/Frescobi/Genutil.pm' => 'blib/man3/Bio::Frescobi::Genutil.3pm',
				      'lib/Bio/Frescobi/BigQuery.pm' => 'blib/man3/Bio::Frescobi::BigQuery.3pm', 
				      'lib/Bio/Frescobi/Sequtil.pm' => 'blib/man3/Bio::Frescobi::Sequtil.3pm',
				      'lib/Bio/Frescobi/CgPg.pm' => 'blib/man3/Bio::Frescobi::CgPg.3pm',
				      'lib/Bio/Frescobi/PGLoad.pm' => 'blib/man3/Bio::Frescobi::PGLoad.3pm',
				      'lib/Bio/Frescobi/PGLock.pm' => 'blib/man3/Bio::Frescobi::PGLock.3pm',
				      'lib/Bio/Frescobi/Dbseq.pm' => 'blib/man3/Bio::Frescobi::Dbseq.3pm',
				      'lib/Bio/Frescobi/PEMatch.pm' => 'blib/man3/Bio::Frescobi::PEMatch.3pm'
				     },
	      dist              => { COMPRESS => 'gzip -1' },
	      ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
	       (ABSTRACT_FROM  => 'lib/Bio/Frescobi/Genutil.pm', # retrieve abstract from module
		AUTHOR         => 'Robert Bruccoleri <bruc@acm.org>') : ()),
	      );

