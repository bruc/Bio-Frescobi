#!/bin/sh

# This test must run with the t directory being the default.

driver=Pg
if [[ "$1" != "" ]]
then
    driver=$1
fi

set -e
function execute_sql ()
{
    local sql="$1"
    local db="$2"

    if [[ "$db" == "" ]]
    then
	db=$dbname
    fi

    if [[ "$driver" == "Pg" ]]
    then
	psql -c "$1" $db
    else
	sqlite3 -separator " | " -header -cmd "$1" $db < /dev/null
    fi
}

source /pg/postgresql-16.4/profile
source $BUILD_PREFIX/binaries/perl_modules/setup.sh
export PGUSER=bruc
export TMPDIR=$(pwd)/tmp
export PAGER=/bin/cat

if [[ "$driver" == "SQLite" ]]
then
    dbname=testest.sqlite
elif [[ "$driver" == "Pg" ]]
then
    dbname=testest
else
    echo Unknown driver = $driver
    exit 1
fi

cp /dev/null load.test.${driver}.log
cp /dev/null filter.test.${driver}.log
cp /dev/null compare.test.${driver}.log
cp /dev/null annotator.test.${driver}.log
cp /dev/null seq.test.${driver}.log
cp /dev/null formatdb.test.${driver}.log

cd ..
scripts=$(pwd)/blib/script
cd t
export PATH=$scripts:$PATH
export PERL5LIB=../blib/lib/${PERL5LIB:+:${PERL5LIB}}

if [[ "$driver" == "Pg" ]]
then
    dropdb $dbname
    createdb $dbname
    psql $dbname <../scripts/create.sql
else
    rm -f $dbname
    sqlite3 $dbname <../scripts/create_sqlite.sql
fi

loadseq.pl -prefix=ld1. -dbname=$dbname -driver=$driver -library=danio -repository=congen dr.seq.1.fasta{,.qual}.triple |& \
    clean_date.pl >>load.test.${driver}.log
execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name"
psql -c "truncate raw_seqs" testest
psql -c "drop sequence seqid_serial;\
         CREATE SEQUENCE seqid_serial \
                start 1 increment 1 maxvalue 2147483647 minvalue 1  cache 1 ;" testest

loadseq.pl -prefix=ld1. -dbname=testest -library=danio -repository=congen dr.seq.1.fasta{,.qual}.triple |& \
    clean_date.pl >>load.test.${driver}.log

execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name" 
execute_sql "delete from raw_seqs"

if [[ "$driver" == "Pg" ]]
then
    dropdb $dbname
    createdb $dbname
    
    psql $dbname <../scripts/create.sql
else
    rm -f $dbname
    sqlite3 $dbname <../scripts/create_sqlite.sql
fi    

execute_sql "
delete from dbs where common_db_name = 'yeast';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('yeast',
	       'testdb/yeast_nrpep',
	       'May 23 18:57 2009',
	       'protein');

delete from dbs where common_db_name = 'HCDNA';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('HCDNA',
	       'testdb/human_cdna.fasta',
	       'Apr 10 17:05 2009',
	       'nucleic');

delete from dbs where common_db_name = 'nrpshmm';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('nrpshmm',
	       'NRPS_HMM',
	       'Aug 21 16:57 2013',
	       'protein');

delete from dbs where common_db_name = 'exonerate_db';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('exonerate_db',
	       'exonerate_db',
	       'Aug 22 19:06 2013',
	       'nucleic');

delete from dbs where common_db_name = 'exonerate_prot_db';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('exonerate_prot_db',
	       'exonerate_prot_db',
	       'Aug 22 19:26 2013',
	       'protein');
"

loadseq.pl -prefix=ld1. -driver=$driver -dbname=$dbname -library=danio -repository=congen dr.seq.1.fasta{,.qual}.single |& clean_date.pl >>load.test.${driver}.log
execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name"
execute_sql "delete from raw_seqs" 
if [[ "$driver" == "Pg" ]]
then
    execute_sql "drop sequence seqid_serial;\
        	 CREATE SEQUENCE seqid_serial \
                 start 1 increment 1 maxvalue 2147483647 minvalue 1  cache 1 ;"
else
    execute_sql "update seqid_serial set i = 0"
fi
loadseq.pl -prefix=ld1. -driver=$driver -dbname=$dbname -library=danio -repository=congen dr.seq.1.fasta{,.qual}.single |& clean_date.pl >>load.test.${driver}.log
execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name" 

rm -rf $TMPDIR/test/filter

update_filtered_seqs.pl -tmp=$TMPDIR/test/filter -driver=$driver -dbname=$dbname -maxproc=1 ${scripts}/identity.pl |& clean_date.pl >>filter.test.${driver}.log

execute_sql "select f.name, f.seqid, d.checksum from filtered_seqs f, seq_data d where f.seqid = d.seqid order by f.name" 

loadseq.pl -prefix=ld2. -driver=$driver -dbname=$dbname -library=danio -repository=congen dr.seq.2.fasta{,.qual} |& clean_date.pl >>load.test.${driver}.log
execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name" $dbname
update_filtered_seqs.pl -tmp=$TMPDIR/test/filter -library=danio -redo -driver=$driver -dbname=$dbname -maxproc=1 ${scripts}/identity.pl |& clean_date.pl >>filter.test.${driver}.log
execute_sql "select f.name, f.seqid, d.checksum, d.length from filtered_seqs f, seq_data d where f.seqid = d.seqid order by f.name" $dbname

loadseq.pl -prefix=ld3. -driver=$driver -dbname=$dbname -library=danio -repository=congen dr.seq.3.fasta{,.qual} |& clean_date.pl >>load.test.${driver}.log
loadseq.pl -prefix=ld4. -driver=$driver -dbname=$dbname -library=danio -repository=congen dr.seq.4.fasta{,.qual} |& clean_date.pl >>load.test.${driver}.log
execute_sql "select r.name, r.seqid, d.checksum from raw_seqs r, seq_data d where r.seqid = d.seqid order by r.name" $dbname

update_filtered_seqs.pl -tmp=$TMPDIR/test/filter -driver=$driver  -dbname=$dbname -maxproc=1 ${scripts}/identity.pl |& clean_date.pl >>filter.test.${driver}.log
execute_sql "select f.name, f.seqid, d.checksum from filtered_seqs f, seq_data d where f.seqid = d.seqid order by f.name" $dbname

rm -rf $TMPDIR/test/blast

loadseq.pl -prefix=p. -seqtype=protein -driver=$driver -dbname=$dbname \
    -library=s_hygro -repository=congen -uniquify -uniquify_separator=_ NC_017765.prot |& \
    clean_date.pl >>load.test.${driver}.log

annotator.pl -dbname=$dbname \
	     -driver=$driver \
             -db=nrpshmm \
             -maxhits=1000 \
             -method=hmm \
             -maxproc=3 \
	     -block=100 \
             -sql="select r.seqid from raw_seqs r where r.library = 's_hygro'" \
             -nofilter \
             -cutoff=1000 \
	     |& clean_date.pl >>annotator.test.${driver}.log
		 

execute_sql "select hits.seqid, hits.common_db_name, hits.method, hits.key, hits.description, hits.best_score, hsps.score, hsps.bits, hsps.expect, hsps.length, hsps.gaps, hsps.num_identical, hsps.num_conserved, hsps.frac_identical, hsps.frac_conserved, hsps.query_strand_positive, hsps.sbjct_strand_positive, hsps.query_start, hsps.query_end, hsps.sbjct_start, hsps.sbjct_end from hits, hsps where hits.hit_id = hsps.hit_id and hits.seqid like 'P%'" | sort

(extract_seq.pl -dbname=$dbname -driver=$driver -sequences -sql="select seqid from seq_data" -type=nucleic >test.seq) >>  seq.test.${driver}.log 2>&1 
(extract_seq.pl -dbname=$dbname -driver=$driver -quality -nosequences -sql="select seqid from seq_data" -type=nucleic >test.seq.qual)  >> seq.test.${driver}.log 2>&1
(extract_seq.pl -dbname=$dbname -driver=$driver -sequences -sql="select seqid from seq_data" -type=protein >test.prot.seq) >>  seq.test.${driver}.log 2>&1 
    
execute_sql "select hits.common_db_name, hits.method, hits.key, hits.description, hits.best_score, hsps.score, hsps.bits, hsps.expect, hsps.length, hsps.gaps, hsps.num_identical, hsps.num_conserved, hsps.frac_identical, hsps.frac_conserved, hsps.query_strand_positive, hsps.sbjct_strand_positive, hsps.query_start, hsps.query_end, hsps.sbjct_start, hsps.sbjct_end, hits.seqid from hits, hsps where hits.hit_id = hsps.hit_id " $dbname | sort
execute_sql "select hits.common_db_name, hits.method, hits.key, hits.description, hits.best_score, hsps.score, hsps.bits, hsps.expect, hsps.length, hsps.gaps, hsps.num_identical, hsps.num_conserved, hsps.frac_identical, hsps.frac_conserved, hsps.query_strand_positive, hsps.sbjct_strand_positive, hsps.query_start, hsps.query_end, hsps.sbjct_start, hsps.sbjct_end from hits, hsps where hits.hit_id = hsps.hit_id " $dbname | sort
execute_sql "select count(*) from hits, hsps where hits.hit_id = hsps.hit_id group by seqid order by count(*)" $dbname
execute_sql "select count(*) from hits, hsps where hits.hit_id = hsps.hit_id group by seqid, hits.hit_id order by count(*) " $dbname
      
(extract_seq.pl \
     -dbname=$dbname \
     -driver=$driver \
     -sql="select seqid from raw_seqs where library = 'danio'" \
     -type=nucleic \
     >danio.seq ) >> seq.test.${driver}.log 2>&1 
danio_path=$(pwd)

formatdb -i danio.seq -p F -o T >> formatdb.test.${driver}.log 2>&1
formatdb -i testdb/yeast_nrpep -p T -o T >> formatdb.test.${driver}.log 2>&1

execute_sql "delete from dbs where common_db_name = 'danio';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('danio',
	       '$danio_path/danio.seq',
	       'Oct 18 18:57 2016',
	       'nucleic');
"

annotator.pl \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=danio \
    -noselfself \
    -method=blastn \
    -cutoff=0.001 \
    -sql="select seqid from raw_seqs where library = 'danio'" |& clean_date.pl >>annotator.test.${driver}.log
execute_sql "select count(*) from hits where seqid = key and common_db_name = 'danio'" $dbname
annotator.pl \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=danio \
    -selfself \
    -method=blastn \
    -cutoff=0.001 \
    -replace \
    -sql="select seqid from raw_seqs where library = 'danio'" |& clean_date.pl >>annotator.test.${driver}.log

execute_sql "select count(*) from hits where seqid = key and common_db_name = 'danio'" $dbname

execute_sql "
delete from hits;
delete from hsps;
delete from annotation_runs;
"
rm -f testdb/human_cdna.fasta.n??
formatdb -i testdb/human_cdna.fasta -p F -o T

annotator.pl \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=HCDNA \
    -method=blastn \
    -cutoff=1.0e-5 \
    -sql="select seqid from raw_seqs where seqid in (select seqid from seq_data where seq_type = 'nucleic')" |& clean_date.pl >>annotator.test.${driver}.log
execute_sql "select seqid, best_score, key, description from hits order by seqid, best_score, key" $dbname

annotator.pl \
    -replace \
    -extraopts="-r 2 -q -10" \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=HCDNA \
    -method=blastn \
    -cutoff=1.0e-5 \
    -sql="select seqid from raw_seqs where seqid in (select seqid from seq_data where seq_type = 'nucleic')" |& clean_date.pl >>annotator.test.${driver}.log

execute_sql "select seqid, best_score, key, description from hits order by seqid, best_score, key" $dbname

annotator.pl \
    -replace \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=HCDNA \
    -method=blastn \
    -cutoff=1.0e-10 \
    -sql="select seqid from raw_seqs where seqid in (select seqid from seq_data where seq_type = 'nucleic')" |& clean_date.pl >>annotator.test.${driver}.log
execute_sql "select seqid, best_score, key, description from hits where method='blastn' order by seqid, best_score, key" $dbname

#     -cutoff=1.0e-5 \
#     -extraopts="-r 2 -q -10" \

rm testdb/human_cdna.fasta.n??
makeblastdb -in testdb/human_cdna.fasta -dbtype nucl  >> formatdb.test.${driver}.log 2>&1 

annotator.pl \
    -replace \
    -maxproc=4 \
    -block=100 \
    -dbname=$dbname \
    -driver=$driver \
    -db=HCDNA \
    -method=blast+blastn \
    -cutoff=1.0e-10 \
    -extraopts="-penalty -3 -task blastn" \
    -sql="select seqid from raw_seqs where seqid in (select seqid from seq_data where seq_type = 'nucleic')" |& clean_date.pl >>annotator.test.${driver}.log
execute_sql "select seqid, best_score, key, description from hits where method='blast+blastn' order by seqid, best_score, key" $dbname

#     -extraopts="-reward 2 -penalty -10" \
    #     -cutoff=1.0e-5 \
    
exit 0
