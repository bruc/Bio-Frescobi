drop table if exists seqid_serial;

create table seqid_serial (i integer);
insert into seqid_serial (i) values (0);

drop table if exists hit_id_serial;

CREATE table  hit_id_serial (i integer);
insert into hit_id_serial (i) values (0);

drop table if exists hit_id_record;

create table hit_id_record (
       hit_id_val       integer,
       create_time      timestamp with time zone);

drop table if exists seq_data;

create table seq_data (
       seqid             text, 
       checksum          char(8),
       length            int4,
       create_time       timestamp with time zone,
       seq_type          text); 

create index seq_data_seqid on seq_data  (seqid);
create index seq_data_checksum on seq_data (checksum);

drop table if exists seq_pieces;

create table seq_pieces (
       seqid       text,
       segnum      int4,
       seq_piece   text,
       qual_piece  text);

create index seq_pieces_seqid on seq_pieces (seqid , segnum );

drop table if exists locks;

create table locks (
       table_name   text,
       pid          text,
       lock_time    timestamp with time zone,
       host         text);

drop table if exists raw_seqs;

create table raw_seqs (
       name       text not null,
       repository text default '',
       library    text default '',
       seqid      text,
       taxon      int4,
       filtered   bool,
       annotation text,
       create_time timestamp with time zone,
       name_type  text default '',
       concatid   text); 

create index raw_seqs_name on raw_seqs (name);
create index raw_seqs_seqid on raw_seqs (seqid);
create index raw_seqs_library on raw_seqs (library);
create index raw_seqs_id on raw_seqs (repository, library, name_type, name);
create index raw_seqs_concatid on raw_seqs (concatid);

drop table if exists libraries;

create table libraries (
       library      text,
       description  text,
       code         int4,
       seqcount     int4,
       assembly_count int4);

insert into libraries (library, description, code, seqcount, assembly_count)
       values('', 'Null library', 0, 0, 0);

drop table if exists duplicate_raw_seqs;

create table duplicate_raw_seqs ( -- Keep track of any duplicates sequences.
       name       text,
       seqid      text);

create index duplicate_raw_seqs_name on duplicate_raw_seqs (name);
create index duplicate_raw_seqs_seqid on duplicate_raw_seqs (seqid);

drop table if exists filtered_seqs;

create table filtered_seqs (
       name       text,           -- Link to raw_seqs table
       seqid      text,           -- Filtered sequence.
       compared   bool,
       create_time timestamp with time zone);

create index filtered_seqs_name on filtered_seqs (name);
create index filtered_seqs_seqid on filtered_seqs (seqid);
       
drop table if exists duplicate_filtered_seqs;

create table duplicate_filtered_seqs ( -- Keep track of any duplicates sequences.
       name       text,
       seqid      text);

create index duplicate_filtered_seqs_name on duplicate_filtered_seqs (name);
create index duplicate_filtered_seqs_seqid on duplicate_filtered_seqs (seqid);

drop table if exists hits;

create table hits (
       seqid           text, 
       common_db_name  text,
       method          text, 
       key             text, 
       description     text, 
       best_score      float8, 
       hit_id          text);

create index hits_seqid on hits (seqid);
create index hits_id on hits (hit_id);
create index hits_db_all on hits (seqid, common_db_name, method);
create index hits_db_key on hits (common_db_name, key);

drop table if exists hsps;

create table hsps (
       hit_id                  text,      
       score                   int, 
       bits                    float8, 
       expect                  float8, 
       length                  int, 
       gaps                    int, 
       num_identical           int, 
       num_conserved           int,
       frac_identical          float8,
       frac_conserved          float8,
       query_strand_positive   bool, 
       sbjct_strand_positive   bool, 
       query_start             int, 
       query_end               int, 
       sbjct_start             int, 
       sbjct_end               int);

create index hsps_id on hsps (hit_id);

drop table if exists annotation_runs;

create table annotation_runs (
       seqid                 text,
       common_db_name        text,
       method                text,
       start_time            timestamp with time zone,
       status                text);

create index annotation_runs_seqid on annotation_runs (seqid);
create index annotation_runs_db_seqid
    on annotation_runs
       (common_db_name,
        method,
        seqid);

drop table if exists annotator_run_count;

CREATE table annotator_run_count (i integer);
insert into annotator_run_count (i) values (0);

drop table if exists annotator_history;

create table annotator_history (
       library_number  int4,            -- linked to contig_run
       common_db_name  text,
       method          text,
       hit_id_value    int4);

create index annotator_history_db
    on annotator_history (library_number, common_db_name, method);

drop table if exists genscan;

create table genscan (
       seqid           text,
       subopt          bool,
       gene            int,
       exon            int,
       ex_type         text,
       plus_strand     bool,
       start           int,
       stop            int,
       frame           int,
       phase           int,
       i_ac            int,
       do_t            int,
       coding_score    int,
       prob            float,
       total_score     float);

create index genscan_seqid on genscan (seqid);

drop table if exists pred_pep;

create table pred_pep (
       n_seqid         text,
       seqid           text,  -- Predicted peptide sequence.
       method          text,
       gene            int);  -- Gene within the seqid. This is relative
                              -- to the gene list from the methods.

create index pred_pep_n_seqid on pred_pep (n_seqid);


drop table if exists dbs;

create table dbs (
       common_db_name        text,
       file_name             text,
       last_update_time      timestamp with time zone,
       data_type             text,
       internal              bool default 'f',
       generation_command    text,
       prepare_command       text);

insert into dbs
       (common_db_name, data_type, internal, generation_command, prepare_command)
       values('internal_filtered',
              'nucleic',
              't',
              'extract_seq.pl -table=^filtered_seqs',
              'formatdb -o T -p F -i');

drop table if exists table_of_comparison_tables;

create table table_of_comparison_tables (
       name            text,
       script          text,
       script_options  text,
       description     text);

drop table if exists comparison_table_elements;

create table comparison_table_elements (
       table_name      text,
       element_name    text,
       element_type    text,
       description     text,
       index_type      text,
       index_func      text);

drop table if exists blast_a_comparisons;

create table blast_a_comparisons
       (name1             text,
        name2             text,
        library_code      int4,
        bits              float8,
        score             float8,
        normalized_score  float8,
        expected          float8,
        matches           int4,
        match_length      int4,
        gaps              int4,
        fake              bool);
       

delete from table_of_comparison_tables where name = 'blast_a';
delete from comparison_table_elements where table_name = 'blast_a';

insert into table_of_comparison_tables
       (name, script, script_options, description)
       values ('blast_a',
               'blast_compare',
               '-min_overlap=30 -gap_penalty=50 -ext_penalty=10 -mismatch_penalty=-30 -match_reward=2 -word_size=0 -db_len=2000000000 ',
               'Blast with parameters to avoid homologies.');

insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'bits', 'float8', 'bit score', 'btree', 'float8_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'score', 'float8', 'blast score', 'btree', 'float8_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'normalized_score', 'float8', 'blast score / length * 100', 'btree', 'float8_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'expected', 'float8', 'Expectation of randomness', 'btree', 'float8_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'matches', 'int4', 'Number of nucleotide matches in overlap region', 'btree', 'int4_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'match_length', 'int4', 'Length of overlap region', 'btree', 'int4_ops');
insert into comparison_table_elements
       (table_name, element_name, element_type,
        description, index_type, index_func)
       values ('blast_a', 'gaps', 'int4', 'Number of gaps in overlap region', 'btree', 'int4_ops');

drop table if exists cluster_serial;

create table cluster_serial (i integer);
insert into cluster_serial (i) values (0);

drop table if exists contig_run;

create table contig_run (
       number           int4,           -- linked to cluster table.
       program          text,           -- assembly program.
       program_options  text,           -- options to assembly program.
       compare_table    text,           -- table from which comparison data came.
       criterion        text,           -- taken from assembler command.
       operator         text,
       cutoff           text,
       create_time      timestamp with time zone,
       where_clause     text,
       reference_cluster int4,
       library          text,
       completed        bool);          -- records whether the assembler run completed
                                        -- and was properly analyzed.   

insert into contig_run (number, program, create_time, completed)
       values(0, 'placeholder', datetime('now', 'localtime'), 'f');

drop table if exists users;

create table users (
       username    text,
       fullname    text,
       email       text);

insert into users (username, fullname, email)
       values ('bruc', 'Robert E. Bruccoleri', 'bruc@acm.org');

drop table if exists curation;

create table curation (
       seqid          text,
       create_time    timestamp with time zone,
       username       text,
       comment        text,
       current        bool,
       invalidation_time timestamp with time zone,
       version        integer);

create index curation_seqid on curation (seqid);
create index curation_user on curation (username);

drop table if exists library_groups;
create table library_groups (
       group_name     text,
       group_desc     text);

drop table if exists library_group_members;
create table library_group_members (
       group_name     text,
       library        text);
