from Bio import Seq, SeqIO

import local_env_variables.env_variables as env

seq_file = env.orthoDB_files.all_seqs_fasta
sqlite_file = env.orthoDB_files.all_seqs_sqlite
records = SeqIO.index_db(sqlite_file, seq_file, 'fasta')
records.close()
