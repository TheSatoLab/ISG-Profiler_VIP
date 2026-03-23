#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -l lmem
#$ -l s_vmem=60G,mem_req=60G
#$ -o /home/kyokok/projects/2507blastx/log/260302
#$ -e /home/kyokok/projects/2507blastx/log/260302

# change directory
base_path="/home/kyokok/projects/2507blastx"
cd ${base_path}

# commands
module load /usr/local/package/modulefiles/blast+/2.15.0

blastn -query ${base_path}/data/260302_virus_hit/virus_hit_all_contigs_w_bl.fasta \
	-db /usr/local/db/blast/ncbi/core_nt \
	-out ${base_path}/output/260302_blastn/virus_hit_all_contigs_w_bl.txt \
	-evalue 1e-4 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe staxids stitle qlen slen qcovs qcovhsp" \
	-num_threads 8
