#!/bin/bash

### the script succeeds from test_single.lsf

#SBATCH --partition=normal
#SBATCH --job-name=caki_slurm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=72G
#SBATCH --output=~/scratch/erv/herv_project/logs/caki_slurm.log
#SBATCH --error=~/scratch/erv/herv_project/logs/caki_slurm.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jliu78@mgh.harvard.edu

cd ~/scratch/erv/herv_project
if [ ! -d "logs" ]; then
	mkdir logs
fi


########### NOTE #############################################
### this should be good for both human data and mouse data, test_single_mouse.lsf can be reference to show only differences
###  are parameters "-REF_GENOME", "-REF_TRANSCRIPT" and "TRANSCRIPTS" of "bash main.sh"

########### FILE INPUT ########################################
### raw_data stores all fastq files named as 'fq' or better 'fq.gz' to be processed, it is suggested by siyi to store gzipped fastq files in the raw_data folder. Please note $name=random_name.fq.gz; ${name#*.} is used to get the postfix of the raw fastq file names, so DO NOT USE random_name.r1.fq.gz as raw fastq file name, but instead use random_name_r1.fq.gz. See below example
### (base)$ name=random_name.fq.gz; echo ${name#*.}
###   fq.gz
### (base)$ name=random_nam.e.fq.gz; echo ${name#*.}
###   e.fq.gz

### if raw_data folder has reletive dir as '../raw_data', input value as "raw_data" to the parameter '-SOURCE' of main.sh

### for subread, store in ../hERV_Work the annotation gtf/gff3 files of both transcriptom and ERV, genome fasta 

### for salmon, if you have index already, store in ../hERV_Work/salmon the index folder named as "$hervname_index", by default it is "package-entities-erv_index"
###             if you need generate index, 
###                         first option is to supply parameter "-REF_GENOME", "-REF_TRANSCRIPT" and "TRANSCRIPTS" of "bash main.sh" with customized link, for example:
###                                            bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
###                                                         -REF_GENOME "https://.*/weird.genome.fa.gz" \
###                                                         -REF_TRANSCRIPT "https://.*/transcripts_of_widely_Used_Normal_OR_Commen_Genes.fa.gz" \
###                                                         -TRANSCRIPTS "https://.*/transcripts_of_erv_or_other_special_features.fa.gz" \
###                                                         -TRANSCRIPTS "https://.*/if_you_have_another_special_features.fa.gz"

###                         2nd option is to store the already downloaded "fa.gz" or "fa" files in ../hERV_Work/ to avoid downloading, their file names need follow the below: 
###                          1. genome fasta can be named as 'GRCh38.p13.genome.fa' which matchs with
###                              " REF_GENOME_LOC='https://.*/GRCh38.p13.genome.fa.gz' " in main.sh;
###                              or named as something else but make sure you change the value of the "REF_GENOME_LOC" variable in main.sh
###                                                              by supplying parameter "-REF_GENOME" of "bash main.sh" with corresponding value, for example:
###                                            bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
###                                                         -REF_GENOME "https://.*/weird.genome.fa.gz"
###                                            to match you genemo fa file you name as "weird.genome.fa"
###                                            !!!siyi got the fasta name by remove everything before the last "/" and the last "/" itself and remove .gz$

###                          2. similarly, transcriptom fasta can be named as 'gencode.v41.transcripts.fa.gz' which matchs with "REF_TRANSCRIPT_LOC='https://.*/gencode.v41.transcripts.fa.gz'" in main.sh;
###                              or named as something else but make sure you change the value of the "REF_TRANSCRIPT_LOC" variable in main.sh
###                                                              by supplying parameter "-REF_TRANSCRIPT" of "bash main.sh" with corresponding value, for example:
###                                            bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
###                                                         -REF_TRANSCRIPT "/weird.transcripts.fa.gz"
###                                            to match you transcriptom fq file you name as "weird.transcripts.fa"

###                          3. similarly, special fasta can be named as 'package-entities-erv.fa.gz' which matchs with "TRANSCRIPT_LOCS='https://.*/package-entities-erv.fa.gz'" in main.sh;
###                              or named as something else but make sure you change the value of the "TRANSCRIPT_LOCS" variable in main.sh
###                                                              by supplying parameter "-TRANSCRIPTS" of "bash main.sh" with corresponding value, for example:
###                                            bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
###                                                         -TRANSCRIPTS "https://.*/weird.transcripts.fa.gz"
###                                            to match you special fq.gz file you name as "weird.transcripts.fa.gz"

###                          3a. if you input multiple special fasta, you have to supply multiple times of parameter "-TRANSCRIPTS" to "bash main.sh", for example:
###                                            bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
###                                                         -TRANSCRIPTS "https://.*/weird.transcripts.fa.gz"
###                                                         -TRANSCRIPTS "https://.*/2nd.transcripts.fa.gz"
###                                                         -TRANSCRIPTS "https://.*/3rd.fa.gz"
###                                            and you can store already downloaded "weird.transcripts.fa", 2nd.transcripts.fa.gz and 3rd.fa.gz in ../hERV_Work
###                                            the codes with first check if there are files corresponding to file names obtained from the value of TRANSCRIPT_LOCS,
###                                            if yes, use them, if no download using the link

### The script will identify files as paired when there are two files named in the way that
###   there is file2 whose name is one of sub(name of file1, pattern="r1", replacement="r2") and sub(name of file1, pattern="R1", replacement="R2")
### see group_fastq() and check_name() in analysis.sh for detail

### parameter 'EXIT_ON_SINGLE' of main.sh has default value of 'false', thus single-end sequenced fq files can be analysis using the same bash command

### for mgh cluster, keep -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1

bash main.sh -PLATFORM "cluster-mgh" -ANALYSIS_STEP "all" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 -REF_GENOME "https://a/GRCm39.primary_assembly.genome.fa" \
-TRANSCRIPTS "https://a/mouse_erv_combined_repeatMasker_on_m39.fa" \
-TRANSCRIPTS "transgene.fa" \
-REF_TRANSCRIPT "https://a/mouse_gencode.vM33.transcripts.fa.gz"

### for local desktop like the imagedeskto, either not supply '-PLATFORM' parameter or supply anything else, for example
### bash main.sh -PLATFORM "anything_different_with_cluster-mgh"

################## OUTPUT FILE ######################################
### all output files will be stored in ../hervwork/results

