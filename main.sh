#!/bin/bash
#expected file organisation:
# /project
#    `----> gencode_name (.gtf)
#    `----> herv_cat_file (.gff3)
#    `----> hERVd_name (dir)
#      `--> all_other_files (.gff3)
#    `----> TCGA_data
#      `--> genome_data (dir)
#        `-> file_id (.bam)
#      `--> transcriptome_data (dir)
#        `-> file_id (.bam)
# git remote set-url origin https://@github.com/<username>/<repositoryname>.gi
#git config credential.helper cache

REF_ANNOTATION_LOC=\
'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz'
REF_GENOME_LOC='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.p13.genome.fa.gz'
hERV_DIR="hERV_misc"
hERV_FILE="hERVd.gff3"
COMB_ANNOTATION="generv.gff3"
ANALYSIS_STEP="all"
PLATFORM="replit" #cluster-mgh 
PAIR_FILE="fq_pairs.csv"
SOURCE_LOC="../test_folder"
CLEAR_TMP="false"
OVER_WRITE="false"
COUNT_METHOD="seperated"
BATCH_NAME="none"

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-REF_ANNOTATION) REF_ANNOTATION_LOC="$2"; shift ;;
		-REF_GENOME) REF_GENOME_LOC="$2"; shift ;;
		-hERV_DIR) hERV_DIR="$2" shift ;;
		-hERV_FILE) hERV_FILE="$2" shift ;;
		-COMB_ANNOTATION) COMB_ANNOTATION="$2" shift ;;
		-SOURCE) SOURCE_LOC="$2" shift ;;
		-ANALYSIS_STEP) ANALYSIS_STEP="$2" shift ;; #steps: index, convert, qc, align, count
		-PLATFORM) PLATFORM="$2" shift ;;
		-COUNT_METHOD) COUNT_METHOD="$2" shift ;;
		-CLEAR_TMP) CLEAR_TMP="true" ;;
		-OVER_WRITE) OVER_WRITE="true" ;;
		-BATCH_NAME) BATCH_NAME="$2" ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done

main_loc="$PWD"

if [ ! -d "../hERV_Work" ]; then 
	mkdir ../hERV_Work 
fi 

cd ../hERV_Work

read -p "starting, enter programme run type: " RUN_TYPE
# RUN_TYPE="-d" #debug
if [ "$RUN_TYPE" = "-d" ]; then #have to keep spaces between square brackets and vars
	read -p "enter debug type: " DEBUG_TYPE
	# DEBUG_TYPE='-f' #function
	if [ "$DEBUG_TYPE" = '-f' ]; then
		read -p "enter function to debug: " FUNC
		${FUNC[0]} ${FUNC[@]:1}
		echo $?
	fi

	if [ "$DEBUG_TYPE" = '-u' ]; then
		read -p "enter file to debug: " FILE
		./$FILE ${FILE[@]:1}
		echo $?
	fi
elif [ "$RUN_TYPE" = "-r" ]; then
	SEQ_TYPE=("RNA-Seq" "WXS")
	( . $main_loc/setup.sh )
	( . $main_loc/downloads.sh )
	if [ $SOURCE_LOC = "tcga" ] && [ ! -d "tcga" ]; then 
		(	. $main_loc/tcga.sh )
	fi
	if [ ! $BATCH_NAME = "none" ]; then 
		(. $main_loc/test_batch.sh)
	else 
		(. $main_loc/analysis.sh)
	fi 
fi
# jq --arg rty "$RUN_TYPE" \
# 	 --arg dty "$DEBUG_TYPE" \
# 	 --arg plf "$PLATFORM" \
# 	 --arg 
# 	 '{"rty": $rty, "dty": $dty, "plf": $plf}'
#{"filters":{"op":"=","content":{"field":"cases.demographic.gender","value":["female"]}},"size":1}
