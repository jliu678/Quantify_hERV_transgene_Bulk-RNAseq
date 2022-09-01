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

clear_tmp() {
	if [ $CLEAR_TMP = "true" ]; then 
		rm -r tmp
		rm $PAIR_FILE
	fi
}

REF_ANNOTATION_LOC=\
'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz'
REF_GENOME_LOC='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.p13.genome.fa.gz'
REF_TRANSCRIPT_LOC='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz'

TRANSCRIPT_LOCS=('https://herv.img.cas.cz/f/package-entities-erv.fa.gz')
ANNOT_DIR="hERV_misc"
ANNOT_LOCS=("https://herv.img.cas.cz/f/package-entities-erv.gff3.gz")

ANALYSIS_STEP="all"
PLATFORM="replit" #cluster-mgh 
PAIR_FILE="fq_pairs.csv"
SOURCE_LOC="../test_folder"
CLEAR_TMP="false"
OVER_WRITE="false"
COUNT_METHOD="seperated"
BATCH_SIZE="none" # number of total batches/forks/processes
MAX_PARALLEL=4 # max number of parallel forks/processes
MEM_SIZE="2000"
THREAD_SIZE="4" # parameter for salmon and subread
ALIGN_METHOD="salmon"
QC_METHOD="fastp"
RUN_TYPE="-r"
CHILD="false" # meant for multuple bsub, not used
EXIT_ON_SINGLE="false"
RUN_MODE="local"

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-REF_ANNOTATION) REF_ANNOTATION_LOC="$2"; shift ;;
		-REF_GENOME) REF_GENOME_LOC="$2"; shift ;;
		-REF_TRANSCRIPT) REF_TRANSCRIPT_LOC="$2" ; shift ;;
		-TRANSCRIPTS) TRANSCRIPT_LOCS+=("$2") ; shift ;;
		-ANNOT_DIR) ANNOT_DIR="$2"; shift ;;
		-ANNOT_LOCS) ANNOT_LOCS="$2"; shift ;;
		-COMB_ANNOTATION) COMB_ANNOTATION="$2"; shift ;;
		-SOURCE) SOURCE_LOC="$2"; shift ;;
		-ANALYSIS_STEP) ANALYSIS_STEP="$2"; shift ;; #steps: index, convert, qc, align, count
		-PLATFORM) PLATFORM="$2"; shift ;;
		-COUNT_METHOD) COUNT_METHOD="$2"; shift ;;
		-BATCH_SIZE) BATCH_SIZE="$2"; shift ;;
		-MAX_PARALLEL) MAX_PARALLEL=$2 ; shift ;;
		-CLEAR_TMP) CLEAR_TMP="true" ;;
		-OVER_WRITE) OVER_WRITE="true" ;;
		-CHILD) CHILD="true" ;;
		-d) RUN_TYPE="-d" ;;
		-r) RUN_TYPE="-r" ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done

# TRANSCRIPT_LOCS=(${TRANSCRIPT_LOCS//;;/ })
# ANNOT_LOCS=(${ANNOT_LOCS//;;/ })

# printf -v "$CONFIG_OPTION" "%s" "$CONFIG_VALUE"

main_loc="$PWD"
chmod -R +x ./


if [ ! -d "../hERV_Work" ]; then 
	mkdir ../hERV_Work 
fi 

cd ../hERV_Work

# read -p "starting, enter programme run type: " RUN_TYPE
if [ "$RUN_TYPE" = "-d" ]; then # have to keep spaces between square brackets and vars
	read -p "enter debug type: " DEBUG_TYPE
	# DEBUG_TYPE='-f' #function
	if [ "$DEBUG_TYPE" = '-f' ]; then
		read -p "enter function to debug: " FUNC
		${FUNC[0]} ${FUNC[@]:1}
		echo $?
	fi

	if [ "$DEBUG_TYPE" = '-u' ]; then
		read -p "enter file to debug: " FILE
		bash ${main_loc}/$FILE ${FILE[@]:1}
		echo $?
	fi
elif [ "$RUN_TYPE" = "-r" ]; then
	if [ $CHILD = "false" ]; then
#		SEQ_TYPE=("RNA-Seq" "WXS") #used for TCGA
		. $main_loc/setup.sh 
		#	conda env list
		. $main_loc/downloads.sh 
		if [ $SOURCE_LOC = "tcga" ] && [ ! -d "tcga" ]; then 
			. $main_loc/tcga.sh
		fi
		if [ ! $BATCH_SIZE = "none" ]; then 
			. $main_loc/batch.sh
		else 
			. $main_loc/analysis.sh
		fi
	elif [ $CHILD = "true" ]; then # remenant from job parallisation attempt
			echo 'this should not run!!'
      . $main_loc/analysis.sh
	fi
	clear_tmp
fi

# jq --arg rty "$RUN_TYPE" \
# 	 --arg dty "$DEBUG_TYPE" \
# 	 --arg plf "$PLATFORM" \
# 	 --arg 
# 	 '{"rty": $rty, "dty": $dty, "plf": $plf}'
#{"filters":{"op":"=","content":{"field":"cases.demographic.gender","value":["female"]}},"size":1}
