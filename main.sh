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
SOURCE_LOC=""

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-REF_ANNOTATION) REF_ANNOTATION_LOC="$2"; shift ;;
		-REF_GENOME) REF_GENOME_LOC="$2"; shift ;;
		-hERVd_DIR) hERVd_DIR="$2" shift ;;
		-hERV_FILE) hERV_FILE="$2" shift ;;
		-COMB_ANNOTATION) COMB_ANNOTATION="$2" shift ;;
		-SOURCE) SOURCE="$2" shift ;;
		-ANALYSIS_STEP) ANALYSIS_STEP="$2" ;; #steps: index, convert, qc, align, count
		-PLATFORM) PLATFORM="$2" ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done

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
elif [ "$RUN_TYPE" = "-r" ]
	SEQ_TYPE=("RNA-Seq" "WXS")
	( . ./setup.sh )
	( . ./downloads.sh )
	if [ $SOURCE = "tcga" ] && [ ! -d "tcga" ]; then 
		(	. ./tcga.sh )
	fi
fi

# jq --arg rty "$RUN_TYPE" \
# 	 --arg dty "$DEBUG_TYPE" \
# 	 --arg plf "$PLATFORM" \
# 	 --arg 
# 	 '{"rty": $rty, "dty": $dty, "plf": $plf}'
#{"filters":{"op":"=","content":{"field":"cases.demographic.gender","value":["female"]}},"size":1}
