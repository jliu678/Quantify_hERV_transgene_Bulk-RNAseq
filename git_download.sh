#!/bin/bash

branch_name="main"
commit_message="fix"

while [ -n "$1" ]; do #setting variables for sub-processes
	case "$1" in
		-branch_name) branch_name="$2"; shift ;;
		-commit_message) commit_message="$2"; shift ;;
#		-REF_TRANSCRIPT) REF_TRANSCRIPT_LOC="$2" ; shift ;;
#		-hERV_TRANSCRIPT) hERV_TRANSCRIPT_LOC="$2" ; shift ;;
#		-hERV_DIR) hERV_DIR="$2"; shift ;;
#		-hERV_FILE) hERV_FILE="$2"; shift ;;
#		-COMB_ANNOTATION) COMB_ANNOTATION="$2"; shift ;;
#		-SOURCE) SOURCE_LOC="$2"; shift ;;
#		-ANALYSIS_STEP) ANALYSIS_STEP="$2"; shift ;; #steps: index, convert, qc, align, count
#		-PLATFORM) PLATFORM="$2"; shift ;;
#		-COUNT_METHOD) COUNT_METHOD="$2"; shift ;;
#		-BATCH_SIZE) BATCH_SIZE="$2"; shift ;;
#		-MAX_PARALLEL) MAX_PARALLEL=$2 ; shift ;;

#		-CLEAR_TMP) CLEAR_TMP="true" ;;
#		-OVER_WRITE) OVER_WRITE="true" ;;
#		-CHILD) CHILD="true" ;;
#		-d) RUN_TYPE="-d" ;;
#		-r) RUN_TYPE="-r" ;;
		*) echo "$1 is not an option"; exit 1 ;;
	esac
	shift
done




git reset --hard
git pull origin $branch_name