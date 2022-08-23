#!/bin/bash

batch_op() {
	if [ ! -d batches ]; then
		mkdir batches
		lines=$(wc -l < "$PAIR_FILE"); 
		# ((lines=$lines/$BATCH_SIZE))
		split -l $BATCH_SIZE "$PAIR_FILE" "batches/batch_file_"
	fi 
}

setup_lsf() {
	if [ $RUN_MODE = "lsf" ]; then 
		for i in batches/*; do 
			pairs="$PAIR_FILE" steps="$ANALYSIS_STEPS" envsubst <test_batch.lsf | echo > "batch_lsf/$(basename $i)"
		done
	fi
}

run_batch() {
	case $RUN_MODE in
		lsf) bsub < "batch_lsf/$(basename $1)" ;;
		local) 
			PAIR_FILE=$i; CHILD="true"
			(trap "kill 0" SIGINT ; . $main_loc/analysis.sh) & 
		;;
	esac
}

# if [ ! -d "logs" ]; then
# 	mkdir logs
# fi 

main() {
	local tmp_ast=$ANALYSIS_STEP
	ANALYSIS_STEP="index,convert"
	. $main_loc/analysis.sh 
	ANALYSIS_STEP="$tmp_ast"
	
	batch_op
	setup_lsf

	my_jobs=()
	local tmp_PAIR_FILE="$PAIR_FILE"

	# group_id=$(ps ax -O pgid | grep main.sh)
	# trap "kill -SIGINT -- -$group_id" SIGINT

	for i in batches/*; do
		while [[ ${#my_jobs[@]} -ge $MAX_PARALLEL ]]; do
			timed_print ${my_jobs[@]}
			sleep 10
			for j in ${!my_jobs[@]}; do 
				output=$(ps -p "${my_jobs[$j]}")
				if [ ! $? -eq 0 ] ; then
					echo ${my_jobs[$j]} is not running
					unset my_jobs[$j] # remove it
				fi
			done
			my_jobs=("${my_jobs[@]}")
		done

		run_batch $i

		# (trap "kill 0" SIGINT ; . $main_loc/analysis.sh) &
		# (. $main_loc/analysis.sh) &
		my_jobs+=( $! )
	done
	CHILD="false"
	PAIR_FILE="$tmp_PAIR_FILE"
	wait ${my_jobs[@]}
	rm -r batches
}

# if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main # "$@"
# fi

# my_jobs=( ) # empty arry
# 	for ((i=0; i <= 10; )) ; do
# 		echo ${my_jobs[@]} still running
# 		if [[ ${#my_jobs[@]} -ge 4 ]]; then
# 			sleep 1
# 			for j in ${!my_jobs[@]}; do 
# 				output=$(ps -p "${my_jobs[$j]}")
# 				if [ ! $? -eq 0 ] ; then
# 					echo ${my_jobs[$j]} is not running
# #          echo ${my_jobs[$j]}
# 					unset my_jobs[$j] # remove it
# 				fi
# 			done
# 			my_jobs=("${my_jobs[@]}")
# 		else
# 			# PAIR_FILE=$i; CHILD="true"
# 			# bsub < (main_loc=$main_loc envsubst <test_batch.lsf)
# 			(trap "kill 0" SIGINT ; sleep 10; echo $! exiting) &
# 			my_jobs+=( $! )
#       ((i++))
# 		fi
# 	done
# wait ${my_jobs[@]}