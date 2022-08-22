#!/bin/bash

batch_op() {
	if [ ! -d batches ]; then
		mkdir batches
		lines=$(wc -l < "$PAIR_FILE"); 
		# ((lines=$lines/$BATCH_SIZE))
		split -l $BATCH_SIZE "$PAIR_FILE" "batches/batch_file_"
	fi 

	# begin=0; len=$BATCH_SIZE; no=0
	# until [[ $((lines-begin)) -le 0 ]]; do
	# 	tail -n +$begin "$PAIR_FILE" | head -$len > "batches/batch_file$no"
	# 	((no++)); ((begin+=$len))
	# done 
}

if [ ! -d "logs" ]; then
	mkdir logs
fi 

main() {
	local tmp_ast=$ANALYSIS_STEP
	ANALYSIS_STEP="index,convert"
	. $main_loc/analysis.sh 
	ANALYSIS_STEP="$tmp_ast"
	
	batch_op
	my_jobs=()
	local tmp_PAIR_FILE="$PAIR_FILE"
	id=$(ps ax -O pgid | grep main.sh)
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

		PAIR_FILE=$i; CHILD="true"
		# bsub < (main_loc=$main_loc envsubst <test_batch.lsf)
		# (trap "kill 0" SIGINT SIGTERM ; . $main_loc/analysis.sh) &
		(. $main_loc/analysis.sh) &
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