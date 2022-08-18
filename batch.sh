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
	ANALYSIS_STEP="convert"
	. $main_loc/analysis.sh 
	ANALYSIS_STEP="$tmp_ast"

	batch_op
	my_jobs=()
	local tmp_PAIR_FILE="$PAIR_FILE"
	for i in batches/*; do
		if [[ ${#my_jobs[@]} -ge $MAX_PARALLEL ]]; then
			timed_print ${my_jobs[@]}
			# sleep 128
			for i in ${!my_jobs[@]}; do 
				if [ ! -d "/proc/${my_jobs[$i]}" ]; then 
					unset my_jobs[$i]
				fi
			done
			my_jobs=("${my_jobs[@]}")
		else
			PAIR_FILE=$i; CHILD="true"
			# bsub < (main_loc=$main_loc envsubst <test_batch.lsf)
			(trap 'kill 0' SIGINT; . $main_loc/analysis.sh) &
			my_jobs+=( $! )
		fi
	done
	CHILD="false"
	PAIR_FILE="$tmp_PAIR_FILE"
	wait ${my_jobs[@]}
	rm -r batches
}

# if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main # "$@"
# fi