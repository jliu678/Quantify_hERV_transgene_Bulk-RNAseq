#!/bin/bash

batch_op() {
	if [ ! -d batches ]; then
		mkdir batches
		lines=$(wc -l < "$PAIR_FILE")
		split -l $BATCH_SIZE "$PAIR_FILE" "batches/batch_file_"
	fi 
	# begin=0; len=$BATCH_SIZE; no=0
	# until [[ $((lines-begin)) -le 0 ]]; do
	# 	tail -n +$begin "$PAIR_FILE" | head -$len > "batches/batch_file$no"
	# 	((no++)); ((begin+=$len))
	# done 
}

main() {
	. $main_loc/analysis.sh "convert"
	batch_op
	my_jobs=()
	for i in batches/*; do
		local batch_name=$i
		(trap 'kill 0' SIGINT; . $main_loc/analysis.sh "$ANALYSIS_STEP" "$batch_name") &
		my_jobs+=($!)
		timed_print ${my_jobs[@]}
	done
	wait ${my_jobs[@]}
	# rm -r batches
}

# if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main # "$@"
# fi