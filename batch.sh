#!/bin/bash

batch_op() {
	mkdir batches
	lines=$(wc -l < "$PAIR_FILE")
	begin=0; len=$BATCH_SIZE; no=0
	until [[ $((lines-begin)) <= 0 ]]; do
		tail -n +$begin "$PAIR_FILE" | head -$len > "batches/batch_file$no"
		(($no++)); ((begin+=$len))
	done 
}

main() {
	batch_op
	for i in "batches/*"; do
		local batch_name=$i
		. ./analysis.sh $batch_name
	done 
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main "$@"
fi