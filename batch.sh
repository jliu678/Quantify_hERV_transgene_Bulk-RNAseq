#!/bin/bash
. $main_loc/timed.sh

batch_op() {
	if [ ! -d batches ]; then
		mkdir batches
		lines=$(wc -l < "$PAIR_FILE"); 
		# ((lines=$lines/$BATCH_SIZE))
		split -l $BATCH_SIZE "$PAIR_FILE" "batches/batch_file_"
	fi 
}

setup_lsf() { # !!!UNTESTED!!! support for submitting each "batch" as its own job
	if [ "$RUN_MODE" = "lsf" ]; then 
		for i in batches/*; do 
			pairs="$PAIR_FILE" steps="$ANALYSIS_STEPS" envsubst <test_batch.lsf | echo > "batch_lsf/$(basename $i)"
		done
	fi
}

run_batch() { 
	case $RUN_MODE in
		lsf) bsub < "batch_lsf/$(basename $1)" ;;
		local) # changes the PAIR_FILE for the forked/child process, add trap to kill when parent is
			PAIR_FILE=$i; CHILD="true"
			(trap "kill 0" SIGINT ; . $main_loc/analysis.sh) & 
			my_jobs+=( $! ) # array, used to keep track of how many jobs being run
		;;
	esac
}

wait_batches() { # wait for existing batches (limited by MAX_PARALLEL) to exit before spawning new ones
	while [[ ${#my_jobs[@]} -ge $MAX_PARALLEL ]]; do
		timed_print ${my_jobs[@]}
		sleep 60
		for j in ${!my_jobs[@]}; do 
			if ! ps -p "${my_jobs[$j]}" ; then # if pid does not exist
				echo ${my_jobs[$j]} is not running
				unset my_jobs[$j] # change the PID to empty string
			fi
		done
		my_jobs=("${my_jobs[@]}") # convert array with empty elements to new consequtive
	done
}

launch_batches() {
	my_jobs=()
	tmp_PAIR_FILE="$PAIR_FILE"

	# adds trap signal for the entire group, haven't used it ('~')
	# group_id=$(ps ax -O pgid | grep main.sh)
	# trap "kill -SIGINT -- -$group_id" SIGINT

	for i in batches/*; do
		wait_batches
		run_batch $i
		# (trap "kill 0" SIGINT ; . $main_loc/analysis.sh) &
		# (. $main_loc/analysis.sh) &
	done
	wait ${my_jobs[@]}
}

check_quant_sf() {
	rtn=1
	for i in /results/salmon/*; do 
		if [ ! -f $i/quant.sf ]; then 
			rm -r $i
			need_to_restart+=($i) # init and used in "check_core_dump()"
			rtn=0
		fi 
	done
	return $rtn 
}

check_core_dump() { # remove core dumps and start over
	if compgen -G "core.*" > /dev/null; then # if files begining with "core." (core dumps) exists
		for i in core.*; do 
			rm $i
		done

		need_to_restart=() # array of .fq files that doesn't have quant.sf files
		if check_quant_sf; then
			for i in $need_to_restart; do # run them again... )'-'(
				wait_batches
				pair_file=$(grep "$i" batches/* | cut -d":" -f1) # get the proper pair file from the 
				run_batch "$pair_file" 
			done
			wait ${my_jobs[@]}
			return 0
		fi
	fi 
	return 1 
}

loop_until_finished() {
	my_jobs=()
	while check_core_dump; do # run the check_core_dump function until there are no more core dumps
		timed_print "there were core dumps :("
	done
}

main() {
	local tmp_ast=$ANALYSIS_STEP
	ANALYSIS_STEP="index,convert"
	. $main_loc/analysis.sh #source and run analysis.sh give ANALYSIS_STEP="index,convert"
	ANALYSIS_STEP="$tmp_ast"
	
	batch_op
	setup_lsf

	launch_batches
	loop_until_finished

	CHILD="false"
	PAIR_FILE="$tmp_PAIR_FILE"
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
# bash main.sh -PLATFORM "cluster-mgh" -ANALYSIS_STEP "all" -SOURCE "raw_data"