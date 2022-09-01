#!/bin/bash
test_module_load() {
echo $1
local cmd_name=$(echo $1 | cut -d'/' -f1)
echo $cmd_name
if ! conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1; then 
              echo "conda has no $cmd_name in cluster-mgh" 
              elif conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1; then
              echo "conda has $cmd_name h" 
              fi
              
}

test_module_load $1


check_core_dump() {
echo ========= $1 $2 ==============
	if [[ $1 -eq 0 ]]; then
		#for i in core.*; do 
		#	rm $i
		#fi
echo psudoForLoopRm

		need_to_restart=()
check_quant_sf=$2

		if [[ $check_quant_sf -eq 0 ]]; then
			for i in {1..7}; do 
				#wait_batches
				#pair_file=$(grep "$i" batches/*)
				#run_batch "$pair_file"
				echo psudoProcessingDumpedFq
			done
			#wait ${my_jobs[@]}
			return 0
		fi
	fi 
 echo "i think 1 returned"
	return 1
 
}

#if there is core dump, $1 = 0; if theres folder wihtout .sf, $2 =0 

#check_core_dump 0 0 && echo 0
#echo $?
#
#check_core_dump 0 1 && echo 0
#echo $?
#
#check_core_dump 1 0 && echo 0
#echo $?
#
#check_core_dump 1 1 && echo 0
echo $?