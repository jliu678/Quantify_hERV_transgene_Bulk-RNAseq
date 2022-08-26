#!/bin/bash

### submit jobs for both normal and bigmem queue and run this bash to automatically kill the one that was run later

nn=0

while true; do 
	for((i=0;i<=120;i+=2)); do
    		printf "%-*s" $((i+1)) '[' | tr ' ' '#'
    		printf "%*s%3dsec\r"  $((120-i))  "]" "$i"
    		sleep 2
	done; echo
	
	((nn+=2))
	
	echo "Current time : $(date +"%T") and $nn mins have past"
	if bjobs | grep -q "RUN"; then
		echo "someone starts running"
		if bjobs | grep "bigmem" | grep -q "RUN"; then 
			bkill "$(bjobs | grep "normal" | cut -d' ' -f1)"
			echo "bigmem is running, normal is killed"
		else 
			bkill "$(bjobs | grep "bigmem" | cut -d' ' -f1)"
			echo "normal is running, bigmem is killed"
		fi
		exit 0
	fi 
done


