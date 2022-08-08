#!/bin/bash

timed_print() {
	echo "$(date +"%D %T"): ${@:3}"
}

timed_download() {
	local name=$(basename $1)

	if [ -f "$name" ]; then
		timed_print "${name} already exists"
	else
		timed_print "downloading ${name}..."
		wget $1
	
		if [[ "$name" == *.gz ]]; then
			timed_unzip "$name"
			name=${name::-3}
		fi
	fi

	if [ ! -n "$2" ]; then 
		"$2"=$name
	fi
}

timed_unzip() {
	timed_print "unzipping $1..."
	gunzip $1
}