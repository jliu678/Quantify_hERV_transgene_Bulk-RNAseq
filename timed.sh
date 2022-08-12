#!/bin/bash

timed_print() {
	echo "$(date +"%D %T"): $@"
}

timed_download() { #download w/ printed time 
	local name=$(basename $1)
	local unzip_name=$(basename $1 .gz)

	if [ -f "$unzip_name" ]; then
		timed_print "${name} already exists"
	else
		timed_print "downloading ${name}..."
		wget --user-agent="Mozilla" $1
	
		if [[ "$name" == *.gz ]]; then
			timed_unzip "$name"
			name="$unzip_name"
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
