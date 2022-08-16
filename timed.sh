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
		if [ $PLATFORM = "cluster-mgh" ]; then
			curl -L --user-agent "Mozilla/5.0" --url $1 --remote-name
		else
			wget --user-agent="Mozilla" $1
		fi
		if [[ "$name" == *.gz ]]; then
			timed_unzip "$name"
			name="$unzip_name"
		fi
	fi

	if [[ $# = 2 ]]; then 
		"$2"=$name
	fi
}

timed_unzip() {
	timed_print "unzipping $1..."
	gunzip $1
}

#curl --user-agent "Mozilla" --url http://opengene.org/fastp/fastp -O 