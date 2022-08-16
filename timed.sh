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
		case "$name" in
			*.gz) timed_gunzip $name ;;
			*.bz2) timed_bunzip2 $name ;;
		esac
	fi

	if [[ $# = 2 ]]; then 
		"$2"=$name
	fi
}

timed_gunzip() {
	timed_print "unzipping $1..."
	gunzip $1
}

timed_bunzip2(){
	timed_print "unzipping $1..."
	bunzip2 $1
}

#curl --user-agent "Mozilla" --url http://opengene.org/fastp/fastp -O 