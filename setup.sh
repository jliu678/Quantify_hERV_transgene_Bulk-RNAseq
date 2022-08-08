#!/bin/bash
#download and setup all the programmes needed for analysis

source timed.sh

setup_conda() {
	if ! command -v conda ; then
		if [ "$PLATFORM" = "cluster-mgh" ]; then
			module load anaconda
			conda init bash
		else
			if [[ ! $(python -V 2>&1) == *"Python 3.9"* ]]; then
				apt-get install python3.9
			fi
			timed_download https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
			bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
			exec bash
		fi
	fi 
}

setup_tools(){
	if [ ! -d "fastp" ] ; then 
		timed_download http://opengene.org/fastp/fastp
		chmod a+x ./fastp
	fi 
	
	#tools that could be found on anaconda
	conda_tools=(
		"samtools=1.15.1"
		"subread=2.0.1"
	)
	
	for i in ${conda_tools[@]}; do
		if [ "$PLATFORM" = "cluster-mgh" ]; then
			local mgh_name = $("$i" | tr = /)
			if [[ $(module avail mgh_name) = *"mgh_name"* ]]; then
				module load $i
			fi 
		fi

		if [ ! command -v $i ] ; then 
			conda install -c bioconda $i
		fi
	done 
}

main(){
	setup_conda
	setup_tools
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main "$@"
fi