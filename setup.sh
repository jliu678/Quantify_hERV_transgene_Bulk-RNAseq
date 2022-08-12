#!/bin/bash
#download and setup all the programmes needed for analysis

. $main_loc/timed.sh

setup_conda() {
	if ! command -v conda ; then
	echo $? ! command -v conda
		if [ "$PLATFORM" = "cluster-mgh" ]; then
			module load anaconda
echo $? module load anaconda
			conda init bash
echo $? conda init bash
		else
			if [[ ! $(python -V 2>&1) == *"Python 3.9"* ]]; then
				apt-get install python3.9
			fi
			timed_download https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
			bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
		fi
		exec bash
	fi 
	if ! conda info --envs | grep "$PLATFORM" >/dev/null 2>&1; then
		#conda create -n $PLATFORM -y && conda activate $PLATFORM -y
conda create -n $PLATFORM -y
echo $? conda create -n $PLATFORM -y
conda activate $PLATFORM -y
echo $? conda activate $PLATFORM -y
	fi
}

setup_tools(){
	if [ ! -f "fastp" ] ; then 
echo $? ! -f "fastp"
		timed_download http://opengene.org/fastp/fastp
		chmod a+x ./fastp
	fi 
	
	#tools that could be found on anaconda
	conda_tools=(
		"samtools=1.15.1"
		"subread=2.0.1"
		"salmon=1.9.0"
	)
	
	for i in ${conda_tools[@]}; do
		# if [ "$PLATFORM" = "cluster-mgh" ]; then
		# 	local mgh_name=$("$i" | tr = /)
		# 	if [[ $(module avail $mgh_name) = *"$mgh_name"* ]]; then
		# 		module load $i
		# 	fi 
		# fi 

		cmd_name=$(echo $i | cut -d'=' -f1)
		if ! conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1 ; then 

echo $? ! conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1

			conda install -c conda-forge -c bioconda -y $i

echo wtf

		fi
	done
}

main(){
	setup_conda
	setup_tools
}

main
