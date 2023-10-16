#!/bin/bash
#download and setup all the programmes needed for analysis

. $main_loc/timed.sh

setup_conda() {
	# get conda if it doesn't exist
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
		fi
		# restart bash bcs conda needs it
		exec bash
	fi 

	# create env if it doesn't exist
	if ! conda info --envs | grep "$PLATFORM" >/dev/null 2>&1; then
		conda create -n $PLATFORM -y
	fi
	# !!!magic!!! (`',^,'`)
	eval "$(conda shell.bash hook)"
	conda activate $PLATFORM
}

setup_tools(){
	if [ ! -f "fastp" ] ; then # fastp needs special downloading
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
		 cmd_name=$(echo $i | cut -d'=' -f1) # get the name w/o the version control stuff
     echo "dealing with $cmd_name"
    
		 if [ "$PLATFORM" = "cluster-mgh" ]; then
		 	local mgh_name=$(echo "$i" | tr = /)
      echo $mgh_name
      if module load "$mgh_name" >/dev/null 2>&1; then 
        echo "module is avail: $mgh_name"
        module list  
		 	elif ! conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1; then 
          echo "cluster-mgh conda has no $i" 
          conda install -c conda-forge -c bioconda -y $i
          echo "cluster-mgh conda installed/passed $cmd_name"   
       fi 
		 elif ! conda list "$cmd_name" | grep "$cmd_name" >/dev/null 2>&1 ; then
			    echo "conda has no $cmd_name in $PLATFORM"  
			    conda install -c conda-forge -c bioconda -y $i
          echo "conda installed/passed $cmd_name in $PLATFORM" 
		fi
	done
}

main(){
	setup_conda
	setup_tools
}

main
