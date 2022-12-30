#!/bin/bash
#main programme run file, bam -(samtools)> fastq -(fastp)> fastq -(ANY)> bam -(ANY)> counts

. $main_loc/timed.sh

# change the source location if it is a relative path
if [ ! ${SOURCE_LOC[0]} = "/" ]; then
	SOURCE_LOC="../$SOURCE_LOC"
fi

# get the file names only and w/o the extensions
SOURCE=$(basename $SOURCE_LOC)
REF_GENOME=$(basename $REF_GENOME_LOC .gz) # for salmon decoy
REF_ANNOTATION=$(basename $REF_ANNOTATION_LOC .gz) # for subread
REF_TRANSCRIPT=$(basename $REF_TRANSCRIPT_LOC .gz) # for salmon quant, fasta file
TRANSCRIPTS=$(basename $(basename $TRANSCRIPT_LOCS .gz) .bz2) # for salmon quant, fasta file

split_fastq(){ #splits pair ended fastq from bam into 2 files
	# check the bam file if it is pair ended
	# get the flags for the file
	local bam_flags=$(samtools view -c -f 1 "${SOURCE_LOC}/$1.bam")
	# convert the flags into words
	local str_flags=$(samtools flags $bam_flags)

	# if they are paired, split the file into the constituent pairs
	# generate the pairs from the flags for each transcript (either ending in 1 or 2)
	# we also zip up the files to save space and for consistency
	if  [[ $str_flags == *"PAIRED"* ]]; then 
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/1$' -A 3 --no-group-separator | \
				gzip > "tmp/${SOURCE}/$1.r1.fq.gz"
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/2$' -A 3 --no-group-separator | \
				gzip > "tmp/${SOURCE}/$1.r2.fq.gz"
		echo -e "$1.r1,$1.r2\n" >> $PAIR_FILE

	#if we should keep single ended files, then we compress the fq and
	elif [ $EXIT_ON_SINGLE = "false" ]; then 
		gzip -c tmp/${SOURCE}/$1.fq > tmp/${SOURCE}/$1.fq.gz
		echo -e "$1,\n" >> $PAIR_FILE 
	fi
	rm "tmp/${SOURCE}/$1.fq"
}

bam_to_fastq(){ #wrapper, function as the name suggests
	if [ "$OVER_WRITE" = "true" ] || [ ! -e "tmp/${SOURCE}/$1.r1.fq.gz" ] || \
		 [ ! -e "tmp/${SOURCE}/$1.r2.fq.gz" ] || \
		 [ ! -e "tmp/${SOURCE}/$1.fq.gz" ]
	then
		samtools bam2fq "${SOURCE_LOC}/$1.bam" > "tmp/${SOURCE}/$1.fq"
		split_fastq $1
	fi
}

get_r2name(){
	local r2id="r2"
	if echo $1 | grep -q "R1"; then r2id="R2"; fi
	echo $1 | sed "s/r1/$r2id/i"
}

check_name(){
	echo $1 | grep -iq "r1"
	local is_r1=$?
	return $is_r1 && [ -e "$(get_r2name $1)" ]
}

group_fastq(){ #group fastq files into pairs
	timed_print "grouping fastq files in tmp/${SOURCE}/"
	files=(tmp/${SOURCE}/*.fq.gz)
	local total=0
	for i in "${!files[@]}"; do #loop w/ index bcs easier
		# local file_name=$(echo "${files[$i]}" | cut -f 1 -d '.') 
		local file_name=${files[$i]%%.*} #get file name w/o extention
		# local file_ext=${files[$i]##*.} # not sure why this is here.. :(

		if ! grep -q "^$(basename $file_name .fq.gz)$" $PAIR_FILE ; then #if the file does not have pair
			if check_name ${files[$i]}; then #if formatted correctly
        timed_print "===!!!looks like formatted corredtly!!!===" 
				echo -e "$(basename $file_name),$(get_r2name $(basename $file_name))" >> $PAIR_FILE
			elif ! echo "${files[$i]}" | grep -iq "r2" ; then #the choice is yours how to deal with single-ended files
				timed_print "compliment to ${files[$i]} not found"
				((total+=1))

				if [ $EXIT_ON_SINGLE = "false" ]; then # should single ended files will be included
					echo -e "$(basename $file_name)," >> $PAIR_FILE
				fi
			fi
		fi
	done

	if [[ $total -ge 1 && $EXIT_ON_SINGLE = "true" ]]; then  # if single ended files result in errors
		exit 1
	fi

	timed_print "grouped fastq files, found total $total files without compliment"
}

mv_fq() {
	if [[ "$OVER_WRITE" = "true"  ||  ! -e "$2" ]]; then
		case "${name#*.}" in 
			fq) gzip -c $1 > $2 ;;
			fq.gz) ln -rs $1 $2 ;;
		esac
	fi
}

get_pairs_all() { #place all files into tmp, group them
	if [[ ! -d "tmp/${SOURCE}" ]]; then 
		mkdir "tmp/${SOURCE}"
	fi 
	touch $PAIR_FILE

	for i in ${SOURCE_LOC}/*; do
		timed_print "moving $i..."
		local name=$(basename $i)
		local only_name=${name%.*} #removes last extention, ie bam or gz
		case "${name#*.}" in #get extention
			bam) bam_to_fastq $only_name ;; 
			fq) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/${name}.gz" ;; #move and zip for smaller size 
			fq.gz) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/${name}" ;; #move and does not change data
		esac
		timed_print "moved $i"
	done 
	group_fastq #groups all other files
}

fastp_qc(){ 
	# results stored in /tmp/{name of the source}/qc/* and all end in *.qc.fq.gz
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/qc/$1.qc.fq.gz" ]; then
		if [[ $# -eq 2 ]]; then # if pair ended
			./fastp -i "tmp/${SOURCE}/$1.fq.gz" -o "tmp/${SOURCE}/qc/$1.qc.fq.gz" \
				-I "tmp/${SOURCE}/$2.fq.gz" -O "tmp/${SOURCE}/qc/$2.qc.fq.gz" \
				-j "results/fastp/$1.json" -h "results/fastp/$1.html"
		elif [[ $# -eq 1 ]]; then # if single ended
			./fastp -i "tmp/${SOURCE}/$1.fq.gz" -o "tmp/${SOURCE}/qc/$1.qc.fq.gz" \
				-j "results/fastp/$1.json" -h "results/fastp/$1.html"
		fi	
	fi
}

qc_all(){ # only supports fastp as of now
	if [[ ! -d "tmp/${SOURCE}/qc" ]]; then 
		mkdir "tmp/${SOURCE}/qc"
	fi

	if [[ ! -d "results/${QC_METHOD}" ]]; then 
		mkdir "results/${QC_METHOD}"
	fi

	timed_print "qc-ing with $QC_METHOD"
	while IFS=, read -r r1 r2; do # loop over pair file, csv
		if [ "$OVER_WRITE" = "true" ] || [ ! -d "results/${ALIGN_METHOD}/$r1" ]; then 
			#if salmon quant already exists

			timed_print "qc-ing ${r1} and ${r2}"
			case $QC_METHOD in 
				fastp) fastp_qc $r1 $r2 ;;
			esac 
		fi
	done < $PAIR_FILE
}

subread_build_index(){
	if [ ! -d subread/${REF_GENOME%.*}_index  ]; then
		timed_print "building subread index @: subread/${REF_GENOME%.*}_index"
		subread-buildindex -o subread/${REF_GENOME%.*}_index $REF_GENOME
	else 
		timed_print "index already exists @: subread/${REF_GENOME%.*}_index"
	fi
}

salmon_build_index(){ # supports multiple transcripts
	# get dir name for the index, which is the transcript names seperated by a hyphen(-)
	# eg. for erv.fa and line.fa, you get: erv-line
	local index_name="salmon/$(IFS=-; echo "${TRANSCRIPTS[*]%.*}")" 
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "${index_name}_index" ]; then 
		timed_print "building salmon index @: ${index_name}_index"

		# salmon requires a decoy, which in this case is the entire reference genome
		if [ "$OVER_WRITE" = "true" ] || [ ! -f "${index_name}_decoys.txt" ]; then
			timed_print "building decoys..."
			grep "^>" "${REF_GENOME}" | cut -d " " -f 1 > "${index_name}_decoys.txt"
			sed -i.bak -e 's/>//g' "${index_name}_decoys.txt"
		fi 

		# salmon requires all the transcripts being quantified to be concat into one file
		if [ "$OVER_WRITE" = "true" ] || [ ! -f "${index_name}_gentrome.fa" ]; then
			timed_print "building gentrome..."
			cat "${TRANSCRIPTS[@]}" "${REF_TRANSCRIPT}" "${REF_GENOME}" > "${index_name}.fa"
		fi

		# adding --gencode flag due to salmon doing some preprocessing for it
		# this flag does nothing for none gencode files (I think)
		salmon index -t ${index_name}.fa -d "${index_name}_decoys.txt" -i "${index_name}_index" --gencode
	fi
}

build_index() {
	case $ALIGN_METHOD in
		subread) subread_build_index ;;
		salmon) salmon_build_index ;;
	esac
}

subread_align(){
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/subread_aligned/$1.bam" ]; then
		timed_print "aligning $1.qc.fq and $2.qc.fq..."
		subread-align -i subread/${REF_GENOME%.*}_index -r "tmp/${SOURCE}/qc/$1.qc.fq" -R "tmp/${SOURCE}/qc/$2.qc.fq" -t 0 -o "tmp/${SOURCE}/subread_aligned/$1.bam" -T $THREAD_SIZE -M MEM_SIZE 
		timed_print "aligned tmp/${SOURCE}/qc/$1.qc.fq and tmp/${SOURCE}/qc/$2.qc.fq"
	fi
}

align_all(){
	case $ALIGN_METHOD in
		subread)
			if [ ! -d "tmp/${SOURCE}/${ALIGN_METHOD}_aligned" ]; then 
				mkdir "tmp/${SOURCE}/${ALIGN_METHOD}_aligned"
			fi ;;
		salmon) ;;
	esac
	
	while IFS=, read -r r1 r2; do
		case $ALIGN_METHOD in 
			subread) subread_align $r1 $r2 ;;
			salmon) return 0 ;; # salmon does not need alignment, so it does nothing at this step
		esac
	done < $PAIR_FILE
}

get_gene_types_comb(){
	gene_types=$( cut -f3 ${COMB_ANNOTATION}.gff3 | grep -v "^#" | sort | uniq )
}

subread_count_comb(){
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "results/subread/fc_comb_$1.tsv" ]; then 
		timed_print "counting $1.bam..."
		featureCounts -a ${COMB_ANNOTATION}.gff3 -o "results/subread/fc_comb_$1.tsv" "tmp/${SOURCE}/subread_aligned/$1.bam" -T 6 -t ${gene_types[@]} -g "ID"
		timed_print "counted $1.bam"
	fi
}

get_gene_types_sep(){
	ref_gene_types=( $(cut -f3 ${REF_ANNOTATION} | grep -v "^#" | sort | uniq) )
	erv_gene_types=( $(cut -f3 ${hERV_FILE} | grep -v "^#" | sort | uniq) )
}

subread_count_sep(){
	local IFS=,
	echo "${ref_gene_types[*]}"
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "results/subread/fc_ref_$1.tsv" ]; then 
		timed_print "counting $1.bam with ${REF_ANNOTATION}..."
		featureCounts -a ${REF_ANNOTATION} -o results/subread/fc_ref_$1.tsv "tmp/${SOURCE}/subread_aligned/$1.bam" -T $THREAD_SIZE -t "${ref_gene_types[*]}" -g ID
		timed_print "counted $1.bam"
	fi 
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "results/subread/fc_erv_$1.tsv" ]; then 
		timed_print "counting $1.bam with ${hERV_FILE}..."
		featureCounts -a ${hERV_FILE} -o results/subread/fc_erv_$1.tsv "tmp/${SOURCE}/subread_aligned/$1.bam" -T $THREAD_SIZE -t "${erv_gene_types[*]}" -g ID
		timed_print "counted $1.bam"
	fi
}

subread_count() {
	case $COUNT_METHOD in
		combined) 
			get_gene_types_comb
			for r1 in $(cut -d, -f1 < ${PAIR_FILE}); do 
				subread_count $r1
			done ;;
		seperated)
			get_gene_types_sep
			for r1 in $(cut -d, -f1 < ${PAIR_FILE}); do 
				subread_count_sep $r1
			done ;;
	esac
}

salmon_quant() {
	# get the same index_name as specified in "salmon_build_index()"
	local index_name="salmon/$(IFS=-; echo "${TRANSCRIPTS[*]%.*}")"
	#salmon outputs into a dir, so we check for that instead
	#results stored in results/salmon/*
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "results/salmon/$1" ]; then
		mkdir results/salmon/$1
		if [[ $# -eq 2 ]]; then # if pair ended
			salmon quant -i "${index_name}_index" -l A -1 "tmp/${SOURCE}/qc/$1.qc.fq.gz" -2 "tmp/${SOURCE}/qc/$2.qc.fq.gz" --validateMappings -o "results/salmon/$1"
			#sleep 1
		elif [[ $# -eq 1 ]]; then # if single ended
			salmon quant -i "${index_name}_index" -l A -r "tmp/${SOURCE}/qc/$1.qc.fq.gz" --validateMappings -o "results/salmon/$1"
			#sleep 1
		fi 
	fi
}

salmon_count() {
	while IFS=, read -r r1 r2; do
		salmon_quant $r1 $r2
	done < $PAIR_FILE
}

count_all(){
	if [ ! -d "results/${ALIGN_METHOD}" ]; then 
		mkdir "results/${ALIGN_METHOD}"
	fi
	
	case $ALIGN_METHOD in 
		subread) subread_count ;;
		salmon) salmon_count
	esac
}

main(){
	timed_print "starting $$"
	timed_print "analysing..."

	# make these essential dirs
	local directories=( "tmp" "$ALIGN_METHOD" "results" )
	for i in ${directories[@]}; do
		if [[ ! -d "$i" ]]; then 
			mkdir $i
		fi
	done

	if [ "$ANALYSIS_STEP" = "all" ]; then 
		ANALYSIS_STEP="index,convert,qc,align,count"
	fi 

	#turn "index,c,q" into ( "index" "c" "q" )
	ANALYSIS_STEP=(${ANALYSIS_STEP//,/ }) 

	for i in ${ANALYSIS_STEP[@]}; do
	  timed_print "$i-ing..."	
		case "$i" in 
			# for the first 2, do not run if they are in a child/forked process
			index) if [ $CHILD = false ]; then build_index; fi ;; 
			convert) if [ $CHILD = false ]; then get_pairs_all; fi ;;
			qc) qc_all ;;
			align) align_all ;;
			count) count_all ;;
		esac 
		timed_print "finished $i"
	done 

	if [[ $CHILD = true &&  $CLEAR_TMP = 'true' ]]; then
		while IFS=, read -r r1 r2; do
			rm "tmp/${SOURCE}/$r1.fq.gz"
			rm "tmp/${SOURCE}/qc/$r1.qc.fq.gz"

			if [ -n $r2 ]; then 
				rm "tmp/${SOURCE}/$r2.fq.gz"
				rm "tmp/${SOURCE}/qc/$r2.qc.fq.gz"
			fi
		done < $PAIR_FILE 
	fi

	# get the uage metrics before exiting
	timed_print "exiting $$ with: "
	top -b -n 1 -u $USER
	/usr/local/bin/pan_quota $HOME
}

main


#subread-align -i siyi_link/subread/hg38_index -r siyi_link/HcxecP3_1_r1.fastq -R siyi_link/HcxecP3_1_r2.fastq -t 0 -o siyi_link/subread/HcxecP3_1_aligned.bam -T 8 -M 16000

#featureCounts -a gencode.v41.chr_patch_hapl_scaff.annotation.gff3 -o siyi_link/subread/results/feature_count_uncombined.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6

#featureCounts -a erv.gff3 -o siyi_link/subread/results/feature_count_erv.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6 -t BED_feature,exon -g ID

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt 
#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt
#cat package-entities-erv.fa gencode.v41.transcripts.fa GRCh38.p13.genome.fa > gentrome.fa
#salmon index -t gentrome.fa -d decoys.txt -p 6 -i index --gencode
#salmon quant -i index -l A -1 HcxecP3_1_r1.fastq -2 HcxecP3_1_r2.fastq --validateMappings -o results

#gffcompare -R -r gencode.v41.chr_patch_hapl_scaff.annotation.gtf -o siyi_link/combine_annotation/herv_annotation erv.gff3 >>warnings.txt 2>&1
#cuffmerge -g gencode.v41.chr_patch_hapl_scaff.annotation.gtf -s GRCh38.p13.genome.fa -o siyi_link/combine_annotation/herv_annotation gff_list.txt

# usage=$(awk '{u=$2+$4; t=$2+$4+$5; if (NR==1){u1=u; t1=t;} else print ($2+$4-u1) * 100 / (t-t1) "%"; }' <(grep 'cpu ' /proc/stat) <(sleep 0.5;grep 'cpu ' /proc/stat))
