#!/bin/bash
#main programme run file, bam -(samtools)> fastq -(fastp)> fastq -(ANY)> bam -(ANY)> counts

. $main_loc/timed.sh

if [ ! ${SOURCE_LOC[0]} = "/" ]; then
	SOURCE_LOC="../$SOURCE_LOC"
fi
SOURCE=$(basename $SOURCE_LOC)
REF_GENOME=$(basename $REF_GENOME_LOC .gz)
REF_ANNOTATION=$(basename $REF_ANNOTATION_LOC .gz)
REF_TRANSCRIPT=$(basename $REF_TRANSCRIPT_LOC .gz)
hERV_TRANSCRIPT=$(basename $(basename $hERV_TRANSCRIPT_LOC .gz) .bz2)
id_types="ID"

split_fastq(){ #splits pair ended fastq from bam into 2 files
	local bam_flags=$(samtools view -c -f 1 "${SOURCE_LOC}/$1.bam")
	local str_flags=$(samtools flags $bam_flags)
	if  [[ $str_flags == *"PAIRED"* ]]; then 
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/1$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r1.fq"
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/2$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r2.fq"
		echo -e "$1.r1,$1.r2" >> $PAIR_FILE
		rm "tmp/${SOURCE}/$1.fq"
	fi
}

bam_to_fastq(){ #wrapper, function as the name suggests
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/$1.r1.fq" ] || [ ! -f "tmp/${SOURCE}/$1.r2.fq" ]; then
		samtools bam2fq "${SOURCE_LOC}/$1.bam" > "tmp/${SOURCE}/$1.fq"
		split_fastq $1
	fi
}

group_fastq(){ #group fastq files into pairs
	timed_print "grouping fastq files in tmp/${SOURCE}/"
	files=(tmp/${SOURCE}/*.fq)
	local total=0
	for i in "${!files[@]}"; do #loop w/ index bcs easier
		local file_name=${files[$i]%.*} #get file name w/o extention
		if ! grep -q "$file_name" $PAIR_FILE ; then #if the file does not have pair
			if [[ ${file_name: -1} = "1" && -f "${file_name::-1}2.fq" ]]; then #if formatted correctly
				echo -e "$(basename $file_name), $(basename ${file_name::-1})2" >> $PAIR_FILE
			elif [[ ! ${file_name: -1} = "2" ]]; then #the choice is yours how to deal with single-ended files
				timed_print "compliment to ${files[$i]} not found"
				total+=1
				#cat "${files[$i]}\n" >> $PAIR_FILE
			fi
		fi
	done
	timed_print "grouped fastq files, found total $total files without compliment"
}

mv_fq() {
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "$2" ]; then
		case "${name#*.}" in 
			fq) mv $1 $2 ;;
			fq.gz) gunzip -c $1 > $2 ;;
		esac
	fi
}

get_pairs_all() { #place all files into tmp, group them
	if [ ! -d "tmp/${SOURCE}" ]; then 
		mkdir "tmp/${SOURCE}"
	fi 
	touch $PAIR_FILE

	for i in ${SOURCE_LOC}/*; do
		timed_print "moving $i..."
		local name=$(basename $i)
		local only_name=${name%.*} #removes last extention, ie bam or gz
		case "${name#*.}" in #get extention
			bam) bam_to_fastq $only_name ;; 
			fq) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/$name" ;; #move because does not modify original data
			fq.gz) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/$only_name" ;; #unzip for uniformity
		esac
		timed_print "moved $i"
	done 
	group_fastq #groups all other files
}

fastp_qc(){ #only works for pair ended as of now
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/qc/$1.qc.fq" ]; then
		if [[ $# -eq 2 ]]; then 
			./fastp -i "tmp/${SOURCE}/$1.fq" -o "tmp/${SOURCE}/qc/$1.qc.fq" \
				-I "tmp/${SOURCE}/$2.fq" -O "tmp/${SOURCE}/qc/$2.qc.fq" \
				-j "results/fastp/$1.json" -h "results/fastp/$1.html"
		else 
			return 1
		fi	
	fi
}

qc_all(){
	if [ ! -d "tmp/${SOURCE}/qc" ]; then 
		mkdir "tmp/${SOURCE}/qc"
	fi

	if [ ! -d "results/${QC_METHOD}" ]; then 
		mkdir "results/${QC_METHOD}"
	fi

	timed_print "qc-ing with $QC_METHOD"
	while IFS=, read -r r1 r2; do
		timed_print "qc-ing ${r1} and ${r2}"
		case $QC_METHOD in 
			fastp) fastp_qc $r1 $r2 ;;
		esac 
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

salmon_build_index(){
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "salmon/${hERV_TRANSCRIPT%.*}_index" ]; then 
		timed_print "building salmon index @: salmon/${hERV_TRANSCRIPT%.*}_index"

		if [ "$OVER_WRITE" = "true" ] || [ ! -f "salmon/decoys.txt" ]; then
			timed_print "building decoys..."
			grep "^>" "${REF_GENOME}" | cut -d " " -f 1 > "salmon/decoys.txt"
			sed -i.bak -e 's/>//g' "salmon/decoys.txt"
		fi 

		if [ "$OVER_WRITE" = "true" ] || [ ! -f "salmon/gentrome.fa" ]; then
			timed_print "building gentrome..."
			cat ${hERV_TRANSCRIPT} ${REF_TRANSCRIPTS} ${REF_GENOME} > "salmon/gentrome.fa"
		fi

		salmon index -t salmon/gentrome.fa -d "salmon/decoys.txt" -i "salmon/${hERV_TRANSCRIPT%.*}_index" --gencode
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
	if [ ! -d "tmp/${SOURCE}/${ALIGN_METHOD}_aligned" ]; then 
		mkdir "tmp/${SOURCE}/${ALIGN_METHOD}_aligned"
	fi
	
	while IFS=, read -r r1 r2; do
		case $ALIGN_METHOD in 
			subread) subread_align $r1 $r2 ;;
			salmon) return 0 ;;
		esac
	done < $PAIR_FILE
}

get_gene_types_comb(){
	gene_types=$( cut -f3 ${COMB_ANNOTATION}.gff3 | grep -v "^#" | sort | uniq )
}

subread_count_comb(){
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "results/subread/fc_comb_$1.tsv" ]; then 
		timed_print "counting $1.bam..."
		featureCounts -a ${COMB_ANNOTATION}.gff3 -o "results/subread/fc_comb_$1.tsv" "tmp/${SOURCE}/subread_aligned/$1.bam" -T 6 -t ${gene_types[@]} -g ${id_types[@]}
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
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "results/salmon/$1" ]; then
		mkdir results/salmon/$1
		salmon quant -i "salmon/${hERV_TRANSCRIPT%.*}_index" -l A -1 "tmp/${SOURCE}/qc/$1.qc.fq" -2 "tmp/${SOURCE}/qc/$2.qc.fq" --validateMappings -o "results/salmon/$1"
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
	timed_print "analysing..."

	if [[ ! -d "tmp" ]]; then 
		mkdir tmp
	fi

	if [[ ! -d "$ALIGN_METHOD" ]]; then
		mkdir $ALIGN_METHOD 
	fi

	if [[ ! -d "results" ]]; then 
		mkdir results
	fi

	if [ "$ANALYSIS_STEP" = "all" ]; then 
			ANALYSIS_STEP="index,convert,qc,align,count"
	fi 

	if [ $# = 1 ]; then 
		ANALYSIS_STEP=$1
	fi 

	CONVERT=true
	if [ $# = 2 ]; then 
		PAIR_FILE=$2
		CONVERT=false
	fi 

	ANALYSIS_STEP=(${ANALYSIS_STEP//,/ })
	for i in ${ANALYSIS_STEP[@]}; do
	  timed_print "$i-ing..."	
		case "$i" in 
			index) build_index ;;
			convert) if [ $CONVERT = true ]; then get_pairs_all; fi ;;
			qc) qc_all ;;
			align) align_all ;;
			count) count_all ;;
		esac 
		timed_print "finished $i"
	done 

	if [ $CLEAR_TMP = 'true' ]; then
		while IFS=, read -r r1 r2; do
			rm "tmp/${SOURCE}/$r1.fq"
			rm "tmp/${SOURCE}/$r2.fq"
			rm "tmp/${SOURCE}/qc/$r1.qc.fq"
			rm "tmp/${SOURCE}/qc/$r2.qc.fq"
		done < $PAIR_FILE 
	fi

	if [[ $# -ge 1 ]]; then 
		exit
	fi
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
