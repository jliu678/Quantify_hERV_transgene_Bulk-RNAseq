#!/bin/bash
#main programme run file, bam -(samtools)> fastq -(fastp)> fastq -(ANY)> bam -(ANY)> counts

. $main_loc/timed.sh

SOURCE=$(basename $SOURCE_LOC)
REF_GENOME=$(basename $REF_GENOME_LOC .gz)
REF_ANNOTATION=$(basename $REF_ANNOTATION_LOC .gz)
id_types=( "ID" )

split_fastq(){ #splits pair ended fastq from bam into 2 files
	local bam_flags=$(samtools view -c -f 1 "${SOURCE_LOC}/$1.bam")
	local str_flags=$(samtools flags $bam_flags)
	if  [[ $str_flags == *"PAIRED"* ]]; then 
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/1$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r1.fq"
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/2$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r2.fq"
		cat "tmp/${SOURCE}/$1.r1\ttmp/${SOURCE}/$1.r2\n" >> $PAIR_FILE
	fi
}

bam_to_fastq(){ #wrapper, function as the name suggests
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/$1.r1.fq" ] || [ ! -f "tmp/${SOURCE}/$1.r2.fq" ]; then
		samtools bam2fq "${SOURCE_LOC}/$1.bam" > "tmp/${SOURCE}/$1.fq"
		split_fastq
	fi
}

group_fastq(){ #group fastq files into pairs
	files=(tmp/${SOURCE}/*.fq)
	local total=0
	for i in "${!files[@]}"; do #loop w/ index bcs easier
		local file_name=${files[$i]%.*} #get file name w/o extention
		if ! grep -Fxq "$file_name" $PAIR_FILE ; then #if the file does not have pair
			if [[ ${file_name: -1} = "1" && -f "${file_name::-1}2.fq" ]]; then #if formatted correctly
				cat "${file_name}\t${file_name::-1}2\n" >> $PAIR_FILE
			else #the choice is yours how to deal with single-ended files
				timed_print "compliment to ${files[$i]} not found"
				total+=1
				#cat "${files[$i]}\n" >> $PAIR_FILE
			fi
		fi
	done
	timed_print "found total $total files without compliment"
}

mv_fq() {
	if [ ! "$OVER_WRITE" = "true" ] || [ ! -f "$2" ]; then
		case "${name#*.}" in 
			fq) mv $1 $2 ;;
			fq.gz) gunzip -c $1 > $2 ;;
		esac
	fi
}

get_pairs_all() { #place all files into tmp, group them
	touch $PAIR_FILE
	for i in "${SOURCE_LOC}/*"; do
		local name=$(basename $i)
		case "${name#*.}" in #get extention
			bam) bam_to_fastq ${name%*.} ;; #handles grouping
			fq) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/${SEQ_NAME}.fq" ;; #move because does not modify original data
			fq.gz) mv_fq "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/${SEQ_NAME}.fq" ;; #unzip for uniformity
		esac
	done 
	group_fastq #groups all other files
}

fastp_qc(){ #only works for pair ended as of now
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "$1.qc.fq" ]; then
		if [[ $# -eq 2 ]]; then 
			fastp -i "$1.fq" -o "$1.qc.fq" \
						-I "$2.fq" -O "$2.qc.fq" \
						-j "$1.json" -h "$1.html"
		fi
	fi
}

qc_all(){
	while IFS="\t" read -r r1 r2; do
		fastp_qc $r1 $r2
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

build_index() {
	subread_build_index
}

subread_align(){
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/$1.bam" ]; then
		timed_print "aligning $1.qc.fq and $2.qc.fq..."
		subread-align -i subread/${REF_GENOME%.*}_index -r "tmp/${SOURCE}/$1.qc.fq" -R "tmp/${SOURCE}/$2.qc.fq" -t 0 -o "tmp/${SOURCE}/$1.bam" -T 8 -M 16000
		timed_print "aligned $1.qc.fq and $2.qc.fq"
	fi
}

align_all(){
	while IFS="\t" read -r r1 r2; do
		subread_align $r1 $r2
	done < $PAIR_FILE
}

get_gene_types_comb(){
	gene_types=$(cut -f3 ${COMB_ANNOTATION}.gff3 | grep -v ^# | sort | uniq)
}

subread_count(){
	if [ ! "$OVER_WRITE" = "true" ] || [ ! -f "subread/results/fc_comb_$1.tsv" ]; then 
		timed_print "counting $1.bam..."
		featureCounts -a ${COMB_ANNOTATION}.gff3 -o subread/results/fc_comb_$1.tsv subread/results/$1.bam -T 6 -t ${gene_types[@]} -g ${id_types[@]}
		timed_print "counted $1.bam"
	fi
}

get_gene_types_sep(){
	ref_gene_types=$(cut -f3 ${REF_ANNOTATION}.gff3 | grep -v ^# | sort | uniq)
	erv_gene_types=$(cut -f3 ${hERV_FILE}.gff3 | grep -v ^# | sort | uniq)
}

subread_count_sep(){
	if [ ! "$OVER_WRITE" = "true" ] || [ ! -f "subread/results/fc_ref_$1.tsv" ]; then 
		timed_print "counting $1.bam with ${REF_ANNOTATION}.gff3..."
		featureCounts -a ${REF_ANNOTATION}.gff3 -o subread/results/fc_ref_$1.tsv subread/results/$1.bam -T 6 -t ${ref_gene_types[@]} -g ${id_types[@]}
		timed_print "counted $1.bam"
	fi 

	if [ ! "$OVER_WRITE" = "true" ] || [ ! -f "subread/results/fc_erv_$1.tsv" ]; then 
		timed_print "counting $1.bam with ${hERV_FILE}.gff3..."
		featureCounts -a ${hERV_FILE}.gff3 -o subread/results/fc_erv_$1.tsv subread/results/$1.bam -T 6 -t ${erv_gene_types[@]} -g ${id_types[@]}
		timed_print "counted $1.bam"
	fi
}

count_all(){
	if [ $COUNT_METHOD = "combined" ]; then
		get_gene_types 
		for r1 in $(cut -d, -f1 < ${PAIR_FILE}); do 
			subread_count $r1
		done 
	elif [ $COUNT_METHOD = "seperated" ]; then 
		get_gene_types_sep
		for r1 in $(cut -d, -f1 < ${PAIR_FILE}); do 
			subread_count_sep $r1
		done
	fi
}



main(){
	if [[ -d "tmp" ]]; then 
		mkdir tmp
	fi

	if [[ ! -d "subread" ]]; then
		mkdir subread 
	fi

	if [[ "$ANALYSIS_STEP" -eq "all" ]]; then 
			ANALYSIS_STEP="index,convert,qc,align,count"
	fi 

	if [ $# = 1 ]; then 
		PAIR_FILE=$1
	fi 

	ANALYSIS_STEP=(${ANALYSIS_STEP//,/})
	for i in ${ANALYSIS_STEP[@]}; do 
		case "$i" in 
			index) build_index ;;
			convert) get_pairs_all ;;
			qc) qc_all ;;
			align) align_all ;;
			count) count all ;;
		esac 
	done 
}

main


#subread-align -i siyi_link/subread/hg38_index -r siyi_link/HcxecP3_1_r1.fastq -R siyi_link/HcxecP3_1_r2.fastq -t 0 -o siyi_link/subread/HcxecP3_1_aligned.bam -T 8 -M 16000

#featureCounts -a gencode.v41.chr_patch_hapl_scaff.annotation.gff3 -o siyi_link/subread/results/feature_count_uncombined.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6

#featureCounts -a erv.gff3 -o siyi_link/subread/results/feature_count_erv.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6 -t BED_feature,exon -g ID

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt 
#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt
#cat package-entities-erv.gff3 gencode.v41.transcripts.fa GRCh38.p13.genome.fa > gentrome.fa
#salmon index -t gentrome.fa -d decoys.txt -p 6 -i index --gencode
#salmon quant -i index -l A -1 HcxecP3_1_r1.fastq -2 HcxecP3_1_r2.fastq --validateMappings -o results

#gffcompare -R -r gencode.v41.chr_patch_hapl_scaff.annotation.gtf -o siyi_link/combine_annotation/herv_annotation erv.gff3 >>warnings.txt 2>&1
#cuffmerge -g gencode.v41.chr_patch_hapl_scaff.annotation.gtf -s GRCh38.p13.genome.fa -o siyi_link/combine_annotation/herv_annotation gff_list.txt
