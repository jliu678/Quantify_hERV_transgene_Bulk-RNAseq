#!/bin/bash
#main programme run file, bam -(samtools)> fastq -(fastp)> fastq -(ANY)> bam -(ANY)> counts

source main.sh

SOURCE=$(basename $SOURCE_LOC)
REF_GENOME=${basename $REF_GENOME_LOC .gz}

split_fastq(){ #splits pair ended fastq from bam into 2 files
	local bam_flags=$(samtools view -c -f 1 "${SOURCE_LOC}/$1.bam")
	local str_flags=$(samtools flags $bam_flags)
	if  [[ $str_flags == *"PAIRED"* ]]; then 
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/1$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r1.fq"
		cat "tmp/${SOURCE}/$1.fq" | grep '^@.*/2$' -A 3 --no-group-separator > \
				"tmp/${SOURCE}/$1.r2.fq"
		cat "tmp/${SOURCE}/$1.r1.fq\ttmp/${SOURCE}/$1.r2.fq\n" >> $PAIR_FILE
	fi
}

bam_to_fastq(){ #wrapper function, as the name suggests
	samtools bam2fq "${SOURCE_LOC}/$1.bam" > "tmp/${SOURCE}/$1.fq"
	split_fastq
}

group_fastq(){ #group fastq files into pairs
	files=(tmp/${SOURCE}/*.fq)
	for i in "${!files[@]}"; do #loop w/ index bcs easier
		local file_name=${files[$i]%.*} #get file name w/o extention
		if ! grep -Fxq "$file_name" $PAIR_FILE ; then #if the file does not have pair
			if [[ ${file_name: -1} = "1" && -f "${file_name::-1}2.fq" ]]; then #if formatted correctly
				cat "${files[$i]}\t${file_name::-1}2.fq\n" >> $PAIR_FILE
			else #the choice is yours how to deal with single-ended files
				echo "compliment to ${files[$i]} not found"
				#cat "${files[$i]}\n" >> $PAIR_FILE
			fi
		fi
	done 
}

get_pairs_all(){ #place all files into tmp, group them
	touch $PAIR_FILE
	for i in "${SOURCE_LOC}/*"; do
		local name=$(basename $i)
		case "${name#*.}" in #get extention
			bam) bam_to_fastq ${name%*.} ;; #handles grouping
			fq) mv "${SOURCE_LOC}/${name}" "tmp/${SOURCE}/${SEQ_NAME}.fq" #move because does not modify original data
			fq.gz) gunzip -c "${SOURCE_LOC}/${name}" > "tmp/${SOURCE}/${SEQ_NAME}.fq" #unzip for uniformity
		esac
	done 
	group_fastq #groups all other files
}

fastp_qc(){ #only works for pair ended as of now
	if [[ $# -eq 2 ]]; then 
		n1=${1%.*}; n2=${2%.*}
		fastp -i "${n1}.fq" -o "${n1}.qc.fq" \
					-I "${n2}.fq" -O "${n2}.qc.fq"
	fi
}

qc_all(){
	while IFS="\t" read -r r1 r2; do
		fastp_qc $r1 $r2
	done < $PAIR_FILE
}

#STAR --runThreadN 10 --runMode genomeGenerate --genomeDir hg38_index --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v41.chr_patch_hapl_scaff.annotation.gtf

#STAR --genomeDir hg38_index --runThreadN 8 --readFilesIn siyi_link/HcxecP3_1_r1.fastq siyi_link/HcxecP3_1_r2.fastq --outFileNamePrefix /siyi_link/results/HcxecP3_1_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard 

#star --runThreadN 8 --runMode alignReads --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts --outFileNamePrefix /siyi_link/results/HcxecP3_1_ --genomeDir hg38_index --sjdbGTFfile $gtfreference --readFilesIn siyi_link/HcxecP3_1_r1.fastq siyi_link/HcxecP3_1_r2.fastq

subread_build_index(){
	mkdir subread && cd subread
	subread-buildindex -o ${REF_GENOME%.*}_index $REF_GENOME
}

subread_align(){
	subread-align -i subread/${REF_GENOME%.*}_index -r "tmp/${SOURCE}/$1.qc.fq" -R "tmp/${SOURCE}/$2.qc.fq" -t 0 -o "tmp/${SOURCE}/$1.bam" -T 8 -M 16000
}

align_all(){
	while IFS="\t" read -r r1 r2; do
		subread_align $r1 $r2
	done < $PAIR_FILE
}

get_gene_types(){
	gene_types=$(cut -f3 ${COMB_ANNOTATION}.gff3 | grep -v ^# | sort | uniq)
	id_types=( "ID" "gene_id" )
}

subread_count(){
	featureCounts -a ${COMB_ANNOTATION}.gff3 -o subread/results/feature_count_$1.tsv subread/results/$1.bam -T 6 -t $gene_types -g $id_types
}

count_all{
	get_gene_types 
	for r1 in (cut -d, -f1 < in.csv); do 
		subread_align $r1
	done 
}

main(){
	ANALYSIS_STEP=(${ANALYSIS_STEP//,/})
	for i in ${ANALYSIS_STEP[@]}; do 
		case "$i" in 
			index) subread_build_index ;;
			convert) get_pairs_all ;;
			qc) qc_all ;;
			align) align_all ;;
			count) count all ;;
		esac 
	done 
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main "$@"
fi


#subread-align -i siyi_link/subread/hg38_index -r siyi_link/HcxecP3_1_r1.fastq -R siyi_link/HcxecP3_1_r2.fastq -t 0 -o siyi_link/subread/HcxecP3_1_aligned.bam -T 8 -M 16000

#featureCounts -a gencode.v41.chr_patch_hapl_scaff.annotation.gff3 -o siyi_link/subread/results/feature_count_uncombined.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6

#featureCounts -a erv.gff3 -o siyi_link/subread/results/feature_count_erv.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6 -t BED_feature,exon -g ID,gene_id

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt 
#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt
#cat package-entities-erv.gff3 gencode.v41.transcripts.fa GRCh38.p13.genome.fa > gentrome.fa
#salmon index -t gentrome.fa -d decoys.txt -p 6 -i index --gencode
#salmon quant -i siyi_link/salmon/index -l A -1 siyi_link/HcxecP3_1_r1.fastq -2 siyi_link/HcxecP3_1_r2.fastq --validateMappings -o siyi_link/salmon/results

#gffcompare -R -r gencode.v41.chr_patch_hapl_scaff.annotation.gtf -o siyi_link/combine_annotation/herv_annotation erv.gff3 >>warnings.txt 2>&1
#cuffmerge -g gencode.v41.chr_patch_hapl_scaff.annotation.gtf -s GRCh38.p13.genome.fa -o siyi_link/combine_annotation/herv_annotation gff_list.txt