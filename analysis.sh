#!/bin/bash
#main programme run file, bam -(samtools)> fastq -(fastp)> fastq -(ANY)> bam -(ANY)> counts

source main.sh

split_fastq(){
	local bam_flags=$(samtools view -c -f 1 ${SEQ_NAME}.bam)
	local str_flags=$(samtools flags $bam_flags)
	if  [[ $str_flags == *"PAIRED"* ]]; then 
		cat "tmp/${CUR_SEQ}/${SEQ_NAME}.fastq" | grep '^@.*/1$' -A 3 --no-group-separator > \
				"tmp/${CUR_SEQ}/${SEQ_NAME}_r1.fastq"
		cat "tmp/${CUR_SEQ}/${SEQ_NAME}.fastq" | grep '^@.*/2$' -A 3 --no-group-separator > \
				"tmp/${CUR_SEQ}/${SEQ_NAME}_r2.fastq"
		SEQ_NAME=("${SEQ_NAME}_r1" "${SEQ_NAME}_r2")
	fi
}

bam_to_fastq(){
	samtools bam2fq "${SOURCE}/${CUR_SEQ}/${SEQ_NAME}.bam" > "tmp/${CUR_SEQ}/${SEQ_NAME}.fastq"
	split_fastq
}

check_dir(){
		for i in 
}

fastp_qc(){
	if [ ${#SEQ_NAME[@]} = 2]; then 
		fastp -i "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.fastq" \
					-o "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.qc.fastq" \
					-I "tmp/${CUR_SEQ}/${SEQ_NAME[1]}.fastq" \
					-O "tmp/${CUR_SEQ}/${SEQ_NAME[1]}.qc.fastq"
	else 
		fastp -i "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.fastq" \
					-o "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.qc.fastq" 
	fi 
}

#STAR --runThreadN 10 --runMode genomeGenerate --genomeDir hg38_index --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v41.chr_patch_hapl_scaff.annotation.gtf

#STAR --genomeDir hg38_index --runThreadN 8 --readFilesIn siyi_link/HcxecP3_1_r1.fastq siyi_link/HcxecP3_1_r2.fastq --outFileNamePrefix /siyi_link/results/HcxecP3_1_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard 

#star --runThreadN 8 --runMode alignReads --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts --outFileNamePrefix /siyi_link/results/HcxecP3_1_ --genomeDir hg38_index --sjdbGTFfile $gtfreference --readFilesIn siyi_link/HcxecP3_1_r1.fastq siyi_link/HcxecP3_1_r2.fastq

subread_build_index(){
	mkdir subread && cd subread
	subread-buildindex -o ${REF_GENOME%.*}_index $REF_GENOME
}

subread_align(){
	subread-align -i subread/${REF_GENOME%.*}_index -r "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.qc.fastq" -R "tmp/${CUR_SEQ}/${SEQ_NAME[1]}.qc.fastq" -t 0 -o "tmp/${CUR_SEQ}/${SEQ_NAME[0]}.bam" -T 8 -M 16000
}

get_gene_types(){
	gene_types=$(cut -f3 ${COMB_ANNOTATION}.gff3 | grep -v ^# | sort | uniq)
	id_types=( "ID" "gene_id" )
}

subread_count(){
	featureCounts -a ${COMB_ANNOTATION}.gff3 -o subread/results/feature_count_${SEQ_NAME}.tsv subread/results/${SEQ_NAME}.bam -T 6 -t $gene_types -g $id_types
}

main(){
	
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main "$@"
fi


#subread-align -i siyi_link/subread/hg38_index -r siyi_link/HcxecP3_1_r1.fastq -R siyi_link/HcxecP3_1_r2.fastq -t 0 -o siyi_link/subread/HcxecP3_1_aligned.bam -T 8 -M 16000

#featureCounts -a gencode.v41.chr_patch_hapl_scaff.annotation.gff3 -o siyi_link/subread/results/feature_count_uncombined.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6

#featureCounts -a erv.gff3 -o siyi_link/subread/results/feature_count_erv.tsv siyi_link/subread/HcxecP3_1_aligned.bam -T 6 -t BED_feature,exon -g ID,gene_id

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

#grep "^>" GRCh38.p13.genome.fa | cut -d " " -f 1 > siyi_link/salmon/decoys.txt 
#sed -i .bak -e 's/>//g' siyi_link/salmon/decoys.txt
#cat gencode.v41.transcripts.fa GRCh38.p13.genome.fa > siyi_link/salmon/gentrome.fa
#salmon index -t siyi_link/salmon/gentrome.fa -d siyi_link/salmon/decoys.txt -p 6 -i siyi_link/salmon/index --gencode
#salmon quant -i siyi_link/salmon/index -l A -1 siyi_link/HcxecP3_1_r1.fastq -2 siyi_link/HcxecP3_1_r2.fastq --validateMappings -o siyi_link/salmon/results

#gffcompare -R -r gencode.v41.chr_patch_hapl_scaff.annotation.gtf -o siyi_link/combine_annotation/herv_annotation erv.gff3 >>warnings.txt 2>&1
#cuffmerge -g gencode.v41.chr_patch_hapl_scaff.annotation.gtf -s GRCh38.p13.genome.fa -o siyi_link/combine_annotation/herv_annotation gff_list.txt