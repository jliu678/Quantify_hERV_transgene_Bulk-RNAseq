#find . -type f -name "*.fastq.gz" -exec cp {} allFastq \;
#find ./allFastq -type f -name "*.fastq.gz" -exec sh -c "mv -v $1 ${1%.fastq.gz}.fq.gz" _ {} \;
#
#is_wrong=0; \
#for i in salmon/* ; do
#	if ! cat $i/quant.sf | grep "ENST"; then 
#		is_wrong=1 
#	fi 
#done

#the below dosen't work, please double quotes vs single quote
find ../raw_data -type f -name "*1.fq.gz" -exec sh -c "mv -v $1 ${1%1.fq.gz}R1.fq.gz" _ {} \;

#but the below work, please double quotes vs single quote
find ../raw_data -type f -name "*1.fq.gz" -exec sh -c 'mv -v $1 ${1%1.fq.gz}R1.fq.gz' _ {} \;
find ../raw_data -type f -name "*2.fq.gz" -exec sh -c 'mv -v $1 ${1%2.fq.gz}R2.fq.gz' _ {} \;

#a=( $(find ../raw_data/ -type f -name "*1.fq.gz") )
for i in ${a[@]}; do echo ${i%.*}; done
for i in ${a[@]}; do echo ${i%%.*}; done
for i in ${a[@]}; do echo ${i%1.*}R1.fq.gz; done
for i in ${a[@]}; do echo ${i%1.fq.gz}R1.fq.gz; done
for i in ${a[@]}; do echo ${i#*.}; done
for i in ${a[@]}; do echo ${i##*.}; done


#=======
find . -type f -name "*.fastq.gz" -exec cp {} allFastq \;
find ./allFastq -type f -name "*.fastq.gz" -exec sh -c "mv -v $1 ${1%.fastq.gz}.fq.gz" _ {} \;

is_wrong=0; \
for i in salmon/* ; do
	if ! cat $i/quant.sf | grep "ENST"; then 
		is_wrong=1 
	fi 
done

awk -F"\t" '{print $1, $2,"\t"$3,"\t"$4,"\t"$5,"\t"$6,"\t"$7,"\t"$8,"\t"$9" gene_type \"$3\";"}' Mmus38.geve.v1.gtf > Mmus38.geve.v1_mod.gtf

grep -v "^##" gencode.vM30.chr_patch_hapl_scaff.annotation.gtf | cut -f9 | sed -n -e 's/^.*gene_type //p' | cut -d';' -f1 | uniq 
>>>>>>> f0ecaa3f531cb9cdd329f5a786561b3d606af207
