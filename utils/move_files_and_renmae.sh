#find . -type f -name "*.fastq.gz" -exec cp {} allFastq \;
#find ./allFastq -type f -name "*.fastq.gz" -exec sh -c "mv -v $1 ${1%.fastq.gz}.fq.gz" _ {} \;
#
#is_wrong=0; \
#for i in salmon/* ; do
#	if ! cat $i/quant.sf | grep "ENST"; then 
#		is_wrong=1 
#	fi 
#done

#the below dosen't work
find ../raw_data -type f -name "*1.fq.gz" -exec sh -c "mv -v $1 ${1%1.fq.gz}R1.fq.gz" _ {} \;

#but the below work
find ../raw_data -type f -name "*1.fq.gz" -exec sh -c 'mv -v $1 ${1%1.fq.gz}R1.fq.gz' _ {} \;
find ../raw_data -type f -name "*2.fq.gz" -exec sh -c 'mv -v $1 ${1%2.fq.gz}R2.fq.gz' _ {} \;

#a=( $(find ../raw_data/ -type f -name "*1.fq.gz") )
for i in ${a[@]}; do echo ${i%.*}; done
for i in ${a[@]}; do echo ${i%%.*}; done
for i in ${a[@]}; do echo ${i%1.*}R1.fq.gz; done
for i in ${a[@]}; do echo ${i%1.fq.gz}R1.fq.gz; done
for i in ${a[@]}; do echo ${i#*.}; done
for i in ${a[@]}; do echo ${i##*.}; done


