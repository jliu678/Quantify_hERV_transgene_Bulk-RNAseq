find . -type f -name "*.fastq.gz" -exec cp {} allFastq \;
find ./allFastq -type f -name "*.fastq.gz" -exec sh -c "mv -v $1 ${1%.fastq.gz}.fq.gz" _ {} \;

is_wrong=0; \
for i in salmon/* ; do
	if ! cat $i/quant.sf | grep "ENST"; then 
		is_wrong=1 
	fi 
done