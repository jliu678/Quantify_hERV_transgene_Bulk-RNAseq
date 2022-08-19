find . -type f -name "*.fastq.gz" -exec cp {} allFastq \;
find ./allFastq -type f -name "*.fastq.gz" -exec sh -c "mv -v $1 ${1%.fastq.gz}.fq.gz" _ {} \;