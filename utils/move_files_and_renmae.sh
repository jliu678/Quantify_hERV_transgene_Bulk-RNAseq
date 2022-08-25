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
