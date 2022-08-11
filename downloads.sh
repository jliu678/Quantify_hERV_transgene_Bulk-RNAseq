#!/bin/bash
#download the databases needed

. $main_loc/timed.sh

check_gencode(){
	timed_download "$REF_ANNOTATION_LOC" "REF_ANNOTATION"
	timed_download "$REF_ANNOTATION_LOC" "REF_GENOME"
}

download_hERV() {
	hERV_download_list=(
		'https://herv.img.cas.cz/f/package-entities-dna.gff3.gz'
		'https://herv.img.cas.cz/f/package-entities-erv.gff3.gz'
		'https://herv.img.cas.cz/f/package-entities-line.gff3.gz'
		'https://herv.img.cas.cz/f/package-entities-rc.gff3.gz'
	)

	mkdir $hERV_DIR && cd $hERV_DIR
	timed_print "downloading hERVd..."

	for i in ${hERV_download_list[@]}; do
		timed_download $i
	done

	cd ..
}

combine_hERV() { #combine into a sinlge file for ease of use 
  cat "$hERV_DIR"/*.gff3 > "$hERV_FILE".gff3
}

check_hERV(){ 
	if [ -d $hERV_DIR ]; then 
	  echo "hERVd not downloaded"
		download_hERV 
		combine_hERV 
	fi
}

combine_annotations(){
	cat $REF_ANNOTATION $hERV_FILE > $COMB_ANNOTATION
}

main(){
	check_gencode
	check_hERV
}

main
