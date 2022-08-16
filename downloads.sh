#!/bin/bash
#download the databases needed

. $main_loc/timed.sh

check_gencode(){
	timed_download "$REF_GENOME_LOC" "REF_GENOME"
	case "$ALIGN_METHOD"
		salmon)
			timed_download "$REF_TRANSCRIPT_LOC" "REF_TRANSCRIPT" ;;
		subread) 
			timed_download "$REF_ANNOTATION_LOC" "REF_ANNOTATION" ;;
	esac
	
}

download_hERV() {
	if [ "$ALIGN_METHOD" = "subread" ]; then
		hERV_download_list=(
			'https://herv.img.cas.cz/f/package-entities-erv.gff3.gz'
		)
	
		mkdir $hERV_DIR && cd $hERV_DIR
		timed_print "downloading hERVd..."
	
		for i in ${hERV_download_list[@]}; do
			timed_download $i
		done
	fi

	if [ "$ALIGN_METHOD" = "salmon" ]; then
		timed_download "$hERV_TRANSCRIPT_LOC" "hERV_TRANSCRIPT"
	fi

	cd ..
}

combine_hERV() { #combine into a sinlge file for ease of use 
  cat "$hERV_DIR"/*.gff3 > "$hERV_FILE"
}

check_hERV(){ 
	if [ ! -d $hERV_DIR ]; then 
	  timed_print "hERVd not downloaded"
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
