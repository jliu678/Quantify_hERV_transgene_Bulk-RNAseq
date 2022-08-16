#!/bin/bash
#download the databases needed

. $main_loc/timed.sh
hERV_TRANSCRIPT=$(basename $hERV_TRANSCRIPT_LOC .gz)

check_gencode(){
	timed_download "$REF_GENOME_LOC" 
	case "$ALIGN_METHOD" in 
		salmon)
			timed_download "$REF_TRANSCRIPT_LOC" ;;
		subread) 
			timed_download "$REF_ANNOTATION_LOC" ;;
	esac
	
}

combine_hERV() { #combine into a sinlge file for ease of use 
  cat "$hERV_DIR"/*.gff3 > "$hERV_FILE"
}

download_hERV_subread() {
	if [ ! -f $hERV_FILE ]; then
		hERV_download_list=(
			'https://herv.img.cas.cz/f/package-entities-erv.gff3.gz'
		)
	
		mkdir $hERV_DIR && cd $hERV_DIR
		timed_print "downloading hERVd..."
	
		for i in ${hERV_download_list[@]}; do
			timed_download $i
		done
		cd ..
	
		combine_hERV
	fi
}

download_hERV_salmon() {
	if [ ! -f $hERV_TRANSCRIPT ]; then
		timed_download "$hERV_TRANSCRIPT_LOC" 
	fi
}

download_mouse_ERV_salmon() {
	if [ ! -f $hERV_TRANSCRIPT ]; then
		timed_download "$hERV_TRANSCRIPT_LOC" 
	fi
	# gffread "$(basename $hERV_TRANSCRIPT_LOC .gz)" -g "$(basename $REF_GENOME_LOC .gz)" -w "$(basename $hERV_TRANSCRIPT_LOC .gtf .gff3 .gz).fa"
}

check_hERV(){ 
	# case "$ANIMAL_TYPE" in 
	# 	human) 
			case "$ALIGN_METHOD" in 
				subread) download_hERV_subread ;;
				salmon) download_hERV_salmon ;;
			esac 
	# 	mouse)
	# 		case "$ALIGN_METHOD" in 
	# 			salmon) download_mouse_ERV_salmon ;;
	# 		esac ;;
	# esac
}

combine_annotations(){
	cat $REF_ANNOTATION $hERV_FILE > $COMB_ANNOTATION
}

main(){
	check_gencode
	check_hERV
}

main
