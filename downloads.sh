#!/bin/bash
#download the databases needed

. $main_loc/timed.sh
TRANSCRIPTS=()
for i in ${TRANSCRIPT_LOC[@]}; do TRANSCRIPTS+=( $(basename $i .gz) ); done
ANNOT_FILES=()
for i in ${ANNOT_LOCS[@]}; do ANNOT_FILES+=( $(basename $i .gz) ); done

check_gencode(){
	timed_download "$REF_GENOME_LOC" 
	case "$ALIGN_METHOD" in 
		salmon)
			timed_download "$REF_TRANSCRIPT_LOC" ;;
		subread) 
			timed_download "$REF_ANNOTATION_LOC" ;;
	esac
	
}

combine_annotations() { #combine into a sinlge file for ease of use 
	local annot_file_name="$(IFS=-; echo "${ANNOT_FILES[*]%.*}")"
  cat "$ANNOT_DIR"/*.gff3 > "$annot_file_name"
}

download_annotations_subread() {
	local annot_file_name="$(IFS=-; echo "${ANNOT_FILES[*]%.*}")"
	if [ ! -f "$annot_file_name" ]; then
	
		mkdir $ANNOT_DIR && cd $ANNOT_DIR
		timed_print "downloading hERVd..."
	
		for i in ${ANNOT_LOCS[@]}; do
			timed_download $i
		done
		cd ..
	
		combine_annotations
	fi
}

download_transcripts_salmon() {
	if [ ! -f $TRANSCRIPTS ]; then
		for i in ${TRANSCRIPT_LOC[@]}; do
			timed_download "$i" 
		done
	fi
}

# download_mouse_ERV_salmon() {
# 	if [ ! -f $hERV_TRANSCRIPT ]; then
# 		timed_download "$hERV_TRANSCRIPT_LOC" 
# 	fi
# 	# gffread "$(basename $hERV_TRANSCRIPT_LOC .gz)" -g "$(basename $REF_GENOME_LOC .gz)" -w "$(basename $hERV_TRANSCRIPT_LOC .gtf .gff3 .gz).fa"
# }

check_transcripts(){ 
	# case "$ANIMAL_TYPE" in 
	# 	human) 
			case "$ALIGN_METHOD" in 
				subread) download_annotation_subread ;;
				salmon) download_transcripts_salmon ;;
			esac 
	# 	mouse)
	# 		case "$ALIGN_METHOD" in 
	# 			salmon) download_mouse_ERV_salmon ;;
	# 		esac ;;
	# esac
}

main(){
	check_gencode
	check_transcripts
}

main
