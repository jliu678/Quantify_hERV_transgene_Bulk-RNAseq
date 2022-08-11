#!/bin/bash

json_2_query() {
  echo "$*" | jq -sRr @uri
}

query_TCGA_seq() { # inputs: # of queries
	local TCGA_generic_json=$(envsubst <../query.json)
  local TCGA_query=$(json_2_query $TCGA_generic_json)
  curl 'https://api.gdc.cancer.gov/files?filters='"${TCGA_query}"'&fields=file_id,data_type,experimental_strategy,cases.case_id&pretty=true&size='"$1" > query_output.json
}

get_seq_IDs(){ # no inputs
	jq '[.data.hits[] | {"exp": .experimental_strategy, "id": {"file_id": .id, "case_id": .cases[0].case_id}}] | reduce .[] as $d (null; .[$d.exp] += [$d.id])' query_output.json > query_ids.json
}

get_TCGA_post() {
	jq -r '.'"$1"'[] | .file_id' query_ids.json > post.txt
	sed -i -e 's/^/ids=/' post.txt
	tr '\n' '&' < post.txt > "${1}_post.txt"
	rm post.txt
}

download_TCGA(){ #input: location of post.txt
	curl --remote-name --remote-header-name --request POST 'https://api.gdc.cancer.gov/data' --data @"${1}_post.txt" --create-dirs -O --output "TCGA/${1}"
}

main(){
	mkdir "tcga" && cd tcga 
	query_TCGA_seq $NUM_FILES
	get_seq_IDs 
	for i in ${SEQ_TYPE[@]}; do
		get_TCGA_post $i
		download_TCGA $i
	done 
	cd ..
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
  main "$@"
fi