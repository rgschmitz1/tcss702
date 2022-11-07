#!/bin/bash
cd $(dirname $0)
. color-prompt.sh

execute_sebs_fn() {
	if curl -H "Content-Type: application/json" -X POST -d "$1" http://localhost:8080/function/sebs; then
		sleep 2
	else
		prompt_error "Failed to execute $1"
		exit 1
	fi
}

prompt_info "SeBS function, dna_visualization"
execute_sebs_fn '{"fn":"dna_visualization", "key":"bacillus_subtilis.fasta", "bucket":"sebs"}'
prompt_info "SeBS function, graph_bfs"
execute_sebs_fn '{"fn":"graph_bfs", "size":10000}'
prompt_info "SeBS function, graph_mst"
execute_sebs_fn '{"fn":"graph_mst", "size":10000}'
prompt_info "SeBS function, graph_pagerank"
execute_sebs_fn '{"fn":"graph_pagerank", "size":10000}'

printf "\nTotal runtime: $SECONDS\n"
