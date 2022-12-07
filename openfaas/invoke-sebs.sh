#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh

invoke_setup $1 "$2" 'sebs' "$3" || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	prompt_info "SeBS function, dna_visualization"
	execute_fn 'sebs' '{"fn":"dna_visualization", "key":"bacillus_subtilis.fasta", "bucket":"sebs"}' "dna_visualization_$i"
done
check_concurrent_fn

for f in graph_bfs graph_mst graph_pagerank; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "SeBS function, $f"
		execute_fn 'sebs' "{\"fn\":\"$f\", \"size\":10000}" "${f}_$i"
	done
	check_concurrent_fn
done


printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
