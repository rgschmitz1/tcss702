#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh $1 "$2" 'sebs' || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "SeBS function, dna_visualization"
	execute_fn 'sebs' '{"fn":"dna_visualization", "key":"bacillus_subtilis.fasta", "bucket":"sebs"}' "dna_visualization_$i"
	prompt_info "SeBS function, graph_bfs"
	execute_fn 'sebs' '{"fn":"graph_bfs", "size":10000}' "graph_bfs_$i"
	prompt_info "SeBS function, graph_mst"
	execute_fn 'sebs' '{"fn":"graph_mst", "size":10000}' "graph_mst_$i"
	prompt_info "SeBS function, graph_pagerank"
	execute_fn 'sebs' '{"fn":"graph_pagerank", "size":10000}' "graph_pagerank_$i"
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
