#!/bin/bash

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh sebs $@

f=dna_visualization
for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	prompt_info "Executing $FN_NAME, $f"
	execute_fn "{\"fn\":\"$f\", \"key\":\"bacillus_subtilis.fasta\", \"bucket\":\"$FN_NAME\"}" "${f}_$i"
done
check_concurrent_fn

for f in graph_bfs graph_mst graph_pagerank; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "Executing $FN_NAME, $f"
		execute_fn "{\"fn\":\"$f\", \"size\":10000}" "${f}_$i"
	done
	check_concurrent_fn
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
