#!/bin/bash

# Support bash version 5.1 and newer only
if [ -z "$BASH_VERSION" ] || ([ ${BASH_VERSINFO[0]} -le 5 ] && [ ${BASH_VERSINFO[1]} -lt 1 ]); then
	printf "Run script using bash-5.1 or newer...\ncurrent version is $BASH_VERSION\n"
	# attempt to run using statically compiled version of bash for x86_64
	$(dirname $0)/../utility/bash-5.2.15 $0 $@
	exit $?
fi

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
