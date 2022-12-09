#!/bin/bash

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh bwa $@

for f in normal tumor; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "Executing $FN_NAME function, align $f sample"
		execute_fn "{\"inputfile\":\"$f.tar.zst\", \"bucket\":\"$FN_NAME\", \"thread_cnt\":2}" "${f}_$i"
	done
	check_concurrent_fn
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
