#!/bin/bash

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh nlp $@

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	for f in preprocess train query; do
		prompt_info "Executing $FN_NAME, $f function"
		execute_fn "{\"fn\":\"${f}\"}" "${f}_$i"
	done
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
