#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh

invoke_setup $1 "$2" 'nlp' || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	for f in preprocess train query; do
		prompt_info "Invoking NLP $f function"
		execute_fn 'nlp' "{\"fn\":\"${f}\"}" "${f}_$i"
	done
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
