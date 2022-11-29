#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh

invoke_setup $1 "$2" 'bwa' "$3" || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	for f in normal tumor; do
		prompt_info "BWA function, align $f sample"
		execute_fn 'bwa' "{\"inputfile\":\"$f.tar.zst\", \"bucket\":\"bwa\"}" "${f}_$i" $CONCURRENT
	done
done

check_concurrent_fn

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
