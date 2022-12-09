#!/bin/bash

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh

# if third optional argument is set to '-c' then concurrent function runs will take place
invoke_setup $1 "$2" 'bwa' "$3" || exit 1

for f in normal tumor; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "BWA function, align $f sample"
		execute_fn 'bwa' "{\"inputfile\":\"$f.tar.zst\", \"bucket\":\"bwa\", \"thread_cnt\":2}" "${f}_$i"
	done
	check_concurrent_fn
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
