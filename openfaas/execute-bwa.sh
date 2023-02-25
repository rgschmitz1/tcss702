#!/bin/bash

# Support bash version 5.1 and newer only
if [ -z "$BASH_VERSION" ] || ([ ${BASH_VERSINFO[0]} -le 5 ] && [ ${BASH_VERSINFO[1]} -lt 1 ]); then
	printf "Run script using bash-5.1 or newer...\ncurrent version is $BASH_VERSION\n"
	exit 1
fi

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh

main -f bwa $@

for f in normal tumor; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "Executing $FN_NAME function, align $f sample"
		execute_fn "{\"inputfile\":\"$f.tar.zst\", \"bucket\":\"$FN_NAME\", \"thread_cnt\":4}" "${f}_$i"
	done
	check_concurrent_fn
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
