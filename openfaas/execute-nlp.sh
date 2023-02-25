#!/bin/bash

# Support bash version 5.1 and newer only
if [ -z "$BASH_VERSION" ] || ([ ${BASH_VERSINFO[0]} -le 5 ] && [ ${BASH_VERSINFO[1]} -lt 1 ]); then
	printf "Run script using bash-5.1 or newer...\ncurrent version is $BASH_VERSION\n"
	exit 1
fi

# Add library to invoke functions
. $(dirname $0)/lib-openfaas.sh

main -f nlp $@

# Create buckets
for ((i=0; i<$ITERATION; i++)); do
	mc mb minio/${FN_NAME}-${i}
done

for f in preprocess train query; do
	for ((i=0; i<$ITERATION; i++)); do
		prompt_info "Iteration $i\n--"
		prompt_info "Executing $FN_NAME, $f function"
		execute_fn "{\"fn\":\"${f}\", \"bucket\":\"${FN_NAME}-${i}\"}" "${f}_$i"
	done
	check_concurrent_fn
done

# Remove all created buckets
for ((i=0; i<$ITERATION; i++)); do
	mc rb --force minio/${FN_NAME}-${i}
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
