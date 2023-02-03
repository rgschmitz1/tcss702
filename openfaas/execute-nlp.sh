#!/bin/bash

# Support bash version 5.1 and newer only
if [ -z "$BASH_VERSION" ] || ([ ${BASH_VERSINFO[0]} -le 5 ] && [ ${BASH_VERSINFO[1]} -lt 1 ]); then
	printf "Run script using bash-5.1 or newer...\ncurrent version is $BASH_VERSION\n"
	# attempt to run using statically compiled version of bash for x86_64
	$(dirname $0)/../utility/bash-5.2.15 $0 $@
	exit $?
fi

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
