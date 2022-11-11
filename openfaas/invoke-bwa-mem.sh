#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh $1 "$2" 'bwa-mem' || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "BWA-MEM function, align normal sample"
	execute_fn 'bwa-mem' '{"inputfile":"normal.tar.zst", "bucket":"bwa-mem"}' "normal_$i"
	prompt_info "BWA-MEM function, align tumor sample"
	execute_fn 'bwa-mem' '{"inputfile":"tumor.tar.zst", "bucket":"bwa-mem"}' "tumor_$i"
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
