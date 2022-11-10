#!/bin/bash
cd $(dirname $0)
. color-prompt.sh

execute_bwa-mem_fn() {
	if curl -H "Content-Type: application/json" -X POST -d "$1" http://localhost:8080/function/bwa-mem; then
		sleep 2
	else
		prompt_error "Failed to execute $1"
		exit 1
	fi
}

prompt_info "BWA-MEM function, align normal sample"
execute_bwa-mem_fn '{"inputfile":"normal.tar.zst", "bucket":"bwa-mem"}'
prompt_info "BWA-MEM function, align tumor sample"
execute_bwa-mem_fn '{"inputfile":"tumor.tar.zst", "bucket":"bwa-mem"}'

printf "\nTotal runtime: $SECONDS\n"
