#!/bin/bash
cd $(dirname $0)
. color-prompt.sh

execute_bwa-mem_fn() {
	if curl -H "Content-Type: application/json" -X POST http://localhost:8080/function/bwa-mem; then
		sleep 2
	else
		prompt_error "Failed to execute $1"
		exit 1
	fi
}

prompt_info "Invoking BWA-MEM function"
execute_bwa-mem_fn

printf "\nTotal runtime: $SECONDS\n"
