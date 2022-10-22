#!/bin/bash
cd $(dirname $0)
. color-prompt.sh

execute_nlp_fn() {
	if curl -H "Content-Type: application/json" -X POST -d "$1" http://localhost:8080/function/nlp; then
		sleep 2
	else
		prompt_error "Failed to execute $1"
		exit 1
	fi
}

prompt_info "Invoking NLP preprocess function"
execute_nlp_fn '{"fn":"p"}'
prompt_info "Invoking NLP training function"
execute_nlp_fn '{"fn":"t"}'
prompt_info "Invoking NLP query function"
execute_nlp_fn '{"fn":"q"}'

printf "\nTotal runtime: $SECONDS\n"
