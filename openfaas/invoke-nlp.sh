#!/bin/bash

cd $(dirname $0)

# Add library to invoke functions
. lib-openfaas.sh $1 "$2" 'nlp' || exit 1

for ((i=0; i<$ITERATION; i++)); do
	prompt_info "Iteration $i\n--"
	prompt_info "Invoking NLP preprocess function"
	execute_fn 'nlp' '{"fn":"p"}' "preprocess_$i"
	prompt_info "Invoking NLP training function"
	execute_fn 'nlp' '{"fn":"t"}' "training_$i"
	prompt_info "Invoking NLP query function"
	execute_fn 'nlp' '{"fn":"q"}' "query_$i"
done

printf "\nTotal iterations: $ITERATION, Total runtime: $SECONDS sec\n"
