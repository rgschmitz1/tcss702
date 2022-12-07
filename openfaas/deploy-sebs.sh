#!/bin/bash

cd $(dirname $0)

. ./lib-openfaas.sh

# Remove function
if [[ -n "$1" && "$1" = '-d' ]]; then
	remove_fn sebs.yml
	exit $?
fi

# Deploy minio
deploy_minio || exit 1

# Deploy function
deploy_fn sebs.yml || exit 1

# Make SeBS bucket
bucket=minio/sebs
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data to SeBS input bucket
key=bacillus_subtilis.fasta
if [ -z "$(mc ls $bucket/$key)" ] && ! mc cp ../data/$key $bucket; then
	prompt_error "failed to copy $key to $bucket"
	exit 1
fi
