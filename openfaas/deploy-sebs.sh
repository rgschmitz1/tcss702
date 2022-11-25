#!/bin/bash

cd $(dirname $0)

. ./lib-openfaas.sh

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
