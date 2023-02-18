#!/bin/bash

. $(dirname $0)/lib-openfaas.sh

main -f sebs $@

# Deploy minio
../utility/setup-minio.sh || exit 1
../utility/pass-minio-secrets.sh || exit 1

# Make SeBS bucket
bucket=minio/$FN_NAME
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data to SeBS input bucket
key=bacillus_subtilis.fasta
if [ -z "$(mc ls $bucket/$key)" ] && ! mc cp ../data/$key $bucket; then
	prompt_error "failed to copy $key to $bucket"
	exit 1
fi

# Deploy function
deploy_fn || exit 1
