#!/bin/bash

. $(dirname $0)/lib-openfaas.sh

main -f nlp $@

# Install Zstandard
../utility/install-zstd.sh || exit $?

# Deploy minio
../utility/setup-minio.sh || exit 1
../utility/pass-minio-secrets.sh || exit 1

# Make buckets
bucket=minio/topic-modeling-us-east-1
mc ls $bucket 2> /dev/null || mc mb $bucket
mc ls ${bucket}-x86-64 2> /dev/null || mc mb ${bucket}-x86-64

# Move objects to input bucket
if [ -z "$(mc ls $bucket)" ]; then
	temp=$(mktemp -d)
	tar -xvf $PWD/../data/news_data.tar.zst -C $temp
	# Move data into bucket for pipeline
	if ! mc cp $temp/* $bucket; then
		prompt_error "failed to copy object(s) to $bucket"
		rm -fr $temp
		exit 1
	fi
	rm -fr $temp
fi

# Deploy function
deploy_fn || exit $?
