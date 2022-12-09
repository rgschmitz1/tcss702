#!/bin/bash

. $(dirname $0)/lib-openfaas.sh

# Remove function
if [[ -n "$1" && "$1" = '-d' ]]; then
	remove_fn nlp.yml
	exit $?
fi

# Install Zstandard
../utility/install-zstd.sh || exit $?

# Deploy minio
deploy_minio || exit 1

# Deploy function
deploy_fn nlp.yml || exit $?

# Make nlp buckets
bucket=minio/topic-modeling-us-east-1
mc ls $bucket 2> /dev/null || mc mb $bucket
mc ls ${bucket}-x86-64 2> /dev/null || mc mb ${bucket}-x86-64

# Move nlp objects to input bucket
if [ -z "$(mc ls $bucket)" ]; then
	# Move data into bucket for nlp pipeline
	temp=$(mktemp -d)
	tar -xvf $PWD/../data/news_data.tar.zst -C $temp
	if ! mc cp $temp/* $bucket; then
		prompt_error "failed to copy nlp object(s) to $bucket"
		rm -fr $temp
		exit 1
	fi
	rm -fr $temp
fi
