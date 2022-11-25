#!/bin/bash

cd $(dirname $0)

. ./lib-openfaas.sh

# Deploy function
deploy_fn nlp.yml || exit 1

# Verify Zstandard is installed
if ! which zstd > /dev/null; then
	sudo apt update && sudo apt install -y zstd
fi

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
