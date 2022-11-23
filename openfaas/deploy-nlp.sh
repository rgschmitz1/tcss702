#!/bin/bash

cd $(dirname $0)

../utility/setup-openfaas.sh
../utility/setup-minio.sh
../utility/install-docker.sh
../utility/pass-minio-secrets.sh || exit $?

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
		printf "\nERROR: failed to copy nlp object(s) to $bucket\n"
		rm -fr $temp
		exit 1
	fi
	rm -fr $temp
fi

# Build, push, and deploy nlp function
faas-cli up -f nlp.yml
exit $?
