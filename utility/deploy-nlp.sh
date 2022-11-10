#!/bin/bash

cd $(dirname $0)

./pass-minio-secrets.sh || exit $?

if ! which zstd 2> /dev/null; then
	sudo apt update && sudo apt install -y zstd
fi

# Make nlp buckets
if [ -z "$(mc ls minio/topic-modeling-us-east-1)" ]; then
	mc mb minio/topic-modeling-us-east-1
	mc mb minio/topic-modeling-us-east-1-x86-64

	# Move data into bucket for nlp pipeline
	TEMP=$(mktemp -d)
	tar -xvf $(dirname $0)/../data/news_data.tar.zst -C $TEMP
	mc cp $TEMP/* minio/topic-modeling-us-east-1
	rm -fr $TEMP
fi

cd ../openfaas

faas-cli up -f nlp.yml
exit $?
