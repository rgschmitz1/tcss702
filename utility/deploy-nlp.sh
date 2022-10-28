#!/bin/bash

$(dirname $0)/pass-minio-secrets.sh || exit $?

# Make nlp buckets
if [ -z "$(mc ls minio/topic-modeling-us-east-1)" ]; then
	mc mb minio/topic-modeling-us-east-1
	mc mb minio/topic-modeling-us-east-1-x86-64

	# Move data into bucket for nlp pipeline
	TEMP=$(mktemp -d)
	tar -xvf $(dirname $0)/../data/news_data.tar.xz -C $TEMP
	mc cp $TEMP/* minio/topic-modeling-us-east-1
	rm -fr $TEMP
fi

cd $(dirname $0)/..

faas-cli up -f nlp.yml
exit $?
