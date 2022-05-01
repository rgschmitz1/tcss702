#!/bin/bash

# Get the access and secret key to gain access to minio
ACCESSKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.accesskey}" | base64 --decode)
SECRETKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.secretkey}" | base64 --decode)

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
sed -i "s|\(.*access_key: \).*|\1${ACCESSKEY}|" nlp.yml
sed -i "s|\(.*secret_key: \).*|\1${SECRETKEY}|" nlp.yml

faas-cli up -f nlp.yml
exit $?
