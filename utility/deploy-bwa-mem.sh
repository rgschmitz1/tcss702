#!/bin/bash

$(dirname $0)/pass-minio-secrets.sh

# Make bwa-mem bucket
bucket=minio/bwa-mem
if [ -z "$(mc ls $bucket)" ]; then
	mc mb $bucket

	# Move data into bucket for bwa-mem
	mc cp ../data/normal.tar.xz $bucket
	mc cp ../data/tumor.tar.xz $bucket
fi

cd $(dirname $0)/..

faas-cli up -f bwb-mem.yml
exit $?
