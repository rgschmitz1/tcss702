#!/bin/bash

$(dirname $0)/pass-minio-secrets.sh

# Make bwa-mem bucket
bucket=minio/bwa-mem
[ -z "$(mc ls $bucket)" ] && mc mb $bucket

# Move data into bucket for bwa-mem
for f in normal tumor; do
	[ -f ../data/$f.tar.xz ] || \
		docker run --rm -v ../data:/data schmitzr1984/tcss702-bwa-mem-reference cp /$f.tar.xz /data
	mc cp ../data/$f.tar.xz $bucket
done

cd $(dirname $0)/..

faas-cli up -f bwb-mem.yml
exit $?
