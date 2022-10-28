#!/bin/bash

$(dirname $0)/pass-minio-secrets.sh || exit $?

# Make bwa-mem bucket
bucket=minio/bwa-mem
[ -z "$(mc ls $bucket)" ] && mc mb $bucket

# Move data into bucket for bwa-mem
for f in normal.tar.xz tumor.tar.xz; do
	[ -f ../data/$f ] || \
		docker run --rm -v ../data:/data schmitzr1984/tcss702-bwa-mem-reference cp /$f /data
	if ! (mc ls $bucket/$f || mc cp ../data/$f $bucket); then
		printf "\nERROR: failed to copy $f to $bucket\n"
		exit 1
	fi
done

cd $(dirname $0)/..

faas-cli up -f bwa-mem.yml
exit $?
