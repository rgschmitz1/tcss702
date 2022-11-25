#!/bin/bash

cd $(dirname $0)

. ./lib-openfaas.sh

# Deploy function
deploy_fn bwa-mem.yml || exit 1

# Make bwa-mem bucket
bucket=minio/bwa-mem
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data into bwa-mem input bucket
for f in normal.tar.zst tumor.tar.zst; do
	[ -n "$(mc ls $bucket/$f)" ] && continue
	if [ ! -f "../data/$f" ]; then
		docker run --rm -v $PWD/../data:/data \
			schmitzr1984/tcss702-bwa-mem-reference cp /$f /data || exit 1
	fi
	if ! mc cp ../data/$f $bucket; then
		prompt_error "failed to copy $f to $bucket"
		exit 1
	fi
done
