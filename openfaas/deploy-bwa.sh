#!/bin/bash

. $(dirname $0)/lib-openfaas.sh bwa $@

# Deploy minio
deploy_minio || exit 1

# Make bucket
bucket=minio/$FN_NAME
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data into input bucket
pushd ../data > /dev/null
for f in normal.tar.zst tumor.tar.zst; do
	[ -n "$(mc ls $bucket/$f)" ] && continue
	if [ ! -f "$f" ]; then
		docker run --rm -v $PWD:/data \
			schmitzr1984/tcss702-bwa-reference \
			cp -v /$f /data || exit 1
	fi
	if ! mc cp $f $bucket; then
		prompt_error "failed to copy $f to $bucket"
		exit 1
	fi
done
popd > /dev/null

# Deploy function
deploy_fn || exit 1
