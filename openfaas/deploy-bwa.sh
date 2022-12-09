#!/bin/bash

. $(dirname $0)/lib-openfaas.sh

# Remove function
if [[ -n "$1" && "$1" = '-d' ]]; then
	remove_fn bwa.yml
	exit $?
fi

# Deploy minio
deploy_minio || exit 1

# Deploy function
deploy_fn bwa.yml || exit 1

# Make bwa bucket
bucket=minio/bwa
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data into bwa input bucket
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
