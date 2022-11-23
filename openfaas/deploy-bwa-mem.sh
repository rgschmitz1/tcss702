#!/bin/bash

cd $(dirname $0)

../utility/setup-openfaas.sh
../utility/setup-minio.sh
../utility/install-docker.sh
../utility/pass-minio-secrets.sh || exit $?

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
		printf "\nERROR: failed to copy $f to $bucket\n"
		exit 1
	fi
done

# Build, push, and deploy bwa-mem function
faas-cli up -f bwa-mem.yml
exit $?
