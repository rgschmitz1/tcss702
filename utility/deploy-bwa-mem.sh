#!/bin/bash

cd $(dirname $0)

./pass-minio-secrets.sh || exit $?

# Make bwa-mem bucket
bucket=minio/bwa-mem
mc ls $bucket || mc mb $bucket

# Move data into bucket for bwa-mem
for f in normal.tar.xz tumor.tar.xz; do
	[ -n "$(mc ls $bucket/$f)" ] && continue
	if [ ! -f "../data/$f" ]; then
		docker run --rm -v ../data:/data schmitzr1984/tcss702-bwa-mem-reference cp /$f /data
	fi
	if ! mc cp ../data/$f $bucket; then
		printf "\nERROR: failed to copy $f to $bucket\n"
		exit 1
	fi
done

cd ..

faas-cli up -f bwa-mem.yml
exit $?
