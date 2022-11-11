#!/bin/bash

cd $(dirname $0)

../utility/pass-minio-secrets.sh || exit $?

# Make SeBS bucket
bucket=minio/sebs
mc ls $bucket 2> /dev/null || mc mb $bucket

# Move sample data to SeBS input bucket
key=bacillus_subtilis.fasta
if [ -z "$(mc ls $bucket/$key)" ] && ! mc cp ../data/$key $bucket; then
	printf "\nERROR: failed to copy $key to $bucket\n"
	exit 1
fi

# Build, push, and deploy SeBS function
faas-cli up -f sebs.yml
exit $?
