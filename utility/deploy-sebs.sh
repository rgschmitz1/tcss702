#!/bin/bash

cd $(dirname $0)

./pass-minio-secrets.sh || exit $?

bucket=minio/sebs
key=bacillus_subtilis.fasta

# Make SeBS bucket
mc ls $bucket || mc mb $bucket

# Move data into bucket for SeBS
if [ -z "$(mc ls $bucket/$key)" ] && ! mc cp ../data/$key $bucket; then
	printf "\nERROR: failed to copy $key to $bucket\n"
	exit 1
fi

cd ../openfaas

faas-cli up -f bwa-mem.yml
exit $?
