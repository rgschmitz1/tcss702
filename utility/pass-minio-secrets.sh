#!/bin/bash

. $(dirname $0)/color-prompt.sh

ACCESSKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.accesskey}" | base64 --decode)
if [ -z "$ACCESSKEY" ]; then
	prompt_error "ERROR: ACCESSKEY variable is empty"
	exit 1
fi
SECRETKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.secretkey}" | base64 --decode)
if [ -z "$SECRETKEY" ]; then
	prompt_error "ERROR: SECRETKEY variable is empty"
	exit 1
fi
if ! (kubectl get secret -n openfaas-fn minio-access-key || \
	kubectl create secret -n openfaas-fn generic minio-access-key --from-literal minio-access-key="$ACCESSKEY"); then
	prompt_error "ERROR: failed to assign minio access key"
	exit 1
fi
if ! (kubectl get secret -n openfaas-fn minio-secret-key || \
	kubectl create secret -n openfaas-fn generic minio-secret-key --from-literal minio-secret-key="$SECRETKEY"); then
	prompt_error "ERROR: failed to assign minio secret key"
	exit 1
fi
