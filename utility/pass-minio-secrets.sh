#!/bin/bash
ACCESSKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.accesskey}" | base64 --decode)
SECRETKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.secretkey}" | base64 --decode)
kubectl create secret generic minio-access-key --from-literal minio-access-key="$ACCESSKEY" --namespace openfaas-fn
kubectl create secret generic minio-secret-key --from-literal minio-secret-key="$SECRETKEY" --namespace openfaas-fn
