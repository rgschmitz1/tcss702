#!/bin/bash

# Install minio with arkade
arkade install minio
# TODO: figure out persistant storage
#    --set persistence.enabled=true
    
# Forward the minio port to your machine
kubectl rollout status -n default deploy/minio
pgrep -f kubectl.*9000 && (pkill -f kubectl.*9000 > /dev/null; sleep 3)
kubectl port-forward -n default svc/minio 9000:9000 &
if [ $? -ne 0 ]; then
	echo "Failed to port-forward minio service"
	exit 1
fi
sleep 1

# Get the access and secret key to gain access to minio
ACCESSKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.accesskey}" | base64 --decode)
SECRETKEY=$(kubectl get secret -n default minio -o jsonpath="{.data.secretkey}" | base64 --decode)

# Add a host
mc config host add minio http://localhost:9000 $ACCESSKEY $SECRETKEY
ret=$?
[ $ret -eq 0 ] && sleep 1
exit $ret
