#!/bin/bash

cd $(dirname $0)

. color-prompt.sh

# Install dependencies
./install-arkade_helm_faas-cli_mc.sh
./install-kubectl.sh

prompt_info "Setting up Minio"
if ! kubectl rollout status --timeout=0s deploy/minio; then
	arkade install minio
	# TODO: figure out persistant storage
	#    --set persistence.enabled=true
	kubectl rollout status deploy/minio
fi
    
# Forward the minio port to your machine
pgrep -f kubectl.*9000 && (pkill -f kubectl.*9000 > /dev/null; sleep 3)
kubectl port-forward svc/minio 9000:9000 &
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
