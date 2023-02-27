#!/bin/bash

minio_port_forward() {
    # Forward the minio port to your machine
	pgrep -f kubectl.*9000 > /dev/null && (pkill -f kubectl.*9000; sleep 3)
    kubectl port-forward svc/minio 9000:9000 &
	if [ $? -ne 0 ]; then
        echo "Failed to port-forward minio service"
        return 1
    fi
    sleep 1

    # Get the access and secret key to gain access to minio
    local accesskey=$(kubectl get secret -n default minio -o jsonpath="{.data.accesskey}" | base64 --decode)
    local secretkey=$(kubectl get secret -n default minio -o jsonpath="{.data.secretkey}" | base64 --decode)

    # Add a host
    mc config host add minio http://localhost:9000 $accesskey $secretkey
    local ret=$?
    [ $ret -eq 0 ] && sleep 1

    return $ret
}

cd $(dirname $0)

. color-prompt.sh

# Install dependencies
./install-arkade_helm_faas-cli_mc.sh
./install-kubectl.sh

if ! kubectl rollout status --timeout=0s deploy/minio 2> /dev/null; then
	prompt_info "Setting up Minio"
	arkade install minio || exit $?
	# TODO: figure out persistant storage
	#    --set persistence.enabled=true
	kubectl rollout status deploy/minio || exit $?
fi

minio_port_forward
exit $?
