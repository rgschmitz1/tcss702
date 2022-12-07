#!/bin/bash
cd $(dirname $0)

# Color prompt for easier info output
. color-prompt.sh

check_exit_status() {
	./$1 && return
	prompt_error $2
	exit 1
}

# Install required dependencies
check_exit_status install-arkade_helm_faas-cli_mc.sh \
	"Encountered an error installing dependencies"
check_exit_status install-docker.sh \
	"Encountered an error installing docker"
check_exit_status install-kind.sh \
	"Encountered an error installing kind"

if [ "$(kind get clusters)" = 'kind' ]; then
	kind export kubeconfig || exit $?
else
	prompt_info "Starting kind cluster (k8s in docker)"
	kind create cluster || exit $?
fi
#prompt_info "Starting kind cluster (k8s in docker) with local container registry"
#check_exit_status setup-kind-with-registry.sh \
#	"Encountered an error setting up local cluster and registry"

# Install openfaas, configure for long running functions
check_exit_status "setup-openfaas.sh -l" \
	"Encountered an error setting up OpenFaaS"

# Install minio
check_exit_status setup-minio.sh \
	"Encountered an error setting up Minio"
