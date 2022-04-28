#!/bin/bash
cd $(dirname $0)

# Color prompt for easier info output
. color-prompt.sh

check_exit_status() {
	if ! ./$1; then
		prompt_error $2
		exit 1
	fi
}

# Check for required binaries
prompt_info "Checking for installed binary dependencies"
if which arkade > /dev/null || curl -sLS https://get.arkade.dev | sudo sh; then
	for bin in kind kubectl helm faas-cli mc; do
		which $bin > /dev/null || \
		(arkade get $bin && sudo install -m 0755 $HOME/.arkade/bin/$bin /usr/local/bin)
	done
	# Create a simlink for faas-cli
	which faas > /dev/null || sudo ln -fs /usr/local/bin/faas-cli /usr/local/bin/faas
else
	prompt_error "Did not find arkade, debug me!"
	exit 1
fi

# Deploy kind (k8s in docker) along with local container registry
prompt_info "Setting up kind cluster with local container registry"
check_exit_status setup-kind-with-registry.sh \
	"Encountered an error setting up local cluster and registry"

# Install openfaas with arkade, configure for long running functions
prompt_info "Installing OpenFaaS into Kubernetes cluster"
check_exit_status setup-openfaas.sh \
	"Encountered an error setting up OpenFaaS"

# Install minio
prompt_info "Setting up Minio"
check_exit_status setup-minio.sh \
	"Encountered an error setting up Minio"

# Setup nlp pipeline buckets and transfer data
prompt_info "Deploy nlp pipeline buckets, data, function"
check_exit_status deploy-nlp.sh \
	"Encountered an error deploying nlp pipeline"
