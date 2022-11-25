#!/bin/bash

set -o pipefail

# Add colorized console prompts
. $(dirname $0)/../utility/color-prompt.sh

export_gateway_url() {
	export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080 \
	&& echo "OpenFaaS gateway: $OPENFAAS_URL" || echo "OpenFaaS gateway not found"
	[ "$OPENFAAS_URL" = ":8080" ] && export OPENFAAS_URL="localhost:8080"
}

# Build, push, and deploy OpenFaaS function
deploy_fn() {
	../utility/setup-openfaas.sh
	../utility/setup-minio.sh
	../utility/install-docker.sh
	../utility/pass-minio-secrets.sh || return $?

	export_gateway_url

	# login faas-cli
	kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | \
	base64 --decode | \
	faas-cli login --username admin --password-stdin || return $?

	local yaml="$1"
	if [ -z "$(docker images -q $(awk '/image:/ {print $NF}' $yaml))" ]; then
		faas-cli up -f $yaml
	else
		read -n 1 -p "docker image is already present locally, enter 'y' to rebuild. " ANS
		if [ "${ANS,}" = 'y' ]; then
			faas-cli up -f $yaml --replace=true --update=false
		elif [ "$(docker inspect -f "{{.RepoDigests}}" $(awk '/image:/ {print $NF}' $yaml))" = "[]" ]; then
			faas-cli push -f $yaml && faas-cli deploy -f $yaml
		else
			faas-cli deploy -f $yaml --replace=true --update=false
		fi
	fi
	return $?
}

# Setup for OpenFaaS function invocation
invoke_setup() {
	# Verify all arguments are set
	if [ -z "$3" ]; then
		cat <<-_USAGE
		Argumentes include
		\$1 - iterations (set to '$1')
		\$2 - cluster type (set to '$2')
		\$3 - function name (set to '$3')
		_USAGE
		exit 1
	fi

	# Verify arguments are all set
	ITERATION=$1
	CLUSTER_TYPE="$2"
	FUNCTION_NAME="$3"

	# Check for OpenFaaS endpoint associated with load balancer URL
	export_gateway_url

	# Create a directory to store logs
	LOG_DIR="logs/$CLUSTER_TYPE/$FUNCTION_NAME"
	[ -d "$LOG_DIR" ] || mkdir -p "$LOG_DIR"
	return $?
}

# Execute OpenFaaS
execute_fn() {
	local fn_name="$1"
	local payload="$2"
	local fn_discription="$3"
	local datetime=$(date +"%Y-%m-%d_%H-%M-%S")

	if curl -s -H "Content-Type: application/json" -X POST -d "$payload" http://$OPENFAAS_URL/function/$fn_name | \
		tee "$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"; then
		sleep 2
	else
		prompt_error "Failed to execute $fn_name $fn_discription"
		exit 1
	fi
}
