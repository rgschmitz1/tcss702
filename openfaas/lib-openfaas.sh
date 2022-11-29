#!/bin/bash

set -o pipefail

# Add colorized console prompts
. $(dirname $0)/../utility/color-prompt.sh

export_gateway_url() {
	export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080
	if [ "$OPENFAAS_URL" = ":8080" ]; then
		export OPENFAAS_URL="localhost:8080"
		echo "OpenFaaS load balancer gateway not found"
	else
		echo "OpenFaaS load balancer gateway: $OPENFAAS_URL"
	fi
}

# Build, push, and deploy OpenFaaS function
deploy_fn() {
	../utility/setup-openfaas.sh || return $?
	../utility/setup-minio.sh || return $?
	../utility/install-docker.sh || return $?
	../utility/pass-minio-secrets.sh || return $?

	export_gateway_url

	# login faas-cli
	kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | \
	base64 --decode | \
	faas-cli login --username admin --password-stdin || return $?

	local yaml="$1"
	local image=$(awk '/image:/ {print $NF}' $yaml)
	if [ -z "$(docker images -q $image)" ]; then
		faas-cli up -f $yaml
	else
		read -n 1 -p "docker image is already present locally, enter 'y' to rebuild. " ANS
		echo
		# faas-cli doesn't seem to update or replace existing functions correctly, we need to remove first
		kubectl get service/${yaml%.*} -n openfaas-fn &> /dev/null && \
		faas-cli remove -f $yaml
		if [ "${ANS,}" = 'y' ]; then
			faas-cli up -f $yaml
		else
			faas-cli push -f $yaml && faas-cli deploy -f $yaml
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
	if [[ -n "$4" && "$4" = '-c' ]]; then
		CONCURRENT=true
	fi


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
	if [ -n "$4" ]; then
		local log="$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"
		curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
		http://$OPENFAAS_URL/function/$fn_name \
		-o "$log" &
		PROCESSES+=($!)
		LOGS+=("$log")
	elif curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
		http://$OPENFAAS_URL/function/$fn_name \
		-o "$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"; then
		cat "$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"
		sleep 2
	else
		prompt_error "Failed to execute $fn_name $fn_discription"
		exit 1
	fi
}

# Check concurrent functions
check_concurrent_fn() {
	local i
	for ((i=0; i<${#PROCESSES[@]}; i++)); do
		if wait ${PROCESSES[$i]}; then
			printf "\nSuccessfully completed process ${PROCESSES[$i]}\nSee ${LOGS[$i]}\n"
		else
			prompt_error "Encountered an error with process ${PROCESSES[$i]}\nSee ${LOGS[$i]}"
		fi
	done
}
