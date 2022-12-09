#!/bin/bash

set -o pipefail

# Add colorized console prompts
. $(dirname $0)/../utility/color-prompt.sh

export_gateway_url() {
	export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080
	if [ "$OPENFAAS_URL" = ":8080" ]; then
		export OPENFAAS_URL="localhost:8080"
		echo "OpenFaaS load balancer endpoint not found"
	fi
	echo "OpenFaas using gateway: $OPENFAAS_URL"
}

deploy_minio() {
	../utility/setup-minio.sh
	return $?
}

# Build, push, and deploy OpenFaaS function
deploy_fn() {
	../utility/setup-openfaas.sh || return $?
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
		faas-cli up -f $yaml || return $?
	else
		read -n 1 -p "docker image is already present locally, enter 'y' to rebuild. " ANS
		echo
		# faas-cli doesn't seem to update or replace existing functions correctly, we need to remove first
		kubectl get service/${yaml%.*} -n openfaas-fn &> /dev/null && \
		faas-cli remove -f $yaml
		if [ "${ANS,}" = 'y' ]; then
			faas-cli up -f $yaml || return $?
		else
			faas-cli push -f $yaml && faas-cli deploy -f $yaml || return 1
		fi
	fi
	echo "Checking deployment status, this might take awhile for some functions"
	local i
	for ((i=0; i<100; i++)); do
		if [ $(kubectl get deployment.app/${yaml%.*} -n openfaas-fn \
				--output="jsonpath={.status.conditions[0].status}") = 'True' ]; then
			return 0
		fi
		sleep 1
	done
	return 1
}

# Remove OpenFaaS function
remove_fn() {
	../utility/setup-openfaas.sh || return $?

	export_gateway_url

	# login faas-cli
	kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | \
	base64 --decode | \
	faas-cli login --username admin --password-stdin || return $?

	local yaml="$1"
	faas-cli remove -f $yaml
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
	local log="$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"
	if [ -n "$CONCURRENT" ]; then
		curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
			http://$OPENFAAS_URL/function/$fn_name -o "$log" &
		PROCESSES+=($!)
		LOGS+=("$log")
	elif curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
			http://$OPENFAAS_URL/function/$fn_name -o "$log"; then
		local status=$(jq -r '.status' $log)
		[ $status -ne 200 ] && prompt_error "function exit status is $status"
		printf "\nSuccessfully executed function\nSee $log\n"
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
		wait ${PROCESSES[$i]}
		local ret=$?
		if [ $ret -eq 0 ]; then
			printf "\nCompleted process ${PROCESSES[$i]}\nSee ${LOGS[$i]}\n"
			local status=$(jq -r '.status' ${LOGS[$i]})
			[ $status -ne 200 ] && prompt_error "function exit status is $status"
		else
			prompt_error "Encountered an error status ($ret) with process ${PROCESSES[$i]}\nSee ${LOGS[$i]}"
		fi
	done
	unset PROCESSES
	unset LOGS
}
