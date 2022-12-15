#!/bin/bash

set -o pipefail

main() {
	cd $(dirname $0)

	# Add colorized console prompts
	. ../utility/color-prompt.sh

	# Deploy OpenFaaS
	if ! ../utility/setup-openfaas.sh; then
		prompt_error "OpenFaaS failed to deploy"
		exit 1
	fi

	# Export OpenFaaS gateway
	export_gateway_url

	# Check that function name is passed
	if [ -z "$1" ]; then
		prompt_error "Function name must be passed as the first positional argument"
		exit 1
	fi
	FN_NAME=$1
	shift

	# Default arguments
	ITERATION=1
	REPLICAS=1
	CONCURRENT=false

	# Parse parameters
	while [ -n "$1" ]; do
		case $1 in
			-i|--iterations)
				ITERATION=$2
				shift
				shift
			;;
			-t|--type)
				CLUSTER_TYPE="$2"
				shift
				shift
			;;
			-c|--concurrent)
				CONCURRENT=true
				shift
			;;
			-d|--delete)
				remove_fn
				exit $?
			;;
			-r|--replicas)
				REPLICAS=$2
				echo "Replicas set to $REPLICAS"
				shift
				shift
			;;
			-h|--help)
				usage 0
			;;
			*)
				prompt_error "Invalid argument passed '$1'"
				usage 1
			;;
		esac
	done
}

usage() {
	cat <<-_OPTIONS
	OPTIONS:
	  -c | --concurrent		Execute functions concurrently (defaults to false)
	  -d | --delete			Remove function
	  -h | --help			Print this usage message then exit
	  -i <int> | --iterations <int>	Number of function iterations (defaults to 1)
	  -r <int> | --replicas <int>	Number of function replicas to spawn (defaults to 1)
	  -t <str> | --type <str>	Cluster type
	_OPTIONS
	exit $1
}

export_gateway_url() {
	export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080
	[ "$OPENFAAS_URL" = ":8080" ] && export OPENFAAS_URL="localhost:8080"
	echo "OpenFaas using gateway: $OPENFAAS_URL"
}

# Login faas-cli
faas_login() {
	kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | \
	base64 --decode | \
	faas-cli login --username admin --password-stdin || return $?
}

# Build, push, and deploy OpenFaaS function
deploy_fn() {
	../utility/install-docker.sh || return $?

	faas_login || return $?

	local yaml="$FN_NAME.yml"
	local image=$(awk '/image:/ {print $NF}' $yaml)
	if [ -z "$(docker images -q $image)" ]; then
		faas-cli up -f $yaml \
			--label com.openfaas.scale.min=$REPLICAS \
			--label com.openfaas.scale.max=$REPLICAS || return $?
	else
		read -n 1 -p "docker image is already present locally, enter 'y' to rebuild. " ANS
		echo
		# faas-cli doesn't seem to update or replace existing functions correctly, we need to remove first
		kubectl get service/$FN_NAME -n openfaas-fn &> /dev/null && \
		faas-cli remove -f $yaml
		if [ "${ANS,}" = 'y' ]; then
			faas-cli up -f $yaml \
				--label com.openfaas.scale.min=$REPLICAS \
				--label com.openfaas.scale.max=$REPLICAS || return $?
		else
			faas-cli push -f $yaml && \
			faas-cli deploy -f $yaml \
				--label com.openfaas.scale.min=$REPLICAS \
				--label com.openfaas.scale.max=$REPLICAS || return $?
		fi
	fi
	echo "Checking deployment status, this might take awhile for some functions"
	local i
	for ((i=0; i<100; i++)); do
		if [ $(kubectl get deployment.app/$FN_NAME -n openfaas-fn \
				--output="jsonpath={.status.conditions[0].status}") = 'True' ]; then
			return 0
		fi
		sleep 1
	done
	prompt_error "Failed checking deployment status of function"
	return 1
}

# Remove OpenFaaS function
remove_fn() {
	faas_login || return $?

	faas-cli remove -f $FN_NAME.yml
	return $?
}

# Execute OpenFaaS
execute_fn() {
	if [ -z "$CLUSTER_TYPE" ]; then
		prompt_error "CLUSTER_TYPE must be passed (e.g. kops)"
		usage 1
	fi

	local payload="$1"
	local fn_discription="$2"
	local datetime=$(date +"%Y-%m-%d_%H-%M-%S")
	local log_dir="../logs/openfaas/$CLUSTER_TYPE/$FN_NAME"
	local log="$log_dir/${datetime}_${FN_NAME}_${fn_discription}.log"

	# Create a directory to store logs
	mkdir -p "$log_dir" || exit $?

	if $CONCURRENT; then
		curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
			http://$OPENFAAS_URL/function/$FN_NAME -o "$log" &
		PROCESSES+=($!)
		LOGS+=("$log")
	elif curl -s -H "Content-Type: application/json" -X POST -d "$payload" \
			http://$OPENFAAS_URL/function/$FN_NAME -o "$log"; then
		local status=$(jq -r '.status' $log)
		if [ -z "$status" ] || [ $status -ne 200 ]; then
			prompt_error "function exit status is $status"
		fi
		printf "\nFunction completed\nSee $log\n"
		sleep 2
	else
		prompt_error "Failed to execute $FN_NAME $fn_discription"
		exit 1
	fi
}

# Check concurrent functions
check_concurrent_fn() {
	$CONCURRENT || return
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

main $@
