#!/bin/bash

cd $(dirname $0)

. color-prompt.sh

# Install dependencies
./install-arkade_helm_faas-cli_mc.sh || exit $?
./install-kubectl.sh || exit $?

# Timeout in minutes
TIMEOUT=10m
REPLICAS=2

[ -z "$1" ] && ELB=true || ELB=false

# The AWS Elastic Load Balancer defaults to a 1 minute timeout, increase timeout
set_elb_idle_timeout() {
	# convert timeout to seconds
	local timeout=$((${TIMEOUT%?}*60))
	local i
	for ((i=0; i<180; i++)); do
		export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080
		if [ "$OPENFAAS_URL" = ":8080" ]; then
			sleep 1
			continue
		fi
		sleep 3
		break
	done
	if [ "$OPENFAAS_URL" = ":8080" ]; then
		prompt_error "Failed to detect external gateway endpoint, try again!"
		return 1
	fi

	# If current timeout matches what we will set, then return
	local load_balancer_name=$(echo $OPENFAAS_URL | sed 's/-.*//')
	local current_timeout=$(aws elb describe-load-balancer-attributes \
		--load-balancer-name $load_balancer_name  | awk '/IdleTimeout/ {print $NF; exit}')
	[ -n "$current_timeout" ] && [ $current_timeout = $timeout ] && return 0

	# Modify load balancer idle timeout for longer running functions
	echo "Wait while OpenFaaS gateway is configured"
	local json="{\"ConnectionSettings\":{\"IdleTimeout\":$timeout}}"
	aws elb modify-load-balancer-attributes \
		--load-balancer-name $load_balancer_name \
		--load-balancer-attributes "$json"
	local ret=$?
	[ $ret -eq 0 ] && sleep 1
	return $ret
}

# Check if openfaas is installed in cluster already
if ! kubectl rollout status --timeout=0s -n openfaas deploy/gateway 2> /dev/null; then
	prompt_info "Setting up OpenFaas"

	# If AWS Elastic Load Balancer is used we need to check for a valid AWS account
	if $ELB && ! aws sts get-caller-identity > /dev/null; then
		prompt_error "aws cli is not configured correctly or AWS_PROFILE is not set."
		echo "Export AWS_PROFILE user (e.g. export AWS_PROFILE=kops) and try again."
		exit 1
	fi

	# Install openfaas with arkade, configure for long running functions
	cmd="arkade install openfaas
		--set gateway.upstreamTimeout=$TIMEOUT
		--set gateway.writeTimeout=$TIMEOUT
		--set gateway.readTimeout=$TIMEOUT
		--set faasnetes.writeTimeout=$TIMEOUT
		--set faasnetes.readTimeout=$TIMEOUT
		--set queueWorker.ackWait=$TIMEOUT
		--set gateway.replicas=$REPLICAS
		--set queueWorker.replicas=$REPLICAS
		--set queueWorker.maxInflight=100"
	# Check if external load balancer is used
	$ELB && cmd+=" --set serviceType=LoadBalancer"
	eval $cmd

	# check that openfaas is deployed
	kubectl rollout status -n openfaas deploy/gateway
fi

# Configure external load balancer timeout
if $ELB; then
	set_elb_idle_timeout
	exit $?
fi
