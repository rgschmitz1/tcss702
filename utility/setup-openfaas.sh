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
	echo "Wait while OpenFaaS gateway is configured"
	for ((i=0; i<180; i++)); do
		export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080
		if [ "$OPENFAAS_URL" = ":8080" ]; then
			sleep 1
			continue
		fi
		echo "Your gateway URL is: $OPENFAAS_URL"
		sleep 3
		break
	done
	if [ "$OPENFAAS_URL" = ":8080" ]; then
		echo "Failed to detect external gateway endpoint, try again!"
		exit 1
	fi
	local json="{\"ConnectionSettings\":{\"IdleTimeout\":$timeout}}"
	aws elb modify-load-balancer-attributes \
		--load-balancer-name $(echo $OPENFAAS_URL | sed 's/-.*//') \
		--load-balancer-attributes "$json"
	return $?
}

# Check if openfaas is installed in cluster already
kubectl rollout status --timeout=0s -n openfaas deploy/gateway 2> /dev/null && exit 0

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

# Configure external load balancer timeout
if $ELB; then
	set_elb_idle_timeout || exit $?
else
	# Forward the gateway to your machine
	pgrep -f kubectl.*8080 > /dev/null && (pkill -f kubectl.*8080; sleep 3)
	kubectl port-forward -n openfaas svc/gateway 8080:8080 &
	if [ $? -ne 0 ]; then
		echo "Failed to port-forward openfaas service"
		exit 1
	fi
fi
sleep 1

# If basic auth is enabled, you can now log into your gateway:
kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | base64 --decode | \
faas-cli login --username admin --password-stdin
exit $?
