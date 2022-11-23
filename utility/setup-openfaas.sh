#!/bin/bash

cd $(dirname $0)

# Install dependencies
./install-arkade_helm_faas-cli_mc.sh
./install-kubectl.sh

# The AWS Elastic Load Balancer defaults to a 1 minute timeout, increase to 10 minutes
set_elb_idle_timeout() {
	export OPENFAAS_URL=$(kubectl get svc -n openfaas gateway-external -o jsonpath='{.status.loadBalancer.ingress[*].hostname}'):8080 \
	&& echo "Your gateway URL is: $OPENFAAS_URL"
	if [ "$OPENFAAS_URL" = ":8080" ]; then
		echo "Failed to detect external gateway endpoint, try again!"
		exit 1
	fi
	aws elb modify-load-balancer-attributes \
		--load-balancer-name $(echo $OPENFAAS_URL|sed 's/-.*//') \
		--load-balancer-attributes "{\"ConnectionSettings\":{\"IdleTimeout\":600}}"
}

TIMEOUT=10m
REPLICAS=2

if ! kubectl rollout status --timeout=0s -n openfaas deploy/gateway; then
	# Install openfaas with arkade, configure for long running functions
	cmd="arkade install openfaas
		--set gateway.upstreamTimeout=$TIMEOUT
		--set gateway.writeTimeout=$TIMEOUT
		--set gateway.readTimeout=$TIMEOUT
		--set faasnetes.writeTimeout=$TIMEOUT
		--set faasnetes.readTimeout=$TIMEOUT
		--set queueWorker.ackWait=$TIMEOUT
		--set gateway.replicas=$REPLICAS
		--set queueWorker.replicas=$REPLICAS"
	# Check if external load balancer is used
	if [ -n "$1" ] && [ "$1" = "-l" ]; then
		cmd+=" --set serviceType=LoadBalancer
			--set operator.create=true"
	fi
	eval $cmd

	# check that openfaas is deployed
	kubectl rollout status -n openfaas deploy/gateway
fi

if [ -n "$1" ] && [ "$1" = "-l" ]; then
	set_elb_idle_timeout
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
PASSWORD=$(kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | base64 --decode)
printf $PASSWORD | faas-cli login --username admin --password-stdin
exit $?
