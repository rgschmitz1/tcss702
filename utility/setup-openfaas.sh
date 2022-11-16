# Install openfaas with arkade, configure for long running functions
TIMEOUT=10m
cmd="arkade install openfaas
    --set gateway.upstreamTimeout=$TIMEOUT
    --set gateway.writeTimeout=$TIMEOUT
    --set gateway.readTimeout=$TIMEOUT
    --set faasnetes.writeTimeout=$TIMEOUT
    --set faasnetes.readTimeout=$TIMEOUT
    --set queueWorker.ackWait=$TIMEOUT"
# Check if external load balancer is used
if [ -n "$1" ] && [ "$1" == "-l" ]; then
	cmd+=" --set serviceType=LoadBalancer
		--set operator.create=true"
fi
eval $cmd

# Forward the gateway to your machine
kubectl rollout status -n openfaas deploy/gateway
pgrep -f kubectl.*8080 > /dev/null && (pkill -f kubectl.*8080; sleep 3)
kubectl port-forward -n openfaas svc/gateway 8080:8080 &
if [ $? -ne 0 ]; then
	echo "Failed to port-forward openfaas service"
	exit 1
fi
sleep 1

# If basic auth is enabled, you can now log into your gateway:
PASSWORD=$(kubectl get secret -n openfaas basic-auth -o jsonpath="{.data.basic-auth-password}" | base64 --decode)
printf $PASSWORD | faas-cli login --username admin --password-stdin
exit $?
