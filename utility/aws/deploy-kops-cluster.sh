#!/bin/bash

# Create k8s cluster configuration
export NAME=tcss702.rgschmitz.com
export AWS_REGION=us-east-2a
export MASTER_SIZE=t3.medium
export NODE_SIZE=t3.large
export CLOUD=aws
export SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
export NETWORK_CNI=calico
export KOPS_STATE_STORE=tcss702-rgschmitz-com-state-store

# This IAM user should be created in AWS with the following permissions:
#   AmazonEC2FullAccess
#   AmazonRoute53FullAccess
#   AmazonS3FullAccess
#   IAMFullAccess
#   AmazonVPCFullAccess
#   AmazonSQSFullAccess
#   AmazonEventBridgeFullAccess
export AWS_PROFILE=tcss702

# you should see a list of all your IAM users using:
#   aws iam list-users

# Because "aws configure" doesn't export these vars for kops to use, we export them now
export AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
export AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)

# An S3 bucket must be created to store the state of the k8s cluster
#   aws s3 mb s3://$KOPS_STATE_STORE
#   aws s3api put-bucket-versioning \
#   	--bucket $KOPS_STATE_STORE \
#   	--versioning-configuration Status=Enabled
#   aws s3api put-bucket-encryption \
#   	--bucket $KOPS_STATE_STORE \
#   	--server-side-encryption-configuration \
#   	'{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'

# Setup AWS Route53 to access cluster via DNS subdomain,
#   ID=$(uuidgen) && aws route53 create-hosted-zone --name $NAME --caller-reference $ID | jq .DelegationSet.NameServers
#
# You will now go to your registrar's page and log in.
# You will need to create a new SUBDOMAIN, and use the 4 NS records received from the above command for the new SUBDOMAIN.
# This MUST be done in order to use your cluster.
# Do NOT change your top level NS record, or you might take your site offline.

# Check if kops configuration file exists, then export the kops state store variable
if [ ! -f $HOME/.kops.yaml ]; then
	echo '---' > $HOME/.kops.yaml
fi
if ! grep -q kops_state_store $HOME/.kops.yaml; then
	echo "kops_state_store: s3://$KOPS_STATE_STORE" >> $HOME/.kops.yaml
fi

# Create cluster configuration
kops create cluster \
	--name $NAME \
	--cloud $CLOUD \
	--zones $AWS_REGION \
	--master-size $MASTER_SIZE \
	--node-size $NODE_SIZE \
	--networking $NETWORK_CNI \
	--ssh-public-key $SSH_PUBLIC_KEY
# Deploy cluster
kops update cluster --name $NAME --yes --admin
# Wait until the cluster is up and ready to use
kops validate cluster --wait 20m

printf "\nThis script took $SECONDS seconds to finish\n"

# ---
# Use the following to delete the k8s cluster
#   kops delete cluster --name $NAME --yes
