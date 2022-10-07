#!/bin/bash

# This IAM user should be created in AWS with the following permissions:
#   AmazonEC2FullAccess
#   AmazonRoute53FullAccess
#   AmazonS3FullAccess
#   IAMFullAccess
#   AmazonVPCFullAccess
#   AmazonSQSFullAccess
#   AmazonEventBridgeFullAccess
export AWS_PROFILE=tcss702

# you should see a list of all your IAM users here
aws iam list-users

# Because "aws configure" doesn't export these vars for kops to use, we export them now
export AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
export AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)

# Setup AWS Route53 to access cluster via DNS subdomain,
#   ID=$(uuidgen) && aws route53 create-hosted-zone --name tcss702.rgschmitz.com --caller-reference $ID | jq .DelegationSet.NameServers
#
# You will now go to your registrar's page and log in.
# You will need to create a new SUBDOMAIN, and use the 4 NS records received from the above command for the new SUBDOMAIN.
# This MUST be done in order to use your cluster.
# Do NOT change your top level NS record, or you might take your site offline.

# An S3 bucket must be created to store the state of the k8s cluster
#   aws s3 mb s3://tcss702-rgschmitz-com-state-store
#   aws s3api put-bucket-versioning \
#   	--bucket tcss702-rgschmitz-com-state-store \
#   	--versioning-configuration Status=Enabled
#   aws s3api put-bucket-encryption \
#   	--bucket tcss702-rgschmitz-com-state-store \
#   	--server-side-encryption-configuration \
#   	'{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'

# Check if kops configuration file exists, then export the kops state store variable
if [ ! -f ~/.kops.yaml ]; then
	echo '---' > ~/.kops.yaml
fi
if ! grep -q kops_state_store ~/.kops.yaml; then
	echo 'kops_state_store: s3://tcss702-rgschmitz-com-state-store' >> ~/.kops.yaml
fi

export NAME=tcss702.rgschmitz.com

# Create k8s cluster configuration
kops create cluster --name ${NAME} --zones us-east-2a --master-size t3.medium --node-size t3.large --cloud aws --networking calico --ssh-public-key ~/.ssh/id_rsa.pub
# Deploy k8s cluster
kops update cluster --name ${NAME} --yes
# Wait until the cluster is up and ready to use
kops validate cluster --wait 10m


# ---
# Use the following to delete the k8s cluster
#   kops delete cluster --name ${NAME} --yes
