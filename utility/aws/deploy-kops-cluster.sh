#!/bin/bash

# kops configuration variables
NAME=tcss702.rgschmitz.com  # Domain was purchased from namecheap.com
CLOUD=aws
REGION=us-east-2
MASTER_SIZE=t3.medium
NODE_SIZE=t3.large
SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
NETWORK_CNI=calico
TIMEOUT=20m
export KOPS_STATE_STORE=s3://tcss702-rgschmitz-com-state-store


# Check for and install dependencies
install_dependencies() {
	if ! which jq > /dev/null; then
		sudo apt update
		sudo apt install -y jq || return $?
	fi
	if ! which kops > /dev/null; then
		curl -LO https://github.com/kubernetes/kops/releases/download/v1.25.2/kops-linux-amd64 && \
		sudo install -o root -g root -m 0755 kops-linux-amd64 /usr/local/bin/kops && \
		rm kops-linux-amd64 || return 1
	fi
	return 0
}


# An IAM user must be created in AWS with permissions to create a cluster
create_kops_iam_user() {
	# IAM policies
	local policy=(\
		AmazonEC2FullAccess \
		AmazonRoute53FullAccess \
		AmazonS3FullAccess \
		IAMFullAccess \
		AmazonVPCFullAccess \
		AmazonSQSFullAccess \
		AmazonEventBridgeFullAccess \
	)

	# Profile name will be grepped from this file incase it changes later
	local profile=$(awk -F'=' '/export AWS_PROFILE/ {print $NF}' $0)

	aws iam create-group --group-name $profile || return $?

	# Attach policies to group
	for p in ${policy[@]}; do
		aws iam attach-group-policy \
			--policy-arn arn:aws:iam::aws:policy/$p \
			--group-name $profile || return $?
	done

	aws iam create-user --user-name $profile || return $?

	aws iam add-user-to-group --user-name $profile --group-name $profile || \
		return $?

	# Create Access Key ID and Secret Access Key
	local tmp=$(mktemp)
	if ! aws iam create-access-key --user-name $profile > $tmp; then
		local ret=$?
		rm $tmp
		return $ret
	fi
	local access_key=$(jq -r .AccessKey.AccessKeyId < $tmp)
	local secret_key=$(jq -r .AccessKey.SecretAccessKey < $tmp)

	# Populate AWS credentials
	cat <<- _CREDENTIALS >> $HOME/.aws/credentials

	[$profile]
	aws_access_key_id=$access_key
	aws_secret_access_key=$secret_key
	_CREDENTIALS

	# Populate AWS config
	cat <<- _CONFIG >> $HOME/.aws/config

	[profile $profile]
	region=$REGION
	_CONFIG

	# Remove temporary AWS config file
	rm $tmp

	# Verify kops user was successfully created
	aws iam get-user --user $profile
	return $?
}


# An S3 bucket must be created to store the state of the k8s cluster
create_kops_state_store() {
	local state_store=$(basename $KOPS_STATE_STORE)
	aws s3 mb s3://$state_store || return $?
	aws s3api put-bucket-versioning \
		--bucket $state_store \
		--versioning-configuration Status=Enabled || return $?
	aws s3api put-bucket-encryption \
		--bucket $state_store \
		--server-side-encryption-configuration \
		'{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'
	return $?
}


# Setup AWS Route53 to access cluster via DNS subdomain,
create_route53_dns() {
	local id=$(uuidgen) && \
		aws route53 create-hosted-zone --name $NAME --caller-reference $id | \
		jq .DelegationSet.NameServers
	local ret=$?
	cat <<-_SUBDOMAIN
	You will need to create a new SUBDOMAIN, and use the 4 NS records received from the above command for the new SUBDOMAIN.
	This MUST be done in order to use your cluster.
	Do NOT change your top level NS record, or you might take your site offline.
	_SUBDOMAIN
	return $ret
}


# Create kubernetes cluster
create_cluster() {
	# Create cluster configuration
	kops create cluster \
		--name $NAME \
		--cloud $CLOUD \
		--zones ${REGION}a \
		--master-size $MASTER_SIZE \
		--node-size $NODE_SIZE \
		--networking $NETWORK_CNI \
		--ssh-public-key $SSH_PUBLIC_KEY
	# Deploy cluster
	kops update cluster --name $NAME --yes --admin
}


# Delete kubernetes cluster
delete_cluster() {
	kops delete cluster --name $NAME --yes
	return $?
}


# Check for script dependencies
if ! install_dependencies; then
	printf "\nERROR: failed to install script dependencies\n"
	exit 1
fi

# Verify kops user has been created
if ! (grep -q kops $HOME/.aws/credentials || create_kops_iam_credentials); then
	printf "\nERROR: failed to generate kops user\n"
	exit 1
else
	export AWS_PROFILE=kops

	# awscli do not export these variables for kops to use
	export AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
	export AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)
fi

# Check for user input
while [ -n "$1" ]; do
	case $1 in
		-d)
			if delete_cluster; then
				exit 0
			else
				printf "\nERROR: Encountered an issue deleting cluster\n"
				exit 1
			fi
		;;
		*)
			printf "\nERROR: invalid argument!\n"
			exit 1
		;;
	esac
done

# Verify kops state store has been created
if ! (aws s3 ls | grep -q $(basename $KOPS_STATE_STORE) || create_kops_state_store); then
	printf "\nERROR: failed to create kops state store in S3\n"
	exit 1
fi

# Route53 DNS needs to be deployed to use kops
if ! (aws route53 list-hosted-zones | grep -q $NAME || create_route53_dns); then
	printf "\nERROR: failed to create Route53\n"
	exit 1
fi

# Wait until the cluster is up and ready to use
if create_cluster && kops validate cluster --wait $TIMEOUT; then
	runtime=$(date -ud "@$SECONDS" "+%M minutes, %S seconds")
	printf "\nSuccessfully deployed cluster in $runtime\n"
else
	printf "\nERROR: Failed to deploy cluster\n"
	delete_cluster
	exit 1
fi
