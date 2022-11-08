#!/bin/bash

# AWS user profile name
PROFILE=kops

# kops configuration variables
NAME=tcss702.rgschmitz.com
CLOUD=aws
REGION=us-east-2
ZONES=us-east-2c
MASTER_SIZE=m5.large
MASTER_COUNT=1
NODE_SIZE=c5n.2xlarge
NODE_COUNT=1
SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
NETWORK_CNI=calico
K8S_VERSION=1.25.3
export KOPS_STATE_STORE=s3://tcss702-rgschmitz-com-state-store

# The cluster will be deleted if this timeout is exceeded during validation
TIMEOUT=20m


# Display script usage/flags for user
usage() {
	cat <<- _USAGE
	Usage: $(basename $0) [Options]

	Description: This script will instantiate a kubernetes cluster on AWS

	Options:
	  -d|--delete, deletes a running cluster and all underlying infrastructure
	  -e|--edit, edit cluster configuration before deployment
	  -f|--filename <cluster-spec.yml>, pass a custom cluster spec
	  -h|--help, print this message
	_USAGE
}


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
	$(dirname $0)/install-awscli.sh
	return $?
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

	aws iam create-group --group-name $PROFILE || return $?

	# Attach policies to group
	for p in ${policy[@]}; do
		aws iam attach-group-policy \
			--policy-arn arn:aws:iam::aws:policy/$p \
			--group-name $PROFILE || return $?
	done

	aws iam create-user --user-name $PROFILE || return $?

	aws iam add-user-to-group --user-name $PROFILE --group-name $PROFILE || \
		return $?

	# Create Access Key ID and Secret Access Key
	local tmp=$(mktemp)
	if ! aws iam create-access-key --user-name $PROFILE > $tmp; then
		local ret=$?
		rm $tmp
		return $ret
	fi
	local access_key=$(jq -r .AccessKey.AccessKeyId < $tmp)
	local secret_key=$(jq -r .AccessKey.SecretAccessKey < $tmp)

	# Populate AWS credentials
	cat <<- _CREDENTIALS >> $HOME/.aws/credentials

	[$PROFILE]
	aws_access_key_id=$access_key
	aws_secret_access_key=$secret_key
	_CREDENTIALS

	# Populate AWS config
	cat <<- _CONFIG >> $HOME/.aws/config

	[profile $PROFILE]
	region=$REGION
	_CONFIG

	# Remove temporary AWS config file
	rm $tmp

	# Verify kops user was successfully created
	aws iam get-user --user $PROFILE
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
	if [ -z "$CLUSTER_SPEC" ]; then
		cmd="kops create cluster
			--name $NAME
			--cloud $CLOUD
			--zones $ZONES
			--master-size $MASTER_SIZE
			--master-count $MASTER_COUNT
			--node-size $NODE_SIZE
			--node-count $NODE_COUNT
			--networking $NETWORK_CNI
			--ssh-public-key $SSH_PUBLIC_KEY
			--kubernetes-version=$K8S_VERSION"
		if [ -n "$_EDIT" ]; then
			CLUSTER_SPEC=$(mktemp --suff=.yml)
			cmd+=" --dry-run -o yaml > $CLUSTER_SPEC"
		fi
		eval $cmd
	fi
	# Edit cluster spec before deploying
	[ -n "$_EDIT" ] && vim $CLUSTER_SPEC
	# Create a ssh publickey if deployed from cluster spec
	if [ -n "$CLUSTER_SPEC" ]; then
		kops create -f $CLUSTER_SPEC
		kops create sshpublickey $NAME -i $SSH_PUBLIC_KEY
	fi
	# Deploy cluster
	local ret
	local start=$(date +%s)
	kops update cluster --name $NAME --yes --admin=8760h
	ret=$?
	local end=$(date +%s)
	_RUNTIME=$(date -ud "@$((end-start))" "+%M minutes, %S seconds")
	return $ret
}


# Delete kubernetes cluster
delete_cluster() {
	kops delete cluster --name $NAME --yes
	# Delete ssh key record
	ssh-keygen -f "$HOME/.ssh/known_hosts" -R "api.$NAME"
	return $?
}


# Check for script dependencies
if ! install_dependencies; then
	printf "\nERROR: failed to install script dependencies\n"
	exit 1
fi

# Verify kops user has been created
if ! (grep -q $PROFILE $HOME/.aws/credentials || create_kops_iam_user); then
	printf "\nERROR: failed to generate $PROFILE user\n"
	exit 1
else
	# awscli does not export these variables for kops to use
	export AWS_PROFILE=$PROFILE
	export AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
	export AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)
fi

# Check for user input
while [ -n "$1" ]; do
	case $1 in
		-d|--delete)
			if delete_cluster; then
				exit 0
			else
				printf "\nERROR: Encountered an issue deleting cluster\n"
				exit 1
			fi
		;;
		-e|--edit)
			_EDIT=1
			shift
		;;
		-f|--filename)
			if [[ -z "$2" || ! -f "$2" ]]; then
				printf "\nERROR: a cluster specification file must be passed with this option\n"
				usage
				exit 1
			fi
			CLUSTER_SPEC="$2"
			shift
			shift
		;;
		-h|--help)
			usage
			exit
		;;
		*)
			printf "\nERROR: invalid argument!\n"
			usage
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
	printf "\nSuccessfully deployed cluster in $_RUNTIME\n"
	[ -n "$CLUSTER_SPEC" ] && echo "Deployed from cluster spec, \"$CLUSTER_SPEC\""
else
	delete_cluster
	printf "\nERROR: Failed to deploy cluster\n"
	exit 1
fi
