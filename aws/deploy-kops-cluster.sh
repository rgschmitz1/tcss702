#!/bin/bash

set -o pipefail

# AWS user profile name
PROFILE=kops

# kops configuration variables
NAME=tcss702.k8s.local
CLOUD=aws
REGION=us-east-2
ZONES=us-east-2c
MASTER_SIZE=m5.large
MASTER_COUNT=1
NODE_SIZE=c5.2xlarge
NODE_COUNT=1
SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
NETWORK_CNI=calico
K8S_VERSION=1.24.10
export KUBECONFIG=$HOME/.kube/kops/clusters/$NAME.config
export KUBECTL_VERSION='v1.24.10'
export KOPS_VERSION='v1.25.3'
export KOPS_STATE_STORE=s3://tcss702-rgschmitz-com-state-store

# The cluster will be deleted if this timeout is exceeded during validation
TIMEOUT=20m
LOG_DIR='../logs/kops'

cd $(dirname $0)

# Add colorized console prompt
. ../utility/color-prompt.sh || exit 1


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
	  -s|--suffix, appended folder on log directory (e.g. logs/kops/<suffix>)
	_USAGE
}


# Check for and install dependencies
install_dependencies() {
	../utility/install-kubectl.sh || return $?
	../utility/install-kops.sh || return $?
	../utility/install-jq.sh || return $?
	../utility/install-awscli.sh
	return $?
}


# An IAM user must be created in AWS with permissions to create a cluster
create_kops_iam_user() {
	grep -q $PROFILE $HOME/.aws/credentials && return 0

	# IAM policies
	local policy=(\
		AmazonEC2FullAccess \
		AmazonRoute53FullAccess \
		AmazonS3FullAccess \
		IAMFullAccess \
		AmazonVPCFullAccess \
		AmazonSQSFullAccess \
		AmazonEventBridgeFullAccess \
		AmazonSSMReadOnlyAccess \
	)

	# Create user and group
	if ! aws iam get-user --user $PROFILE &> /dev/null; then
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
	fi

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

	# Give some time for credentials to be activated
	# for some reason this does not seem to happen quickly
	sleep 10

	return 0
}


# An S3 bucket must be created to store the state of the k8s cluster
create_kops_state_store() {
	# Check if bucket is already created and if we have permission to access it
	aws s3api head-bucket --bucket $(basename $KOPS_STATE_STORE) && return 0

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
	aws route53 list-hosted-zones | grep -q $NAME && return 0

	local id=$(uuidgen) && \
		aws route53 create-hosted-zone --name $NAME --caller-reference $id | \
			jq .DelegationSet.NameServers
	local ret=$?
	cat <<-_SUBDOMAIN
	You will need to create a new SUBDOMAIN, and use the 4 NS records received from the above command for the new SUBDOMAIN.
	This MUST be done in order to use your cluster.
	Do NOT change your top level NS record, or you might take your site offline.

	_SUBDOMAIN
	read -n 1 -p "Press anykey to proceed."
	return $ret
}


# Create kubernetes cluster
create_cluster() {
	mkdir -p $(dirname $KUBECONFIG)

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
	kops update cluster --name $NAME --yes --admin=8760h | tee $LOG
	return $?
}


# Delete kubernetes cluster
delete_cluster() {
	[ -n "$SUFFIX" ] && LOG_DIR+="/${SUFFIX}"
	if [ -z "$CLUSTER_SPEC" ]; then
		LOG=$LOG_DIR/$(date +"%Y-%m-%d_%H-%M-%S")_kops_shutdown.log
	else
		LOG=$LOG_DIR/$(date +"%Y-%m-%d_%H-%M-%S")_$(echo $CLUSTER_SPEC | sed 's|.*/\(.*\)\.yml|\1|')_shutdown.log
	fi
	local start=$(date +%s)
	kops delete cluster --name $NAME --yes
	local ret=$?
	local end=$(date +%s)
	local runtime=$(date -d "@$((end-start))" "+%M minutes, %S seconds")
	printf "\nShutdown cluster in $runtime\n" | tee -a $LOG

	return $ret
}


# Check for script dependencies
if ! install_dependencies; then
	prompt_error "failed to install script dependencies"
	exit 1
fi

# Verify kops user has been created
if ! create_kops_iam_user; then
	prompt_error "failed to generate $PROFILE user"
	exit 1
fi
# awscli does not export these variables for kops to use
export AWS_PROFILE=$PROFILE
export AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
export AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)

# Check for user input
while [ -n "$1" ]; do
	case $1 in
		-d|--delete)
			DELETE=true
			shift
		;;
		-e|--edit)
			_EDIT=1
			shift
		;;
		-f|--filename)
			if [[ -z "$2" || ! -f "$2" ]]; then
				prompt_error "a cluster specification file must be passed with this option"
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
		-s|--suffix)
			SUFFIX="$2"
			shift
			shift
		;;
		*)
			prompt_error "invalid argument!"
			usage
			exit 1
		;;
	esac
done

if [ -n "$DELETE" ] && $DELETE; then
	if delete_cluster; then
		exit 0
	else
		prompt_error "encountered an issue deleting cluster"
		exit 1
	fi
fi

# Verify kops state store has been created
if ! create_kops_state_store; then
	prompt_error "failed to create kops state store in S3"
	exit 1
fi

# Route53 DNS needs to be deployed to use kops unless using gossip based DNS
#if ! create_route53_dns; then
#	prompt_error "failed to create Route53"
#	exit 1
#fi

# Create log for kops deployment
[ -n "$SUFFIX" ] && LOG_DIR+="/${SUFFIX}"
mkdir -p $LOG_DIR
if [ -z "$CLUSTER_SPEC" ]; then
	LOG=$LOG_DIR/$(date +"%Y-%m-%d_%H-%M-%S")_kops_deployment.log
else
	LOG=$LOG_DIR/$(date +"%Y-%m-%d_%H-%M-%S")_$(echo $CLUSTER_SPEC | sed 's|.*/\(.*\)\.yml|\1|')_deployment.log
fi

# Wait until the cluster is up and ready to use
start=$(date +%s)
if create_cluster && kops validate cluster --wait $TIMEOUT | tee -a $LOG; then
	end=$(date +%s)
	runtime=$(date -d "@$((end-start))" "+%M minutes, %S seconds")
	printf "\nSuccessfully deployed cluster in $runtime\n" | tee -a $LOG
	[ -n "$CLUSTER_SPEC" ] && echo "Deployed from cluster spec, \"$CLUSTER_SPEC\"" | tee -a $LOG
	prompt_info "\nTo use your cluster, run the following:"
	prompt_info "export KUBECONFIG=$KUBECONFIG\nexport AWS_PROFILE=$PROFILE"
else
	prompt_error "Failed to deploy cluster"
	exit 1
fi
