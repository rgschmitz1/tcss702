#!/bin/bash

set -o pipefail

# AWS user profile name
PROFILE=eks

# EKS configuration variables
NAME=tcss702-eks
REGION=us-east-2
ZONES=us-east-2b,us-east-2c
NODE_SIZE=c5.2xlarge
NODE_COUNT=1
SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
#NETWORK_CNI=calico
K8S_VERSION=1.23
# Kubernetes and kubectl version must be within one minor version of each other
export KUBECTL_VERSION='v1.23.9'

cd $(dirname $0)

# Add colorized console prompt
. ../utility/color-prompt.sh || exit 1


# Display script usage/flags for user
usage() {
	cat <<- _USAGE
	Usage: $(basename $0) [Options]

	Description: This script will instantiate an EKS kubernetes cluster on AWS

	Options:
	  -d|--delete, deletes a running cluster and all underlying infrastructure
	  -e|--edit, edit cluster configuration before deployment
	  -f|--filename <cluster-spec.yml>, pass a custom cluster spec
	  -h|--help, print this message
	_USAGE
}


# Check for and install dependencies
install_dependencies() {
	../utility/install-kubectl.sh || return $?
	../utility/install-eksctl.sh || return $?
	../utility/install-jq.sh || return $?
	../utility/install-awscli.sh
	return $?
}


# An IAM user must be created in AWS with permissions to create a cluster
create_eks_iam_user() {
	# minimum IAM policies per,
	# https://eksctl.io/usage/minimum-iam-policies
	local managed_policies=(\
		AmazonEC2FullAccess \
		AWSCloudFormationFullAccess \
	)
	# self-managed IAM policies
	local policies=(\
		EksAllAccess \
		IamLimitedAccess \
	)
	# get AWS account ID
	local aws_id=$(aws sts get-caller-identity | jq -r .Account)

	aws iam create-group --group-name $PROFILE || return $?

	# Attach managed policies to group
	for p in ${managed_policies[@]}; do
		aws iam attach-group-policy \
			--policy-arn arn:aws:iam::aws:policy/$p \
			--group-name $PROFILE || return $?
	done

	# Attach self-managed policies to group
	for p in ${policies[@]}; do
		# Check if policy has already been created
		if ! aws iam get-policy --policy-arn arn:aws:iam::$aws_id:policy/$p 2> /dev/null; then
			# policy json
			local json=$(sed "s/<account_id>/$aws_id/" eks/$p.json)
			aws iam create-policy \
				--policy-name $p \
				--policy-document "$json"
		fi
		aws iam attach-group-policy \
			--policy-arn arn:aws:iam::$aws_id:policy/$p \
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


# Create an EKS cluster
create_custer() {
	local ret
	local start=$(date +%s)
	eksctl create cluster \
		--name $NAME \
		--region $REGION \
		--zones $ZONES \
		--instance-types $NODE_SIZE \
		--nodes $NODE_COUNT \
		--version $K8S_VERSION \
		--ssh-access=true \
		--ssh-public-key=$SSH_PUBLIC_KEY \
		--auto-kubeconfig \
		--spot | tee $LOG
	ret=$?
	local end=$(date +%s)
	RUNTIME=$(date -ud "@$((end-start))" "+%M minutes, %S seconds")
	return $ret
}


# Delete EKS cluster
delete_cluster() {
	eksctl delete cluster --name $NAME --region $REGION
	return $?
}


# Check for script dependencies
if ! install_dependencies; then
	prompt_error "failed to install script dependencies"
	exit 1
fi

# Verify eks user has been created
if ! (grep -q $PROFILE $HOME/.aws/credentials || create_eks_iam_user); then
	prompt_error "failed to generate $PROFILE user"
	exit 1
else
	export AWS_PROFILE=$PROFILE
fi

# Check for user input
while [ -n "$1" ]; do
	case $1 in
		-d|--delete)
			if delete_cluster; then
				exit 0
			else
				prompt_error "Encountered an issue deleting cluster"
				exit 1
			fi
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
		*)
			prompt_error "invalid argument!"
			usage
			exit 1
		;;
	esac
done

# Create log for EKS deployment
log_dir=logs/eks
[ -d "$log_dir" ] || mkdir -p $log_dir
LOG=$log_dir/$(date +"%Y-%m-%d_%H-%M-%S")_eks.log

# Create a k8s cluster on AWS EKS
if create_custer; then
	printf "\nSuccessfully deployed cluster in $RUNTIME\n" | tee -a $LOG
	[ -n "$CLUSTER_SPEC" ] && echo "Deployed from cluster spec, \"$CLUSTER_SPEC\"" | tee -a $LOG
else
	prompt_error "failed to deploy cluster"
	exit 1
fi
