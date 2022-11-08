#!/bin/bash

# AWS user profile name
PROFILE=eks

# kops configuration variables
NAME=tcss702-eks.rgschmitz.com
REGION=us-east-2c
#MASTER_SIZE=m5.large
#MASTER_COUNT=1
NODE_SIZE=c5n.2xlarge
NODE_COUNT=1
#SSH_PUBLIC_KEY=$HOME/.ssh/id_rsa.pub
#NETWORK_CNI=calico
#K8S_VERSION=1.23

# The cluster will be deleted if this timeout is exceeded during validation
#TIMEOUT=20m

# Check for and install dependencies
install_dependencies() {
	if ! which jq > /dev/null; then
		sudo apt update
		sudo apt install -y jq || return $?
	fi
	if ! which eksctl > /dev/null; then
		curl --silent --location "https://github.com/weaveworks/eksctl/releases/latest/download/eksctl_$(uname -s)_amd64.tar.gz" | tar xz -C /tmp
		sudo mv /tmp/eksctl /usr/local/bin
		eksctl version 2> /dev/null || return $?
	fi
	$(dirname $0)/install-awscli.sh
	return 0
}

# An IAM user must be created in AWS with permissions to create a cluster
create_eks_iam_user() {
	# minimum IAM policies per,
	# https://eksctl.io/usage/minimum-iam-policies
	local managed_policies=(\
		AmazonEC2FullAccess \
		AWSCloudFormationFullAccess \
	)
	local policies=(\
		EksAllAccess \
		IamLimitedAccess \
	)

	# get AWS account ID
	local aws_id=$(aws sts get-caller-identity | jq -r .Account)

	# check if EksAllAccess policy has been created
	if ! aws iam get-policy --policy-arn arn:aws:iam::$aws_id:policy/EksAllAccess 2> /dev/null; then
		# create policy json
		local eks_all_access=$(cat <<-_POLICY
			{
			    "Version": "2012-10-17",
			    "Statement": [
			        {
			            "Effect": "Allow",
			            "Action": "eks:*",
			            "Resource": "*"
			        },
			        {
			            "Action": [
			                "ssm:GetParameter",
			                "ssm:GetParameters"
			            ],
			            "Resource": [
			                "arn:aws:ssm:*:$aws_id:parameter/aws/*",
			                "arn:aws:ssm:*::parameter/aws/*"
			            ],
			            "Effect": "Allow"
			        },
			        {
			             "Action": [
			               "kms:CreateGrant",
			               "kms:DescribeKey"
			             ],
			             "Resource": "*",
			             "Effect": "Allow"
			        },
			        {
			             "Action": [
			               "logs:PutRetentionPolicy"
			             ],
			             "Resource": "*",
			             "Effect": "Allow"
			        }
			    ]
			}
		_POLICY
		)

		aws iam create-policy \
			--policy-name EksAllAccess \
			--policy-document "$eks_all_access"
	fi

	# check if IamLimitedAccess policy has been created
	if ! aws iam get-policy --policy-arn arn:aws:iam::$aws_id:policy/IamLimitedAccess 2> /dev/null; then
		# Create policy json
		local iam_limited_access=$(cat <<-_POLICY
			{
			    "Version": "2012-10-17",
			    "Statement": [
			        {
			            "Effect": "Allow",
			            "Action": [
			                "iam:CreateInstanceProfile",
			                "iam:DeleteInstanceProfile",
			                "iam:GetInstanceProfile",
			                "iam:RemoveRoleFromInstanceProfile",
			                "iam:GetRole",
			                "iam:CreateRole",
			                "iam:DeleteRole",
			                "iam:AttachRolePolicy",
			                "iam:PutRolePolicy",
			                "iam:ListInstanceProfiles",
			                "iam:AddRoleToInstanceProfile",
			                "iam:ListInstanceProfilesForRole",
			                "iam:PassRole",
			                "iam:DetachRolePolicy",
			                "iam:DeleteRolePolicy",
			                "iam:GetRolePolicy",
			                "iam:GetOpenIDConnectProvider",
			                "iam:CreateOpenIDConnectProvider",
			                "iam:DeleteOpenIDConnectProvider",
			                "iam:TagOpenIDConnectProvider",
			                "iam:ListAttachedRolePolicies",
			                "iam:TagRole",
			                "iam:GetPolicy",
			                "iam:CreatePolicy",
			                "iam:DeletePolicy",
			                "iam:ListPolicyVersions"
			            ],
			            "Resource": [
			                "arn:aws:iam::$aws_id:instance-profile/eksctl-*",
			                "arn:aws:iam::$aws_id:role/eksctl-*",
			                "arn:aws:iam::$aws_id:policy/eksctl-*",
			                "arn:aws:iam::$aws_id:oidc-provider/*",
			                "arn:aws:iam::$aws_id:role/aws-service-role/eks-nodegroup.amazonaws.com/AWSServiceRoleForAmazonEKSNodegroup",
			                "arn:aws:iam::$aws_id:role/eksctl-managed-*"
			            ]
			        },
			        {
			            "Effect": "Allow",
			            "Action": [
			                "iam:GetRole"
			            ],
			            "Resource": [
			                "arn:aws:iam::$aws_id:role/*"
			            ]
			        },
			        {
			            "Effect": "Allow",
			            "Action": [
			                "iam:CreateServiceLinkedRole"
			            ],
			            "Resource": "*",
			            "Condition": {
			                "StringEquals": {
			                    "iam:AWSServiceName": [
			                        "eks.amazonaws.com",
			                        "eks-nodegroup.amazonaws.com",
			                        "eks-fargate.amazonaws.com"
			                    ]
			                }
			            }
			        }
			    ]
			}
		_POLICY
		)

		aws iam create-policy \
			--policy-name IamLimitedAccess \
			--policy-document "$iam_limited_access"
	fi

	aws iam create-group --group-name $PROFILE || return $?

	# Attach managed policies to group
	for p in ${managed_policies[@]}; do
		aws iam attach-group-policy \
			--policy-arn arn:aws:iam::aws:policy/$p \
			--group-name $PROFILE || return $?
	done
	# Attach self-managed policies to group
	for p in ${policies[@]}; do
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
	region=${REGION%?}
	_CONFIG

	# Remove temporary AWS config file
	rm $tmp

	# Verify kops user was successfully created
	aws iam get-user --user $PROFILE
	return $?
}

# Check for script dependencies
if ! install_dependencies; then
	printf "\nERROR: failed to install script dependencies\n"
	exit 1
fi

# Verify eks user has been created
if ! (grep -q $PROFILE $HOME/.aws/credentials || create_eks_iam_user); then
	printf "\nERROR: failed to generate $PROFILE user\n"
	exit 1
else
	export AWS_PROFILE=$PROFILE
fi

echo $?
