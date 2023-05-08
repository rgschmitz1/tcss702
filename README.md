# TCSS702 - Computer Science and Systems capstone project

## Purpose
This repository contains scripts for setting up a development environment to nvestigate open-source cloud-native application solutions using Kubernetes (K8s).

## Directory structure
| Directory | Description
| --- | ---
| aws | shell scripts for deployment of Kubernetes clusters using 'kops' and 'eksctl'.
| data | payload files used by OpenFaaS functions.
| logs | json logs captured from serverless application analytics framework (SAAF) running in OpenFaaS functions.
| notebooks | Jupyter notebooks used for data analysis/graphing.
| openfaas | OpenFaaS function use cases, including deployment/execution shell scripts.
| src | source code for bash 5.1.x (not needed for distributions using bash 5.1 or newer)
| utility | shell scripts for installation and setup of prerequsite utilites.

## Kubernetes environments
* __kind__ - Kubernetes in docker, used for local development
* __kOps__ - For unmanaged Kubernetes clusters, AWS (Amazon Web Services) and GCE (Google Cloud Platform) are currently officially supported, with DigitalOcean, Hetzner and OpenStack in beta support, and Azure in alpha
* __eks__ - For managed Kubernetes clusters using AWS Elastic Kubernetes Service (EKS)

## Setup
1. Clone the repository and navigate to the root of the directory
```
git clone https://github.com/rgschmitz1/tcss702.git ~/tcss702
cd ~/tcss702
```

2. Install AWS CLI using the following script
```
./utility/install-awscli.sh
```

3. Install kubectl using the following script (this is used to control the Kubernetes cluster)
```
./utility/install-kubectl.sh
```

### kOps setup
1. Install the kOps CLI using the following script
```
./utility/install-kops.sh
```

2. The AWS setup for kOps can be found in the offical documentation here, https://kops.sigs.k8s.io/getting_started/aws
    - Setup a `kops` group and user using the offical documentation.
    - Setup a state store in Amazon S3 using the offical documentation (note that every S3 bucket needs to have a unique name, e.g. `tcss702-rgschmitz-com-state-store`).

3. The following helper script can be used to deploy a Kubernetes cluster using kOps
> __NOTE:__ The deploy script assumes: an SSH key is configured, an AWS user is configured with AWS CLI with profile name `kops`, the (S3) state store (`KOPS_STATE_STORE`) is setup.
```
./aws/deploy-kops-cluster.sh
```
