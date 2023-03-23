# TCSS702 - Computer Science and Systems capstone project

## Purpose
This repository contains scripts for setting up a development environment to nvestigate open-source cloud-native application solutions using Kubernetes (K8s).

# Directory structure
| Directory | Description
| --- | ---
| aws | shell scripts for deployment of Kubernetes clusters using 'kops' and 'eksctl'.
| data | payload files used by OpenFaaS functions.
| logs | json logs captured from serverless application analytics framework (SAAF) running in OpenFaaS functions.
| notebooks | Jupyter notebooks used for data analysis/graphing.
| openfaas | OpenFaaS function use cases, including deployment/execution shell scripts.
| src | source code for bash 5.1.x (not needed for distributions using bash 5.1 or newer)
| utility | shell scripts for installation and setup of prerequsite utilites.

# Kubernetes environements
* __kind__ - Kubernetes in docker, used for local development
* __kOps__ - For unmanaged Kubernetes clusters, AWS (Amazon Web Services) and GCE (Google Cloud Platform) are currently officially supported, with DigitalOcean, Hetzner and OpenStack in beta support, and Azure in alpha
* __eks__ - For managed Kubernetes clusters using AWS Elastic Kubernetes Service (EKS)
