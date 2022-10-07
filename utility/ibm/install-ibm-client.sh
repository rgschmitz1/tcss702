#!/bin/bash

# Install ibmcloud and plubin for managed kubernetes
curl -fsSL https://clis.cloud.ibm.com/install/linux | sh
ibmcloud plugin install container-service
