#!/bin/bash
ibmcloud login -u schmitzr1984@gmail.com -r us-south && \
ibmcloud ks cluster config --cluster tcss702-bobs-cluster || \
echo Failed to login to ibmcloud
