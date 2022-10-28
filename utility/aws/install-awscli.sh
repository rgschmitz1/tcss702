#!/bin/bash

# Check if unzip is installed
which unzip > /dev/null || sudo apt update && sudo apt install unzip

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip"
cd /tmp
unzip awscliv2.zip
sudo ./aws/install
