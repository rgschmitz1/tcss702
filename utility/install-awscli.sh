#!/bin/bash

# Check if awscli is already installed
which aws > /dev/null && exit 0

# Check if unzip is installed
if ! which unzip > /dev/null; then
	sudo apt update
	sudo apt install -y unzip
fi

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip"
cd /tmp
unzip awscliv2.zip
sudo ./aws/install
exit $?
