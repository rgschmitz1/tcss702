#!/bin/bash
# install kind (v0.17.0)
if ! which kind > /dev/null; then
	curl -Lo ./kind "https://kind.sigs.k8s.io/dl/v0.17.0/kind-linux-amd64"
	chmod +x ./kind
	sudo mv ./kind /usr/local/bin/kind
fi
