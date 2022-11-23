#!/bin/bash

cd $(dirname $0)

# Install arkade
if ! which arkade > /dev/null; then
	curl -SLfs https://get.arkade.dev | sudo sh
fi

# Install helm
if ! which helm > /dev/null; then
	curl -fsSL https://raw.githubusercontent.com/helm/helm/main/scripts/get-helm-3 | sudo bash
fi

# Install faas-cli
if ! which faas-cli > /dev/null; then
	curl -sSL https://cli.openfaas.com | sudo -E sh
fi

# Install minio clent
if ! which mc > /dev/null; then
	wget https://dl.min.io/client/mc/release/linux-amd64/mc
	sudo install -o root -g root -m 0755 mc /usr/local/bin/mc
	rm mc
fi
