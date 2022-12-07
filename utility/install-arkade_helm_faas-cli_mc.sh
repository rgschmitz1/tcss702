#!/bin/bash

cd $(dirname $0)
. color-prompt.sh

if ! which arkade > /dev/null; then
	prompt_info "Installing arkade"
	curl -SLfs https://get.arkade.dev | sudo sh || exit 1
fi

if ! which helm > /dev/null; then
	prompt_info "Installing helm"
	curl -fsSL https://raw.githubusercontent.com/helm/helm/main/scripts/get-helm-3 | sudo bash || exit 1
fi

if ! which faas-cli > /dev/null; then
	prompt_info "Installing faas-cli"
	curl -sSL https://cli.openfaas.com | sudo -E sh || exit 1
fi

which mc > /dev/null && exit 0
prompt_info "Installing minio client"
curl https://dl.min.io/client/mc/release/linux-amd64/mc && \
chmod +x mc && \
sudo mv mc /usr/local/bin/mc
exit $?
