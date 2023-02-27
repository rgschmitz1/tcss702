#!/bin/bash
which kind > /dev/null && exit 0
cd $(dirname $0)
. color-prompt.sh
prompt_info "Installing kind (v0.17.0)"
curl -sLo kind "https://kind.sigs.k8s.io/dl/v0.17.0/kind-linux-amd64" && \
chmod +x kind && \
sudo mv kind /usr/local/bin/kind
exit $?
