#!/bin/bash

cd $(dirname $0)
. color-prompt.sh

[ -z "$KOPS_VERSION" ] && KOPS_VERSION='v1.23.4'
which kops > /dev/null && \
	[ "v$(kops version | awk '{print $2}')" = "$KOPS_VERSION" ] && \
	exit 0
prompt_info "Installing kops ($KOPS_VERSION)"
curl -LO https://github.com/kubernetes/kops/releases/download/$KOPS_VERSION/kops-linux-amd64 && \
chmod +x kops-linux-amd64 && \
sudo mv kops-linux-amd64 /usr/local/bin/kops
exit $?
