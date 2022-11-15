#!/bin/bash
# Install specified kops version
[ -z "$KOPS_VERSION" ] && KOPS_VERSION='v1.23.4'
if ! which kops > /dev/null || \
	[ "v$(kops version | awk '{print $2}')" != "$KOPS_VERSION" ]; then
	curl -LO https://github.com/kubernetes/kops/releases/download/$KOPS_VERSION/kops-linux-amd64 && \
	sudo install -o root -g root -m 0755 kops-linux-amd64 /usr/local/bin/kops && \
	rm kops-linux-amd64 || exit 1
fi
