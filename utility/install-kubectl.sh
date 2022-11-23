#!/bin/bash
# Check if jq is installed, used to parse kubectl version
$(dirname $0)/install-jq.sh
# Install specified kubectl version
if ! which kubectl > /dev/null || \
	[[ -n "$KUBECTL_VERSION" && "$(kubectl version --output=json 2> /dev/null | jq -r .clientVersion.gitVersion)" != "$KUBECTL_VERSION" ]]; then
	[ -z "$KUBECTL_VERSION" ] && KUBECTL_VERSION=$(curl -Ls https://dl.k8s.io/release/stable.txt)
	curl -LO "https://dl.k8s.io/release/$KUBECTL_VERSION/bin/linux/amd64/kubectl"
	echo $(curl -Ls "https://dl.k8s.io/$KUBECTL_VERSION/bin/linux/amd64/kubectl.sha256") kubectl | \
		sha256sum -c && \
	sudo install -o root -g root -m 0755 kubectl /usr/local/bin/kubectl || \
	echo "checksum incorrect for kubectl"
	rm kubectl
fi
