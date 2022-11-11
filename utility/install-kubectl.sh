#!/bin/bash
# Install kubectl
[ -z "$K8S_VERSION" ] && K8S_VERSION=$(curl -Ls https://dl.k8s.io/release/stable.txt)
if ! which kubectl > /dev/null || \
	[ "$(kubectl version --output=json | jq -r .clientVersion.gitVersion)" != "$K8S_VERSION" ]; then
	curl -LO "https://dl.k8s.io/release/$K8S_VERSION/bin/linux/amd64/kubectl"
	echo $(curl -Ls "https://dl.k8s.io/$K8S_VERSION/bin/linux/amd64/kubectl.sha256") kubectl | \
		sha256sum -c && \
	sudo install -o root -g root -m 0755 kubectl /usr/local/bin/kubectl
	rm kubectl
fi
