#!/bin/bash

which docker > /dev/null && exit 0

. $(dirname $0)/color-prompt.sh
prompt_info "Installing docker"

sudo apt-get update && sudo apt-get install -y \
	apt-transport-https \
	ca-certificates \
	curl \
	gnupg-agent \
	software-properties-common || exit $?

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - || exit $?
case $(uname -m) in
	x86_64)
		arch='amd64'
		;;
	aarch64)
		arch='arm64'
		;;
	*)
		echo "ERROR: arch not supported by this script"
		exit 1
		;;
esac
sudo add-apt-repository \
	"deb [arch=$arch] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" || exit $?

# Install docker
sudo apt-get update && sudo apt-get install -y \
	docker-ce docker-ce-cli containerd.io || exit $?

# Verify docker is working
sudo docker run --rm hello-world || exit $?
sudo docker rmi hello-world:latest || exit $?

# Setup so that Docker can be run without sudo
sudo usermod -aG docker `whoami`
exit $?
