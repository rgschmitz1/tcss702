#!/bin/bash

# Remove cached versions of tools
rm -rf $HOME/.arkade

# Install arkade
curl -SLfs https://get.arkade.dev | sudo sh

# Install faas-cli
curl -sSL https://cli.openfaas.com | sudo -E sh

# Install helm
curl -fsSL https://raw.githubusercontent.com/helm/helm/main/scripts/get-helm-3 | sudo bash

# Install minio clent
wget https://dl.min.io/client/mc/release/linux-amd64/mc
sudo install -o root -g root -m 0755 mc /usr/local/bin/mc
rm mc

# Install kubectl
curl -LO "https://dl.k8s.io/release/$(curl -L -s https://dl.k8s.io/release/stable.txt)/bin/linux/amd64/kubectl"
curl -LO "https://dl.k8s.io/$(curl -L -s https://dl.k8s.io/release/stable.txt)/bin/linux/amd64/kubectl.sha256"
echo "$(cat kubectl.sha256)  kubectl" | sha256sum --check
sudo install -o root -g root -m 0755 kubectl /usr/local/bin/kubectl
rm kubectl*
