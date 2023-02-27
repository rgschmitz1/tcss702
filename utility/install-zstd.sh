#!/bin/bash
which zstd > /dev/null && exit 0
. $(dirname $0)/color-prompt.sh
prompt_info "Installing Zstandard"
sudo apt update && sudo apt install -y zstd
exit $?
