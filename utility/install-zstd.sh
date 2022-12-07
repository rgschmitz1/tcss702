#!/bin/bash
. $(dirname $0)/color-prompt.sh
which zstd > /dev/null && exit 0
prompt_info "Installing Zstandard"
sudo apt update && sudo apt install -y zstd
exit $?
