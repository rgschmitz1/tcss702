#!/bin/bash
. $(dirname $0)/color-prompt.sh
which jq > /dev/null && exit 0
prompt_info "Installing jq"
sudo apt update && sudo apt install -y jq
exit $?
