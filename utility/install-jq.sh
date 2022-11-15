#!/bin/bash
if ! which jq > /dev/null; then
	sudo apt update && \
	sudo apt install -y jq
	exit $?
fi
