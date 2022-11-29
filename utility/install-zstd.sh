#!/bin/bash

# Verify Zstandard is installed
if ! which zstd > /dev/null; then
	sudo apt update && sudo apt install -y zstd || exit 1
fi
