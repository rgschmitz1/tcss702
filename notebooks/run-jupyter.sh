#!/bin/bash

cd $(dirname $0)/..

. utility/color-prompt.sh

./utility/install-docker.sh || exit $?

if [ -z "$(docker ps -f "name=notebook" -q)" ]; then
	prompt_info "Starting Jupyter as daemon"
	if ! docker run --rm --name notebook -v $PWD:/home/jovyan/work -p 8888:8888 -d jupyter/scipy-notebook:python-3.10.8; then
		prompt_error "Failed to start Jupyter"
		exit 1
	fi
	sleep 3
fi

# Get Jupyter URL
url=$(docker logs notebook 2> /dev/stdout | awk '/127.0.0.1/ {print $NF; exit}')
if [ -z "$url" ]; then
	prompt_error "Failed to get Jupyter URL"
	exit 1
fi
prompt_info "Connect to Jupyter using:"
echo $url

# If xdg-open is installed, open Jupyter URL automatically
which xdg-open > /dev/null && xdg-open $url
