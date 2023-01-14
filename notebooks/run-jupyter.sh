#!/bin/bash

trap cleanup INT

cleanup() {
	prompt_info "Entering cleanup function"
	docker compose down
	exit $1
}

cd $(dirname $0)

. ../utility/color-prompt.sh

../utility/install-docker.sh || exit $?

if [ -z "$(docker ps -f "name=notebook" -q)" ]; then
	prompt_info "Starting Jupyter as daemon"
	docker compose up
else
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
fi
