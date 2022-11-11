#!/bin/bash

set -o pipefail

# Verify all arguments are set
if [ -z "$3" ]; then
	cat <<-_USAGE
	Argumentes include
	\$1 - iterations (set to '$1')
	\$2 - cluster type (set to '$2')
	\$3 - function name (set to '$3')
	_USAGE
	exit 1
fi

# Verify arguments are all set
ITERATION=$1
CLUSTER_TYPE="$2"
FUNCTION_NAME="$3"

# Create a directory to store logs
LOG_DIR="logs/$CLUSTER_TYPE/$FUNCTION_NAME"
[ -d "$LOG_DIR" ] || mkdir -p "$LOG_DIR"

# Add colorized console prompts
. $(dirname $0)/../utility/color-prompt.sh || exit 1

# Execute OpenFaaS
execute_fn() {
	local fn_name="$1"
	local payload="$2"
	local fn_discription="$3"
	local datetime=$(date +"%Y-%m-%d_%H:%M:%S")

	if curl -H "Content-Type: application/json" -X POST -d "$payload" http://localhost:8080/function/$fn_name | \
		tee "$LOG_DIR/${datetime}_${fn_name}_${fn_discription}.log"; then
		sleep 2
	else
		prompt_error "Failed to execute $fn_name $fn_discription"
		exit 1
	fi
}
