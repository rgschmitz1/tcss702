#!/bin/bash
[ -z "$1" ] && cluster_name=tcss702-test-bobs || cluster_name=$1
start=$(date +%s)
ibmcloud ks cluster create classic --zone dal10 --flavor b3c.4x16 --hardware shared --public-vlan 3229352 --private-vlan 3229354  --workers 1 --name $cluster_name --version 1.22.9 --public-service-endpoint
cnt=0
until ibmcloud ks worker ls --cluster $cluster_name | grep Ready
do
	echo Counter: $cnt
	((cnt++))
done
end=$(date +%s)
runtime=$((end-start))
ibmcloud ks cluster config --cluster $cluster_name
echo "Start time: $start"
echo "End time: $end"
echo "Total time: $runtime"
