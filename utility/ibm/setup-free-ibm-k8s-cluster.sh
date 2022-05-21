#!/bin/bash
[ -z "$1" ] && cluster_name=myfreecluster || cluster_name=$1
start=$(date +%s)
ibmcloud ks cluster create classic --name $cluster_name
cnt=0
until ibmcloud ks worker ls --cluster $cluster_name | grep Ready
do
        echo not ready counter: $cnt
        let cnt++
done
end=$(date +%s)
runtime=$((end-start))
#ibmcloud ks cluster config --cluster $cluster_name
echo "Start time: $start"
echo "End time: $end"
echo "Total time: $runtime"
