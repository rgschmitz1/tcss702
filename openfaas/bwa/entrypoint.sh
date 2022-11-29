#!/bin/bash
# extract bwa index files
tar -xf /home/app/function/genome/GRCh38.d1.vd1_BWA.tar.zst --skip-old-files -C /tmp &
# create empty refernece file, used by bwa for filename only
touch /tmp/GRCh38.d1.vd1.fa &
# start of-watchdog app
exec "$@"
