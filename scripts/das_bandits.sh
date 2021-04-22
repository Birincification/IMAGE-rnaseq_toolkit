#!/bin/bash

echo $@
params=("$@")

watch pidstat -dru -hl '>>' $log/bandits_salmon-$(date +%s).pidstat & wid=$!

/home/scripts/bandits/das_bandits.R ${params[@]}

kill -15 $wid