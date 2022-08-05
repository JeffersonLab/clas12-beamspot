#!/bin/bash

F=( $(scontrol -o show job $1 | grep -i failed | cut -d= -f2 | cut -d\  -f1) )

for i in ${F[@]}
do
echo $i
#scontrol requeue $i
done


