#! /usr/bin/bash
# USAGE: qdel.sh <queue> <status> (eg. "qdel.sh cygno R")

queue=$1
status=$2
user=$USER

qstat | grep "$status $queue" | grep $user | awk '{print $1}' | awk -F '.' '{print "qdel " $1}'
