#! /usr/bin/bash

DIR=$1
runs=`ls $DIR/reco*_run*_chunk*.root | awk -F "run0" '{print $2}' | awk -F "_" '{print $1}' | uniq`
dim=`ls $DIR/reco*_run*_chunk*.root | awk -F "run0" '{print $2}' | awk -F "_" '{print $2}' | head -n1` # assumes all are uniform (either 2D or 3D)
commandsfile='hadd.sh'
touch $commandsfile
for r in $runs; do
    echo "hadd -k -f reco_run0$r""_$dim.root reco*_run0$r""_$dim""_chunk*.root" >> $commandsfile
done
cat $commandsfile

echo "NOW HADDING FOR REAL..."
#source $commandsfile
echo "done."

