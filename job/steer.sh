#!/bin/bash

source /group/clas12/packages/setup.sh
module load clas12/pro

BASEDIR=/work/clas12/fbossu/beamspot/beamspot/
OUTDIR=$(pwd)/output

cp -v $BASEDIR/BeamSpot.java .
cp -v $BASEDIR/run.sh .

ARG=($@)
echo $#
echo ${ARG[@]}
echo ${ARG[@]:1:$#}
./run.sh -B 1 -Z 5 -O BS_$1 ${ARG[@]:1:$#}

mkdir -p $OUTDIR
mv BS_* $OUTDIR
