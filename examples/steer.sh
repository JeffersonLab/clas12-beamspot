#!/bin/bash

source /group/clas12/packages/setup.sh
module load clas12/pro

BASEDIR=/w/hallb-scifs17exp/clas12/fbossu/beamspot/BeamSpot
OUTDIR=$BASEDIR/output

cp -v $BASEDIR/BeamSpot.java .
cp -v $BASEDIR/run.sh .

./run.sh -B 1 -O BS_$1 $2$1*hipo

mkdir -p $OUTDIR
mv BS_* $OUTDIR
