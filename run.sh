#!/bin/bash

if [[ -z $CLAS12DIR ]];
then 
export CLAS12DIR=$HOME/software/Clas12/coatjava/
fi

echo $CLAS12DIR

rm *.class
javac -Xlint:deprecation -cp "$CLAS12DIR/lib/clas/*:$CLAS12DIR/lib/plugins/*" BeamSpot.java
java -Xmx512m -Xms512m -cp ".:$CLAS12DIR/lib/clas/*:$CLAS12DIR/lib/plugins/*" BeamSpot $*
