#!/bin/bash

if [[ -z $CLAS12DIR ]];
then 
    echo CLAS12DIR environment variable must be set to the path of your COATJAVA installation.
    echo Could not find COATJAVA.  Aborting.
    exit
fi

d="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

rm -f $d/*.class
javac -Xlint:deprecation -cp "$CLAS12DIR/lib/clas/*:$CLAS12DIR/lib/plugins/*" $d/BeamSpot.java
java -Xmx512m -Xms512m -cp "$d:$CLAS12DIR/lib/clas/*:$CLAS12DIR/lib/plugins/*" BeamSpot $*

