#!/bin/sh

CPATH=.
for i in $( ls ./lib/*.jar ); do
        CPATH="$CPATH:$i"
done

echo "Using dynamic classpath: $CPATH"

echo "----------------------------------"


java -version

echo "----------------------------------"

java -cp $CPATH:bin br.embrapa.cnptia.gpbc.plc.ProteinDescriptorCalculation $1 $2 $3 $4 $5 $6 $7 $8
