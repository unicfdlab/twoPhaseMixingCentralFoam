#!/bin/bash

rm -rf log
PROCS=`ls -d process*`

for PROC in $PROCS
do
    rm -rf $PROC&
done

rm -rf 0.*

