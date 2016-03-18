#!/bin/bash

rm -rf log

blockMesh | tee -a log
twoPhaseMixingCentralFoam | tee -a log
sample | tee -a log
