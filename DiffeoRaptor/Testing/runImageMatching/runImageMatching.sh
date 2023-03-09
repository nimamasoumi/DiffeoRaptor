#!/bin/bash
execPath="BINARY_PATH"
sourcePath="FILE_PATH"
targetPath="FILE_PATH"
numStep=10
lpower=6
truncX=16
truncY=16
truncZ=16 # if 2D, pass z as 1
sigma=0.03
alpha=3.0
gamma=1.0
maxIter=200
stepSizeGD=5.0e-2
mType=0; # 0: Host; 1: DEVICE
 
${execPath}/ImageMatchingTest ${sourcePath} ${targetPath} ${numStep} ${lpower} ${truncX} ${truncY} ${truncZ} ${sigma} ${alpha} ${gamma} ${maxIter} ${stepSizeGD} ${mType}
