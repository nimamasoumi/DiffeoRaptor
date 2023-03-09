#! /bin/bash

execPath="<point to the built Testing directory"

numStep=10
lpower=6
truncX=16
truncY=16
truncZ=16 # if 2D, pass z as 1
sigma=0.03 # 0.03
alpha=3.0
gamma=1.0

maxIter0=1;maxIter1=30;maxIter2=30;
stepSizeGD0=0.1;stepSizeGD1=0.1;stepSizeGD2=0.1;
TotalPatch0=500;TotalPatch1=2000;TotalPatch2=4000;
corrFactor0=1000;corrFactor1=20000;corrFactor2=120000;
rrPatchSize0=3;rrPatchSize1=3;rrPatchSize2=3;

outlier_TH=0.1;
numLevel=3;

caseNo='<patient case number>';

DataDir="<point to your data directory>";
ResDir="<point to the result directory>";
ResDir2="<point to the result directory where you want to process the results>";

sourcePath="<point to your source file>";
sourcePathSegStruct="<point to the structurl segmentation of your source file>";
sourcePathSegClassic="<point to the tissue segmentation of your source file>";

targetPath="<point to your template file>";
targetPathSegStruct="<point to the structurl segmentation of your template file>";
targetPathSegClassic="<point to the tissue segmentation of your template file>";

if [ -f "$ResDir/deformIm.nii.gz" ]
then
    rm "$ResDir/deformIm.nii.gz"
fi 

if [ -f "$ResDir/collectedInfoLev*.txt" ]
then
	rm "$ResDir/collectedInfoLev*.txt"
fi

${execPath}/DiffRapO3 ${sourcePath} ${targetPath} ${sourcePathSegStruct} ${targetPathSegStruct} ${sourcePathSegClassic} ${targetPathSegClassic} ${numStep} ${lpower} ${truncX} ${truncY} ${truncZ} ${sigma} ${alpha} ${gamma} ${maxIter0} ${maxIter1} ${maxIter2} ${stepSizeGD0} ${stepSizeGD1} ${stepSizeGD2} ${TotalPatch0} ${TotalPatch1} ${TotalPatch2} ${corrFactor0} ${corrFactor1} ${corrFactor2} ${outlier_TH} ${rrPatchSize0} ${rrPatchSize1} ${rrPatchSize2} ${numLevel}


