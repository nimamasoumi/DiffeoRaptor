#! /bin/bash

execPath="<point to the built Testing directory>"

methodC="diff";
DataDir="<point to your data directory>";
ResDir="<point to your result directory>";

targetSegStruct="<point to the structural segmentation of target>";
targetSegClassic="<point to the tissue segmentation of target>";
outFileDirStruct="<.csv file to save structural Dice>";
outFileDirClassic="<.csv file to save tissue Dice>";

if [ -f $outFileDirStruct ]
then
	rm $outFileDirStruct
fi
if [ -f $outFileDirClassic ]
then
	rm $outFileDirClassic
fi

sourceSegStruct="<point to the structural segmentation of initial source>";
sourceSegClassic="<point to the tissue segmentation of initial source>";

sourceSegStructDeform="<point to the structural segmentation of deformed source>";
sourceSegClassicDeform="<point to the tissue segmentation of deformed source>";

jacD="<point to the determinant of Jacobian file>";

${execPath}/ProcessO3 ${sourceSegStruct} ${targetSegStruct} ${sourceSegStructDeform} ${sourceSegClassic} ${targetSegClassic} ${sourceSegClassicDeform} $jacD $methodC $outFileDirStruct $outFileDirClassic ${oasisCase[$i]}


