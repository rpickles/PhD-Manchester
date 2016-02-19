#!/bin/sh

if [ "$1" == "" -o "$1" == "-h" ]
    then
    echo ""
    echo "*** Double_t plotVBFDMresults_ConstantDimension(TString variable=Mjj, TString phasespace=VBFDM, TString Dimension=D5a, TString Norm=Normalised)"
    echo ""
    exit
else 
    cmd="root -l -b -q -x plotVBFDMresults_ConstantDimension.C+'(\"$1\",\"$2\",\"$3\",\"$4\")'"
    eval $cmd 
fi

exit
