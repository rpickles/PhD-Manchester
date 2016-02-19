#!/bin/sh

if [ "$1" == "" -o "$1" == "-h" ]
    then
    echo ""
    echo "*** Double_t plotVBFDMresults_ConstantMass(TString variable=Mjj, TString phasespace=VBFDM, TString Mass=Mass10, TString Norm=Normalised)"
    echo ""
    exit
else 
    cmd="root -l -b -q -x plotVBFDMresults_ConstantMass.C+'(\"$1\",\"$2\",\"$3\",\"$4\")'"
    eval $cmd 
fi

exit
