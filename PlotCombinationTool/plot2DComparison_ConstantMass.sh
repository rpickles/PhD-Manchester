#!/bin/sh

if [ "$1" == "" -o "$1" == "-h" ]
    then
    echo ""
    echo "*** Double_t plot2DComparison_ConstantMass(TString Dimension=D5a, TString variable1=Mjj, TString variable2=DeltaEta, TString phasespace=VBFDM, TString Mass=Mass10, TString Norm=Normalised)"
    echo ""
    exit
else 
    cmd="root -l -b -q -x plot2DComparison_ConstantMass.C+'(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\")'"
    eval $cmd 
fi

exit
