#!/bin/sh

if [ "$1" == "" -o "$1" == "-h" ]
    then
    echo ""
    echo "*** Double_t statisticsTest(TString variable=Mjj, TString PS=VBFDM, TString Dimension=D5a, TString StatOption=chi2 )"
    echo ""
    exit
else 
    cmd="root -l -b -q -x statisticsTest.C+'(\"$1\",\"$2\",\"$3\",\"$4\")'"
    eval $cmd 
fi

exit
