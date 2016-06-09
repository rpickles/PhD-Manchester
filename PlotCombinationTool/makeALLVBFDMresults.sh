#!/bin/sh

# variable, PS, Mass (Constant Mass)

for variable in DeltaEta Mjj DeltaPhi Etmiss Jet1PT #Jet2PT NumJets Jet1Eta Jet2Eta
do
    for PS in VBFDM Monojet #VBFDM_OR_Monojet #VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet_HighPt VBFZ_Baseline VBFZ_HighMass VBFZ_Search
    do
	for Mass in Mass.1 Mass1 Mass10 Mass100 Mass1000
	do
	    for Norm in Normalised
	    do
		./plotVBFDMresults_ConstantMass.sh $variable $PS $Mass $Norm
	    done
	done
    done
done
:'
for variable in Mjj DeltaPhi Etmiss # Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi Etmiss
do
    for PS in VBFDM #Monojet #VBFDM_OR_Monojet #VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet_HighPt VBFZ_Baseline VBFZ_HighMass VBFZ_Search
    do
        for Dimension in D5a D7a D6a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs Higgs750
        do
            for Norm in Absolute #Normalised
            do
	    ./plotVBFDMresults_ConstantDimension.sh $variable $PS $Dimension $Norm
            done
        done
    done
done
'
: '
for Mass in Mass.1 Mass1 Mass10 Mass100 Mass1000 #Mass10000
do
./plotVBFDMresults_ConstantMass.sh Count_for_All CrossSection $Mass Absolute
done
'
