#!/bin/sh

# variable, PS, Mass (Constant Mass)

for variable in Mjj Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi
do
    for PS in VBFZ_Baseline VBFZ_HighMass VBFZ_Search VBFDM VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet VBFDM_OR_Monojet_HighPt 
    do
	for Mass in Mass10 Mass100 Mass1000
	do
	    for Norm in Normalised Absolute
	    do
		./plotVBFDMresults_ConstantMass.sh $variable $PS $Mass $Norm
	    done
	done
    done
done

for variable in Mjj Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi
do
    for PS in VBFZ_Baseline VBFZ_HighMass VBFZ_Search VBFDM VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet VBFDM_OR_Monojet_HighPt
    do
        for Dimension in D5a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs
        do
            for Norm in Normalised Absolute
            do
	    ./plotVBFDMresults_ConstantDimension.sh $variable $PS $Dimension $Norm
            done
        done
    done
done

for Mass in Mass10 Mass100 Mass1000
do
./plotVBFDMresults_ConstantMass.sh Count_for_All CrossSection $Mass Absolute
done

