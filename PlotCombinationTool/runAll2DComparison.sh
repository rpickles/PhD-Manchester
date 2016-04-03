#!/bin/sh

# variable, PS, Mass (Constant Mass)

for variable1 in DeltaEta #Mjj Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaPhi Etmiss
do
    for variable2 in Mjj #Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi Etmiss
    do
	for PS in VBFDM #VBFZ_Baseline VBFZ_HighMass VBFZ_Search VBFDM VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet VBFDM_OR_Monojet_HighPt
	do
            for Dimension in D5a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs
            do
		for Mass in Mass10 Mass100 Mass1000
		do
		    for Norm in Normalised #Absolute
		    do
			./plot2DComparison_ConstantMass.sh $Dimension $variable1 $variable2 $PS $Mass $Norm
		    done
		done
	    done
        done
    done
done
