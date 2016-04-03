#!/bin/sh

# variable, PS, Mass (Constant Mass)

for variable in Mjj # Jet1PT Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi Etmiss
do
    for PS in VBFDM #VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet VBFDM_OR_Monojet_HighPt VBFZ_Baseline VBFZ_HighMass VBFZ_Search
    do
	for Dimension in D5a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs
	do
	    for statTest in chi2
	    do
                ./statisticsTest.sh $variable $PS $Dimension $statTest
	    done
	done
    done
done

