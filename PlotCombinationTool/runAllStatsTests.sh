#!/bin/sh

# variable, PS, Mass (Constant Mass)

for variable in Mjj DeltaPhi Etmiss Jet1PT DeltaEta #Jet2PT NumJets Jet1Eta Jet2Eta DeltaEta DeltaPhi
do
    for PS in VBFDM #VBFDM_OR_Monojet # VBFDM_100 Monojet Monojet_HighPt VBFDM_OR_Monojet_HighPt VBFZ_Baseline VBFZ_HighMass VBFZ_Search
    do
	for Dimension in D5a D7a D5a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs Higgs750
	do
	    for statTest in chi2
	    do
		for binError in 0.02 # 0.2 #Lumi2016 Lumi2015 Lumi2023 #0.02
		do
                    ./statisticsTest.sh $variable $PS $Dimension $statTest $binError
		done	
	    done
	done
    done
done


