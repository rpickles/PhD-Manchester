#!/bin/sh

for Dimension in D5a # D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs Higgs750
do
    for Norm in Absolute #Normalised
    do
	for Mass in .1 1 10 100 1000 
	do
	    gunzip ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/M$Mass/VBFDM_Dijet_${Dimension}_${Mass}GeV/Events/run_01/unweighted_events.lhe.gz
	    cp ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/M$Mass/VBFDM_Dijet_${Dimension}_${Mass}GeV/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
	    cd ~/Documents/contrib/lhef2hepmc/ 
	    ./lhef2hepmc unweighted_events.lhe VBFDM_Dijet_${Dimension}_${Mass}GeV.hepmc
	    mkdir -p ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	    cp ~/Documents/contrib/lhef2hepmc/VBFDM_Dijet_${Dimension}_${Mass}GeV.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	    cp ~/Documents/Rivet_Analyses/MC_VBFDM/removebits.py ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	    cd ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	    rivet -a MC_VBFDM_${Norm} VBFDM_Dijet_${Dimension}_${Mass}GeV.hepmc
	    python removeBits.py Rivet.yoda
	    yoda2root Rivet.yoda
	done
    done
done

: '
for Dimension in EWK QCD
do
    for Norm in Absolute #Normalised
    do
	gunzip ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/$Dimension/Dijet_${Dimension}_SM/Events/run_01/unweighted_events.lhe.gz
	cp ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/$Dimension/Dijet_${Dimension}_SM/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
	cd ~/Documents/contrib/lhef2hepmc/
	./lhef2hepmc unweighted_events.lhe Dijet_${Dimension}_SM.hepmc  
	cp ~/Documents/contrib/lhef2hepmc/Dijet_${Dimension}_SM.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
	cd ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
	rivet -a MC_VBFDM_${Norm} Dijet_${Dimension}_SM.hepmc
	python removeBits.py Rivet.yoda
	yoda2root Rivet.yoda
    done
done


for Dimension in zvv # zmumu 
do
    for Norm in Absolute #Normalised
    do
	gunzip ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/Scaling/Dijet_${Dimension}_SM/Events/run_01/unweighted_events.lhe.gz
	cp ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/Scaling/Dijet_${Dimension}_SM/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
	cd ~/Documents/contrib/lhef2hepmc/  
	./lhef2hepmc unweighted_events.lhe Dijet_${Dimension}_SM.hepmc  
        cp ~/Documents/contrib/lhef2hepmc/Dijet_${Dimension}_SM.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
        cd ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
        rivet -a MC_VBFDM_${Norm} Dijet_${Dimension}_SM.hepmc
        python removeBits.py Rivet.yoda
        yoda2root Rivet.yoda
    done
done
'