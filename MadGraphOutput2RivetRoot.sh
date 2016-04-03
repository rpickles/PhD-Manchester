#!/bin/sh


for Dimension in D5a D5b D5c D5d D6a D6b D7a D7b D7c D7d Higgs
do
    for Norm in Normalised Absolute
    do
	if $Dimension = Higgs; then
    #gunzip ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/$Norm/ppchichijj__1GeV_n100000/Events/run_01/unweighted_events.lhe.gz                                          
#cp ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/M1/ppchichijj_$Dimension_1GeV_n100000/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/      
#cd ~/Documents/contrib/lhef2hepmc/                                               
#./lhef2hepmc unweighted_events.lhe ppchichijj_$Dimension_1GeV_unweighted.hepmc                                                                                
	    cp ~/Documents/contrib/lhef2hepmc/ppchichijj_${Dimension}_1GeV_unweighted.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass1/$Norm/
	    cd ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass1/$Norm/
	    rivet -a MC_VBFDM_${Norm} ppchichijj_${Dimension}_1GeV_unweighted.hepmc
	    python removeBits.py Rivet.yoda
	    yoda2root Rivet.yoda
       
	fi
	for Mass in 10 100 1000 
	do
#gunzip ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/$Norm/ppchichijj__$MassGeV_n100000/Events/run_01/unweighted_events.lhe.gz
#cp ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/M$Mass/ppchichijj_$Dimension_$MassGeV_n100000/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
#cd ~/Documents/contrib/lhef2hepmc/ 
#./lhef2hepmc unweighted_events.lhe ppchichijj_$Dimension_$MassGeV_unweighted.hepmc
	cp ~/Documents/contrib/lhef2hepmc/ppchichijj_${Dimension}_${Mass}GeV_unweighted.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	cd ~/Documents/Rivet_Analyses/MC_VBFDM/$Dimension/Mass${Mass}/$Norm/
	rivet -a MC_VBFDM_${Norm} ppchichijj_${Dimension}_${Mass}GeV_unweighted.hepmc
	python removeBits.py Rivet.yoda
	yoda2root Rivet.yoda
	done
    done
done

for Dimension in EWK QCD
do
    for Norm in Normalised Absolute
    do
	#gunzip ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/$Dimension/ppchichijj_n100000/Events/run_01/unweighted_events.lhe.gz
#cp ~/Documents/MG5_aMC_v2_3_2_2/$Dimension/M$Mass/ppchichijj_$Dimension_$MassGeV_n100000/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
#cd ~/Documents/contrib/lhef2hepmc/
#./lhef2hepmc unweighted_events.lhe ppchichijj_$Dimension_$MassGeV_unweighted.hepmc  
	cp ~/Documents/contrib/lhef2hepmc/${Dimension}_Background_unweighted.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
	cd ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
	rivet -a MC_VBFDM_${Norm} ${Dimension}_Background_unweighted.hepmc
	python removeBits.py Rivet.yoda
	yoda2root Rivet.yoda
    done
done

for Dimension in zmumu zvv
do
    for Norm in Normalised Absolute
    do
#        gunzip ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/Scaling/ppzjj_${Dimension}_n10000/Events/run_01/unweighted_events.lhe.gz
	#cp ~/Documents/MG5_aMC_v2_3_2_2/Backgrounds/Scaling/ppzjj_${Dimension}_n10000/Events/run_01/unweighted_events.lhe ~/Documents/contrib/lhef2hepmc/
	#cd ~/Documents/contrib/lhef2hepmc/  
	#./lhef2hepmc unweighted_events.lhe ppzjj_${Dimension}_Background_unweighted.hepmc  
        cp ~/Documents/contrib/lhef2hepmc/ppzjj_${Dimension}_Background_unweighted.hepmc ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
        cd ~/Documents/Rivet_Analyses/MC_VBFDM/Backgrounds/$Dimension/$Norm/
        rivet -a MC_VBFDM_${Norm} ppzjj_${Dimension}_Background_unweighted.hepmc
        python removeBits.py Rivet.yoda
        yoda2root Rivet.yoda
    done
done