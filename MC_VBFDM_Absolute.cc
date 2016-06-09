// MC Analysis: Vector Boson Fusion Dark Matter Events

// -*- C++ -*-

/// Initialise and register projections here
#include <boost/assign/list_of.hpp>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

double Count_Total = 0;                            
double Count_VBFZ_Baseline = 0;                                                                 
double Count_VBFZ_HighMass = 0;
double Count_VBFZ_Search = 0;                                                                   
double Count_VBFDM = 0;
double Count_VBFDM_100 = 0;
double Count_Monojet = 0;
double Count_Monojet_HighPt = 0;
double Count_VBFDM_OR_Monojet = 0; 
double Count_VBFDM_OR_Monojet_HighPt = 0;
bool greaterthan200 = false;
double diLeptonPt = 5.;
int muonSize = 6;
double weight = 0.;

namespace Rivet {

  using namespace Cuts;

  class MC_VBFDM_Absolute : public Analysis {
  public:

    /// Constructor
    MC_VBFDM_Absolute()
      : Analysis("MC_VBFDM_Absolute")
    {   }

   public:

    /// Initialize                                                                                                               
      void init() {

	const FinalState fs;

	Cut cuts = (Cuts::pT > 7*GeV) & (Cuts::abseta < 2.5);
        ZFinder zfinder_mu(fs, cuts, PID::MUON, 66*GeV, 116*GeV, 0.1);
        addProjection(zfinder_mu, "ZFinder_mu");

	VetoedFinalState had_fs;
	had_fs.addVetoOnThisFinalState(zfinder_mu);			
	FastJets jets(had_fs, FastJets::ANTIKT, 0.4);
	jets.useInvisibles(true);
	addProjection(jets, "Jets");

	MissingMomentum met(fs);
        addProjection(met, "MET");

	// Bins:
	vector<double> bins_Mjj      = { 200., 400., 600., 1000., 2000., 3000., 4000. };
	vector<double> bins_DeltaEta = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9. };
	vector<double> bins_Jet1PT   = { 50., 100., 150., 200., 250., 300., 400., 500., 600. };
	vector<double> bins_Eta      = { -6., -4., -2., -1., 0., 1., 2., 4., 6. };
	vector<double> bins_DeltaPhi = { 0., pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi };
	vector<double> bins_MET      = { 200., 250., 300., 350., 500., 700., 1000., 1400. };
	vector<double> bins_NumJets  = { 0., 1., 2., 3., 4. };
	vector<double> bins_Count    = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9. };

	/// Histograms:             	
	//Mjj	
	_h["Mjj_PS_VBFZ_Baseline"] = bookHisto1D("Mjj_PS_VBFZ_Baseline", bins_Mjj);
	_h["Mjj_PS_VBFZ_HighMass"] = bookHisto1D("Mjj_PS_VBFZ_HighMass", bins_Mjj);
	_h["Mjj_PS_VBFZ_Search"] = bookHisto1D("Mjj_PS_VBFZ_Search", bins_Mjj);
	_h["Mjj_PS_VBFDM"] = bookHisto1D("Mjj_PS_VBFDM", bins_Mjj);
	_h["Mjj_PS_VBFDM_100"] = bookHisto1D("Mjj_PS_VBFDM_100", bins_Mjj);
	_h["Mjj_PS_Monojet"] = bookHisto1D("Mjj_PS_Monojet", bins_Mjj);
	_h["Mjj_PS_Monojet_HighPt"] = bookHisto1D("Mjj_PS_Monojet_HighPt", bins_Mjj);
	_h["Mjj_PS_VBFDM_OR_Monojet"] = bookHisto1D("Mjj_PS_VBFDM_OR_Monojet", bins_Mjj);
	_h["Mjj_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("Mjj_PS_VBFDM_OR_Monojet_HighPt", bins_Mjj);

	//PT of Jet 1
	_h["Jet1PT_PS_VBFZ_Baseline"] = bookHisto1D("Jet1PT_PS_VBFZ_Baseline", bins_Jet1PT);
	_h["Jet1PT_PS_VBFZ_HighMass"] = bookHisto1D("Jet1PT_PS_VBFZ_HighMass", bins_Jet1PT);
	_h["Jet1PT_PS_VBFZ_Search"] = bookHisto1D("Jet1PT_PS_VBFZ_Search", bins_Jet1PT);
	_h["Jet1PT_PS_VBFDM"] = bookHisto1D("Jet1PT_PS_VBFDM", bins_Jet1PT);
	_h["Jet1PT_PS_VBFDM_100"] = bookHisto1D("Jet1PT_PS_VBFDM_100", bins_Jet1PT);
	_h["Jet1PT_PS_Monojet"] = bookHisto1D("Jet1PT_PS_Monojet", bins_Jet1PT);
	_h["Jet1PT_PS_Monojet_HighPt"] = bookHisto1D("Jet1PT_PS_Monojet_HighPt", bins_Jet1PT);
	_h["Jet1PT_PS_VBFDM_OR_Monojet"] = bookHisto1D("Jet1PT_PS_VBFDM_OR_Monojet", bins_Jet1PT);
	_h["Jet1PT_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("Jet1PT_PS_VBFDM_OR_Monojet_HighPt", bins_Jet1PT);

	//PT of Jet 2
	_h["Jet2PT_PS_VBFZ_Baseline"] = bookHisto1D("Jet2PT_PS_VBFZ_Baseline", bins_Jet1PT);
	_h["Jet2PT_PS_VBFZ_HighMass"] = bookHisto1D("Jet2PT_PS_VBFZ_HighMass", bins_Jet1PT);
	_h["Jet2PT_PS_VBFZ_Search"] = bookHisto1D("Jet2PT_PS_VBFZ_Search", bins_Jet1PT);
	_h["Jet2PT_PS_VBFDM"] = bookHisto1D("Jet2PT_PS_VBFDM", bins_Jet1PT);
	_h["Jet2PT_PS_VBFDM_100"] = bookHisto1D("Jet2PT_PS_VBFDM_100", bins_Jet1PT);
	_h["Jet2PT_PS_Monojet"] = bookHisto1D("Jet2PT_PS_Monojet", bins_Jet1PT);
	_h["Jet2PT_PS_Monojet_HighPt"] = bookHisto1D("Jet2PT_PS_Monojet_HighPt", bins_Jet1PT);
        _h["Jet2PT_PS_VBFDM_OR_Monojet"] = bookHisto1D("Jet2PT_PS_VBFDM_OR_Monojet", bins_Jet1PT);
	_h["Jet2PT_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("Jet2PT_PS_VBFDM_OR_Monojet_HighPt", bins_Jet1PT);

	//Number of Jets
        _h["NumJets_PS_VBFZ_Baseline"] = bookHisto1D("NumJets_PS_VBFZ_Baseline", bins_NumJets);
        _h["NumJets_PS_VBFZ_HighMass"] = bookHisto1D("NumJets_PS_VBFZ_HighMass", bins_NumJets);
        _h["NumJets_PS_VBFZ_Search"] = bookHisto1D("NumJets_PS_VBFZ_Search", bins_NumJets);
        _h["NumJets_PS_VBFDM"] = bookHisto1D("NumJets_PS_VBFDM", bins_NumJets);
        _h["NumJets_PS_VBFDM_100"] = bookHisto1D("NumJets_PS_VBFDM_100", bins_NumJets);
        _h["NumJets_PS_Monojet"] = bookHisto1D("NumJets_PS_Monojet", bins_NumJets);
        _h["NumJets_PS_Monojet_HighPt"] = bookHisto1D("NumJets_PS_Monojet_HighPt", bins_NumJets);
        _h["NumJets_PS_VBFDM_OR_Monojet"] = bookHisto1D("NumJets_PS_VBFDM_OR_Monojet", bins_NumJets);
        _h["NumJets_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("NumJets_PS_VBFDM_OR_Monojet_HighPt", bins_NumJets);

	//Eta Jet 1
	_h["Jet1Eta_PS_VBFZ_Baseline"] = bookHisto1D("Jet1Eta_PS_VBFZ_Baseline", bins_Eta);
        _h["Jet1Eta_PS_VBFZ_HighMass"] = bookHisto1D("Jet1Eta_PS_VBFZ_HighMass", bins_Eta);
        _h["Jet1Eta_PS_VBFZ_Search"] = bookHisto1D("Jet1Eta_PS_VBFZ_Search", bins_Eta);
        _h["Jet1Eta_PS_VBFDM"] = bookHisto1D("Jet1Eta_PS_VBFDM", bins_Eta);
        _h["Jet1Eta_PS_VBFDM_100"] = bookHisto1D("Jet1Eta_PS_VBFDM_100", bins_Eta);
        _h["Jet1Eta_PS_Monojet"] = bookHisto1D("Jet1Eta_PS_Monojet", bins_Eta);
        _h["Jet1Eta_PS_Monojet_HighPt"] = bookHisto1D("Jet1Eta_PS_Monojet_HighPt", bins_Eta);
        _h["Jet1Eta_PS_VBFDM_OR_Monojet"] = bookHisto1D("Jet1Eta_PS_VBFDM_OR_Monojet", bins_Eta);
        _h["Jet1Eta_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("Jet1Eta_PS_VBFDM_OR_Monojet_HighPt", bins_Eta);

	//Eta Jet 2                                                                                                               
        _h["Jet2Eta_PS_VBFZ_Baseline"] = bookHisto1D("Jet2Eta_PS_VBFZ_Baseline", bins_Eta);
        _h["Jet2Eta_PS_VBFZ_HighMass"] = bookHisto1D("Jet2Eta_PS_VBFZ_HighMass", bins_Eta);
        _h["Jet2Eta_PS_VBFZ_Search"] = bookHisto1D("Jet2Eta_PS_VBFZ_Search", bins_Eta);
        _h["Jet2Eta_PS_VBFDM"] = bookHisto1D("Jet2Eta_PS_VBFDM", bins_Eta);
        _h["Jet2Eta_PS_VBFDM_100"] = bookHisto1D("Jet2Eta_PS_VBFDM_100", bins_Eta);
        _h["Jet2Eta_PS_Monojet"] = bookHisto1D("Jet2Eta_PS_Monojet", bins_Eta);
        _h["Jet2Eta_PS_Monojet_HighPt"] = bookHisto1D("Jet2Eta_PS_Monojet_HighPt", bins_Eta);
        _h["Jet2Eta_PS_VBFDM_OR_Monojet"] = bookHisto1D("Jet2Eta_PS_VBFDM_OR_Monojet", bins_Eta);
        _h["Jet2Eta_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("Jet2Eta_PS_VBFDM_OR_Monojet_HighPt", bins_Eta);

	//Delta Eta
        _h["DeltaEta_PS_VBFZ_Baseline"] = bookHisto1D("DeltaEta_PS_VBFZ_Baseline", bins_DeltaEta);
        _h["DeltaEta_PS_VBFZ_HighMass"] = bookHisto1D("DeltaEta_PS_VBFZ_HighMass", bins_DeltaEta);
        _h["DeltaEta_PS_VBFZ_Search"] = bookHisto1D("DeltaEta_PS_VBFZ_Search", bins_DeltaEta);
        _h["DeltaEta_PS_VBFDM"] = bookHisto1D("DeltaEta_PS_VBFDM", bins_DeltaEta);
        _h["DeltaEta_PS_VBFDM_100"] = bookHisto1D("DeltaEta_PS_VBFDM_100", bins_DeltaEta);
        _h["DeltaEta_PS_Monojet"] = bookHisto1D("DeltaEta_PS_Monojet", bins_DeltaEta);
        _h["DeltaEta_PS_Monojet_HighPt"] = bookHisto1D("DeltaEta_PS_Monojet_HighPt", bins_DeltaEta);
        _h["DeltaEta_PS_VBFDM_OR_Monojet"] = bookHisto1D("DeltaEta_PS_VBFDM_OR_Monojet", bins_DeltaEta);
        _h["DeltaEta_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("DeltaEta_PS_VBFDM_OR_Monojet_HighPt", bins_DeltaEta);

	//Delta Phi
        _h["DeltaPhi_PS_VBFZ_Baseline"] = bookHisto1D("DeltaPhi_PS_VBFZ_Baseline", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFZ_HighMass"] = bookHisto1D("DeltaPhi_PS_VBFZ_HighMass", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFZ_Search"] = bookHisto1D("DeltaPhi_PS_VBFZ_Search", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFDM"] = bookHisto1D("DeltaPhi_PS_VBFDM", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFDM_100"] = bookHisto1D("DeltaPhi_PS_VBFDM_100", bins_DeltaPhi);
        _h["DeltaPhi_PS_Monojet"] = bookHisto1D("DeltaPhi_PS_Monojet", bins_DeltaPhi);
        _h["DeltaPhi_PS_Monojet_HighPt"] = bookHisto1D("DeltaPhi_PS_Monojet_HighPt", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFDM_OR_Monojet"] = bookHisto1D("DeltaPhi_PS_VBFDM_OR_Monojet", bins_DeltaPhi);
        _h["DeltaPhi_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("DeltaPhi_PS_VBFDM_OR_Monojet_HighPt", bins_DeltaPhi);

	//Missing Energy
        _h["MET_PS_VBFZ_Baseline"] = bookHisto1D("MET_PS_VBFZ_Baseline", bins_MET);
        _h["MET_PS_VBFZ_HighMass"] = bookHisto1D("MET_PS_VBFZ_HighMass", bins_MET);
        _h["MET_PS_VBFZ_Search"] = bookHisto1D("MET_PS_VBFZ_Search", bins_MET);
        _h["MET_PS_VBFDM"] = bookHisto1D("MET_PS_VBFDM", bins_MET);
        _h["MET_PS_VBFDM_100"] = bookHisto1D("MET_PS_VBFDM_100", bins_MET);
        _h["MET_PS_Monojet"] = bookHisto1D("MET_PS_Monojet", bins_MET);
        _h["MET_PS_Monojet_HighPt"] = bookHisto1D("MET_PS_Monojet_HighPt", bins_MET);
        _h["MET_PS_VBFDM_OR_Monojet"] = bookHisto1D("MET_PS_VBFDM_OR_Monojet", bins_MET);
        _h["MET_PS_VBFDM_OR_Monojet_HighPt"] = bookHisto1D("MET_PS_VBFDM_OR_Monojet_HighPt", bins_MET);

	//Cross Section
	_h["Count_for_All_PS_CrossSection"] = bookHisto1D("Count_for_All_PS_CrossSection", bins_Count);
      }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      weight = event.weight();
      //      MSG_WARNING("weight : "<< weight<<"");
      ++Count_Total;
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");       
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      const ZFinder& zfinder_mu = applyProjection<ZFinder>(event, "ZFinder_mu");
      if (!(zfinder_mu.constituents().size() == 0 || zfinder_mu.constituents().size() == 2)) vetoEvent;
      //      if (zfinder_mu.constituents().size() != 2){ vetoEvent; } 

      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(Cuts::pT >25.0*GeV && Cuts::absrap <2.5)) {
	if ( fabs(jet.momentum().rapidity()) > 4.4 ) continue;
	/*	if (zfinder_mu.constituents().size() == 2){
	  for ( int i = 0 ; i < 2 ; ++i ){
	    if ( fabs(deltaR(jet, zfinder_mu.constituents()[i])) < 0.4 ) continue;
	    }
	    }*/
	jets.push_back(jet.momentum());
      }

      if (zfinder_mu.constituents().size() == 2 ){
	if( zfinder_mu.constituents()[0].pT()>zfinder_mu.constituents()[1].pT() && zfinder_mu.constituents()[0].pT()<25*GeV){ vetoEvent; }
	else if (zfinder_mu.constituents()[1].pT()>zfinder_mu.constituents()[0].pT() && zfinder_mu.constituents()[1].pT()<25*GeV){ vetoEvent; }
	//if(zfinder_mu.constituents()[0].abseta()>2.5 || zfinder_mu.constituents()[1].abseta()>2.5){ vetoEvent; }
      }
    
      vector<FourMomentum> muons;
      if (zfinder_mu.constituents().size() == 2 ){
        for ( int i = 0 ; i < 2 ; ++i ){
          const FourMomentum& zmom = zfinder_mu.constituents()[i].momentum();
          muons.push_back(zmom);
        }
      }

      muonSize = muons.size();

      if (jets.size() < 1) { vetoEvent; }

      double Mjj, Jet1PT, Jet2PT, Jet1Eta, Jet2Eta, DeltaEta, dphi, DeltaPhi1, DeltaPhi, Countforhist, JC, LC, MET;
      
      if (jets.size()>=3){
	JC = ( jets[2].rapidity() - ( (jets[0].rapidity() + jets[1].rapidity()) /2 ) ) / (jets[0].rapidity() - jets[1].rapidity());
	  }
      if (jets.size()>=2){
	Mjj = FourMomentum(jets[0]+jets[1]).mass();
	Jet1PT = jets[0].pT();     
	Jet2PT = jets[1].pT();
	Jet1Eta = jets[0].eta();
	Jet2Eta = jets[1].eta();
	DeltaEta = Jet1Eta-Jet2Eta;
	dphi = fabs(jets[0].phi() - jets[1].phi());
	DeltaPhi1 = ( dphi<=pi ) ? dphi/pi : (2.*pi-dphi)/pi;
	DeltaPhi = DeltaPhi1*pi;
	//	LC = ( met.rapidity() - ( (jets[0].rapidity() + jets[1].rapidity()) /2 ) ) / (jets[0].rapidity() - jets[1].rapidity());
      }
      else{
	Mjj = 0*GeV;
	Jet1PT = jets[0].pT();
	Jet2PT = 0*GeV;
	Jet1Eta = jets[0].eta();
	Jet2Eta = 0*GeV;
	DeltaEta = 0;
	dphi = 0;
	DeltaPhi = 0;
	JC = 0;
	LC = 0;
      }

      if(jets.size()!=1 && Jet2PT<=50*GeV) { vetoEvent; }

      size_t NumJets = jets.size();

      if ( muonSize == 2 ){
      const FourMomentum dilepton = zfinder_mu.constituents()[0].momentum() + zfinder_mu.constituents()[1].momentum();
      MET = dilepton.pT();
      MSG_WARNING("Dilepton pT   = " << dilepton.pT()/GeV << " GeV");
      //if (diLeptonPt>200*GeV){ greaterthan200 = true; }
      }
      else{
      MET = met.vectorEt().mod();
      }
      
      if (MET>200*GeV){ greaterthan200 = true; }

      bool PS_VBFZ_Baseline = ( Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4  );
      bool PS_VBFZ_HighMass = ( Mjj>1000*GeV && Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 );
      bool PS_VBFZ_Search = ( Mjj>250*GeV && Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 );
      bool PS_VBFDM = ( Mjj>200*GeV && Jet1PT>80*GeV && Jet2PT>50*GeV && NumJets>=2 && abseta<4.4 && greaterthan200 );
      bool PS_VBFDM_100 = ( Mjj>200*GeV && Jet1PT>100*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 && greaterthan200 );
      bool PS_Monojet = ( Jet1PT>120*GeV && NumJets>=1 && abseta<4.4 && greaterthan200 );
      bool PS_Monojet_HighPt = ( Jet1PT>250*GeV && NumJets>=1 && greaterthan200 );
      bool PS_VBFDM_OR_Monojet = ( PS_VBFDM || PS_Monojet );
      bool PS_VBFDM_OR_Monojet_HighPt = ( PS_VBFDM || PS_Monojet_HighPt );

      if (PS_VBFZ_Baseline) {
	Countforhist = 0;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
	_h["Mjj_PS_VBFZ_Baseline"]->fill( Mjj , weight );
	_h["Jet1PT_PS_VBFZ_Baseline"]->fill( Jet1PT, weight );
	_h["Jet2PT_PS_VBFZ_Baseline"]->fill( Jet2PT, weight );
	_h["NumJets_PS_VBFZ_Baseline"]->fill( NumJets, weight );
	_h["Jet1Eta_PS_VBFZ_Baseline"]->fill( Jet1Eta, weight );
	_h["Jet2Eta_PS_VBFZ_Baseline"]->fill( Jet2Eta, weight );
	_h["DeltaEta_PS_VBFZ_Baseline"]->fill( DeltaEta, weight );
	_h["DeltaPhi_PS_VBFZ_Baseline"]->fill( DeltaPhi, weight );
	_h["MET_PS_VBFZ_Baseline"]->fill( MET, weight );
	++Count_VBFZ_Baseline;
      }
      if (PS_VBFZ_HighMass) {
	Countforhist = 1;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Mjj_PS_VBFZ_HighMass"]->fill( Mjj, weight );
        _h["Jet1PT_PS_VBFZ_HighMass"]->fill( Jet1PT, weight );
        _h["Jet2PT_PS_VBFZ_HighMass"]->fill( Jet2PT, weight );
        _h["NumJets_PS_VBFZ_HighMass"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFZ_HighMass"]->fill( Jet1Eta, weight );
	_h["Jet2Eta_PS_VBFZ_HighMass"]->fill( Jet2Eta, weight );       
	_h["DeltaEta_PS_VBFZ_HighMass"]->fill( DeltaEta, weight );
        _h["DeltaPhi_PS_VBFZ_HighMass"]->fill( DeltaPhi, weight );
        _h["MET_PS_VBFZ_HighMass"]->fill( MET, weight );
	++Count_VBFZ_HighMass;
      }
      if (PS_VBFZ_Search) {
	Countforhist = 2;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Mjj_PS_VBFZ_Search"]->fill( Mjj, weight );
        _h["Jet1PT_PS_VBFZ_Search"]->fill( Jet1PT, weight );
        _h["Jet2PT_PS_VBFZ_Search"]->fill( Jet2PT, weight );
        _h["NumJets_PS_VBFZ_Search"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFZ_Search"]->fill( Jet1Eta, weight );
        _h["Jet2Eta_PS_VBFZ_Search"]->fill( Jet2Eta, weight );
	_h["DeltaEta_PS_VBFZ_Search"]->fill( DeltaEta, weight );
        _h["DeltaPhi_PS_VBFZ_Search"]->fill( DeltaPhi, weight );
        _h["MET_PS_VBFZ_Search"]->fill( MET, weight );
	++Count_VBFZ_Search;
      }
      if (PS_VBFDM) {
	Countforhist = 3;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Mjj_PS_VBFDM"]->fill( Mjj, weight );
        _h["Jet1PT_PS_VBFDM"]->fill( Jet1PT, weight );
        _h["Jet2PT_PS_VBFDM"]->fill( Jet2PT, weight );
        _h["NumJets_PS_VBFDM"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFDM"]->fill( Jet1Eta, weight );
	_h["Jet2Eta_PS_VBFDM"]->fill( Jet2Eta, weight );        
	_h["DeltaEta_PS_VBFDM"]->fill( DeltaEta, weight );
        _h["DeltaPhi_PS_VBFDM"]->fill( DeltaPhi, weight );
        _h["MET_PS_VBFDM"]->fill( MET, weight );
	++Count_VBFDM;
      }
      if (PS_VBFDM_100) {
        Countforhist = 4;
        _h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Mjj_PS_VBFDM_100"]->fill( Mjj, weight );
        _h["Jet1PT_PS_VBFDM_100"]->fill( Jet1PT, weight );
        _h["Jet2PT_PS_VBFDM_100"]->fill( Jet2PT, weight );
        _h["NumJets_PS_VBFDM_100"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFDM_100"]->fill( Jet1Eta, weight );
        _h["Jet2Eta_PS_VBFDM_100"]->fill( Jet2Eta, weight );
        _h["DeltaEta_PS_VBFDM_100"]->fill( DeltaEta, weight );
        _h["DeltaPhi_PS_VBFDM_100"]->fill( DeltaPhi, weight );
        _h["MET_PS_VBFDM_100"]->fill( MET, weight );
        ++Count_VBFDM_100;
      }
      if (PS_Monojet) {
	Countforhist = 5;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Jet1PT_PS_Monojet"]->fill( Jet1PT, weight );
        _h["NumJets_PS_Monojet"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_Monojet"]->fill( Jet1Eta, weight );
        _h["MET_PS_Monojet"]->fill( MET, weight );       
	if(NumJets>1){
	  _h["Mjj_PS_Monojet"]->fill( Mjj, weight );
	  _h["Jet2PT_PS_Monojet"]->fill( Jet2PT, weight );
	  _h["Jet2Eta_PS_Monojet"]->fill( Jet2Eta, weight );
	  _h["DeltaEta_PS_Monojet"]->fill( DeltaEta, weight );
	  _h["DeltaPhi_PS_Monojet"]->fill( DeltaPhi, weight );
	}
	++Count_Monojet;
      }
      if (PS_Monojet_HighPt) {
        Countforhist = 6;
        _h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Jet1PT_PS_Monojet_HighPt"]->fill( Jet1PT, weight );
        _h["NumJets_PS_Monojet_HighPt"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_Monojet_HighPt"]->fill( Jet1Eta, weight );
        _h["MET_PS_Monojet_HighPt"]->fill( MET, weight );
        if(NumJets>1){
          _h["Mjj_PS_Monojet_HighPt"]->fill( Mjj, weight );
          _h["Jet2PT_PS_Monojet_HighPt"]->fill( Jet2PT, weight );
          _h["Jet2Eta_PS_Monojet_HighPt"]->fill( Jet2Eta, weight );
          _h["DeltaEta_PS_Monojet_HighPt"]->fill( DeltaEta, weight );
          _h["DeltaPhi_PS_Monojet_HighPt"]->fill( DeltaPhi, weight );
        }
        ++Count_Monojet_HighPt;
      }
      if (PS_VBFDM_OR_Monojet) {
	Countforhist = 7;
	_h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Jet1PT_PS_VBFDM_OR_Monojet"]->fill( Jet1PT, weight );
        _h["NumJets_PS_VBFDM_OR_Monojet"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFDM_OR_Monojet"]->fill( Jet1Eta, weight );
        _h["MET_PS_VBFDM_OR_Monojet"]->fill( MET, weight );

	if(NumJets>1){
	  _h["Mjj_PS_VBFDM_OR_Monojet"]->fill( Mjj, weight );
	  _h["Jet2PT_PS_VBFDM_OR_Monojet"]->fill( Jet2PT, weight );
	  _h["Jet2Eta_PS_VBFDM_OR_Monojet"]->fill( Jet2Eta, weight );
	  _h["DeltaEta_PS_VBFDM_OR_Monojet"]->fill( DeltaEta, weight );
	  _h["DeltaPhi_PS_VBFDM_OR_Monojet"]->fill( DeltaPhi, weight );
	  //	  _h["MjjvsDeltaEta_PS_VBFDM_OR_Monojet"]->fill( DeltaEta, Mjj, weight );
	}
	++Count_VBFDM_OR_Monojet;
      }
      if (PS_VBFDM_OR_Monojet_HighPt) {
        Countforhist = 8;
        _h["Count_for_All_PS_CrossSection"]->fill( Countforhist );
        _h["Jet1PT_PS_VBFDM_OR_Monojet_HighPt"]->fill( Jet1PT, weight );
        _h["NumJets_PS_VBFDM_OR_Monojet_HighPt"]->fill( NumJets, weight );
        _h["Jet1Eta_PS_VBFDM_OR_Monojet_HighPt"]->fill( Jet1Eta, weight );
        _h["MET_PS_VBFDM_OR_Monojet_HighPt"]->fill( MET, weight );

        if(NumJets>1){
          _h["Mjj_PS_VBFDM_OR_Monojet_HighPt"]->fill( Mjj, weight );
          _h["Jet2PT_PS_VBFDM_OR_Monojet_HighPt"]->fill( Jet2PT, weight );
          _h["Jet2Eta_PS_VBFDM_OR_Monojet_HighPt"]->fill( Jet2Eta, weight );
          _h["DeltaEta_PS_VBFDM_OR_Monojet_HighPt"]->fill( DeltaEta, weight );
          _h["DeltaPhi_PS_VBFDM_OR_Monojet_HighPt"]->fill( DeltaPhi, weight );
        }
        ++Count_VBFDM_OR_Monojet_HighPt;
      }
    }

    /// Normalise, scale and otherwise manipulate histograms here                     
    void finalize() {
      
      cout<<"\n\nTotal Count = "<<Count_Total;
      cout<<"\nCount with VBFZ_Baseline cuts = "<<Count_VBFZ_Baseline;
      cout<<"\nCount with VBFZ_HighMass cuts = "<<Count_VBFZ_HighMass;
      cout<<"\nCount with VBFZ_Search cuts = "<<Count_VBFZ_Search;
      cout<<"\nCount with VBFDM cuts = "<<Count_VBFDM;
      cout<<"\nCount with VBFDM (Jet1Pt>100GeV) cuts = "<<Count_VBFDM_100;
      cout<<"\nCount with Monojet cuts = "<<Count_Monojet;
      cout<<"\nCount with Monojet_HighPt cuts = "<<Count_Monojet_HighPt;
      cout<<"\nCount with VBFDM_OR_Monojet cuts = "<<Count_VBFDM_OR_Monojet<<endl;
      cout<<"\nCount with VBFDM_OR_Monojet_HighPt cuts = "<<Count_VBFDM_OR_Monojet_HighPt<<endl;
      cout<<"muon size = "<<muonSize<<endl;

      // scale(_h_YYYY, crossSection()/sumOfWeights());

      typedef map<string, Histo1DPtr>::iterator it_type;
      const double sf = crossSection() / sumOfWeights();
      for(it_type it = _h.begin(); it != _h.end(); ++it)  scale(it->second, sf);

      /*
      typedef map<string, Histo1DPtr>::iterator in_type;
      for(in_type in = _h.begin(); in != _h.end(); ++in)  normalize(in->second);
      */
    }

  private:
    // Data members like post-cuts event weight counters go here
    map<string, Histo1DPtr> _h;
   
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_VBFDM_Absolute);
 

}
