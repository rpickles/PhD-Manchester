// MC Analysis: Vector Boson Fusion Dark Matter Events

// -*- C++ -*-

/// Initialise and register projections here
#include <boost/assign/list_of.hpp>
#include "Rivet/Analysis.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
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

namespace Rivet {

  using namespace Cuts;

  class MC_VBFDM_Absolute : public Analysis {
  public:

    /// Constructor
    MC_VBFDM_Absolute()
      : Analysis("MC_VBFDM_Absolute")
    {   }

   public:
   
    std::vector<double> bins_Mjj_PS_VBFZ_Baseline;
    std::vector<double> bins_Mjj_PS_VBFZ_HighMass;
    std::vector<double> bins_Mjj_PS_VBFZ_Search;
    std::vector<double> bins_Mjj_PS_VBFDM;
    std::vector<double> bins_Mjj_PS_Monojet;
    std::vector<double> bins_Mjj_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_Jet1PT_PS_VBFZ_Baseline;
    std::vector<double> bins_Jet1PT_PS_VBFZ_HighMass;
    std::vector<double> bins_Jet1PT_PS_VBFZ_Search;
    std::vector<double> bins_Jet1PT_PS_VBFDM;
    std::vector<double> bins_Jet1PT_PS_Monojet;
    std::vector<double> bins_Jet1PT_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_Jet2PT_PS_VBFZ_Baseline;
    std::vector<double> bins_Jet2PT_PS_VBFZ_HighMass;
    std::vector<double> bins_Jet2PT_PS_VBFZ_Search;
    std::vector<double> bins_Jet2PT_PS_VBFDM;
    std::vector<double> bins_Jet2PT_PS_Monojet;
    std::vector<double> bins_Jet2PT_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_NumJets_PS_VBFZ_Baseline;
    std::vector<double> bins_NumJets_PS_VBFZ_HighMass;
    std::vector<double> bins_NumJets_PS_VBFZ_Search;
    std::vector<double> bins_NumJets_PS_VBFDM;
    std::vector<double> bins_NumJets_PS_Monojet;
    std::vector<double> bins_NumJets_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_Jet1Eta_PS_VBFZ_Baseline;
    std::vector<double> bins_Jet1Eta_PS_VBFZ_HighMass;
    std::vector<double> bins_Jet1Eta_PS_VBFZ_Search;
    std::vector<double> bins_Jet1Eta_PS_VBFDM;
    std::vector<double> bins_Jet1Eta_PS_Monojet;
    std::vector<double> bins_Jet1Eta_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_Jet2Eta_PS_VBFZ_Baseline;
    std::vector<double> bins_Jet2Eta_PS_VBFZ_HighMass;
    std::vector<double> bins_Jet2Eta_PS_VBFZ_Search;
    std::vector<double> bins_Jet2Eta_PS_VBFDM;
    std::vector<double> bins_Jet2Eta_PS_Monojet;
    std::vector<double> bins_Jet2Eta_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_DeltaEta_PS_VBFZ_Baseline;
    std::vector<double> bins_DeltaEta_PS_VBFZ_HighMass;
    std::vector<double> bins_DeltaEta_PS_VBFZ_Search;
    std::vector<double> bins_DeltaEta_PS_VBFDM;
    std::vector<double> bins_DeltaEta_PS_Monojet;
    std::vector<double> bins_DeltaEta_PS_VBFDM_OR_Monojet;
    std::vector<double> bins_DeltaPhi_PS_VBFZ_Baseline;
    std::vector<double> bins_DeltaPhi_PS_VBFZ_HighMass;
    std::vector<double> bins_DeltaPhi_PS_VBFZ_Search;
    std::vector<double> bins_DeltaPhi_PS_VBFDM;
    std::vector<double> bins_DeltaPhi_PS_Monojet;
    std::vector<double> bins_DeltaPhi_PS_VBFDM_OR_Monojet;


    /// Initialize                                                                                                               
      void init() {
	//	FinalState(double mineta, double maxeta, double minpt=0.0*GeV);
	const FinalState fs;
			
	FastJets jets(fs, FastJets::ANTIKT, 0.4);
	jets.useInvisibles();
	addProjection(jets, "Jets");

	MissingMomentum met(FinalState(-5.0,5.0,0*GeV));
        addProjection(met, "MET");

	//Assigning bin widths
	//Mjj
	bins_Mjj_PS_VBFZ_Baseline = boost::assign::list_of<double>(0)(500)(615)(695)(785)(880)(985)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	bins_Mjj_PS_VBFZ_HighMass = boost::assign::list_of<double>(1000)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	bins_Mjj_PS_VBFZ_Search = boost::assign::list_of<double>(250)(500)(615)(695)(785)(880)(985)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	bins_Mjj_PS_VBFDM = boost::assign::list_of<double>(250)(500)(615)(695)(785)(880)(985)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	bins_Mjj_PS_Monojet = boost::assign::list_of<double>(250)(500)(615)(695)(785)(880)(985)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	bins_Mjj_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(250)(500)(615)(695)(785)(880)(985)(1100)(1200)(1350)(1500)(2000)(2700)(4000);
	//PT of Jet 1
        bins_Jet1PT_PS_VBFZ_Baseline = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet1PT_PS_VBFZ_HighMass = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet1PT_PS_VBFZ_Search = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet1PT_PS_VBFDM = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
	bins_Jet1PT_PS_Monojet = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet1PT_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
	//PT of Jet 2
        bins_Jet2PT_PS_VBFZ_Baseline = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet2PT_PS_VBFZ_HighMass = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet2PT_PS_VBFZ_Search = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet2PT_PS_VBFDM = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
	bins_Jet2PT_PS_Monojet = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
        bins_Jet2PT_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(45)(55)(75)(95)(115)(140)(165)(190)(220)(250)(350)(480)(685);
	//Number of Jets
	bins_NumJets_PS_VBFZ_Baseline = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	bins_NumJets_PS_VBFZ_HighMass = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	bins_NumJets_PS_VBFZ_Search = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	bins_NumJets_PS_VBFDM = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	bins_NumJets_PS_Monojet = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	bins_NumJets_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(0)(1)(2)(3)(4);
	//Eta of Jet 1
	bins_Jet1Eta_PS_VBFZ_Baseline = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet1Eta_PS_VBFZ_HighMass = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet1Eta_PS_VBFZ_Search = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet1Eta_PS_VBFDM = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet1Eta_PS_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet1Eta_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	//Eta of Jet 2
	bins_Jet2Eta_PS_VBFZ_Baseline = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet2Eta_PS_VBFZ_HighMass = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet2Eta_PS_VBFZ_Search = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet2Eta_PS_VBFDM = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet2Eta_PS_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_Jet2Eta_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	//Delta Eta
	bins_DeltaEta_PS_VBFZ_Baseline = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_DeltaEta_PS_VBFZ_HighMass = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_DeltaEta_PS_VBFZ_Search = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_DeltaEta_PS_VBFDM = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_DeltaEta_PS_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	bins_DeltaEta_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(-4)(-3)(-2)(-1)(0)(1)(2)(3)(4);
	//Delta Phi
	bins_DeltaPhi_PS_VBFZ_Baseline = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);
        bins_DeltaPhi_PS_VBFZ_HighMass = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);
        bins_DeltaPhi_PS_VBFZ_Search = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);
        bins_DeltaPhi_PS_VBFDM = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);
	bins_DeltaPhi_PS_Monojet = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);
        bins_DeltaPhi_PS_VBFDM_OR_Monojet = boost::assign::list_of<double>(0)(0.6)(0.75)(0.8)(0.85)(0.9)(0.95)(1);

	/// Histograms:             
	//_h_ngapjets_1D_inclusive[i] = bookHisto1D("ngapjets_1D_inclusive"+tags[i], bins_ngapjets_1D_inclusive, "Njets in gap, inclusive region"+labels[i]);
	
	//Mjj	
 	_hist_Mjj_PS_VBFZ_Baseline = bookHisto1D("Mjj_PS_VBFZ_Baseline", 50, 0, 5000, "Mjj_PS_VBFZ_Baseline", "Mjj, GeV", "Probability");
	_hist_Mjj_PS_VBFZ_HighMass = bookHisto1D("Mjj_PS_VBFZ_HighMass", 50, 1000, 5000,"Mjj_PS_VBFZ_HighMass", "Mjj, GeV", "Probability");
	_hist_Mjj_PS_VBFZ_Search = bookHisto1D("Mjj_PS_VBFZ_Search", 50, 250, 5000, "Mjj_PS_VBFZ_Search", "Mjj, GeV", "Probability");
	_hist_Mjj_PS_VBFDM = bookHisto1D("Mjj_PS_VBFDM", 50, 250, 5000, "Mjj_PS_VBFDM", "Mjj, GeV", "Probability");
        _hist_Mjj_PS_VBFDM_100 = bookHisto1D("Mjj_PS_VBFDM_100", 50, 250, 5000, "Mjj_PS_VBFDM_100", "Mjj, GeV", "Probability");	
	_hist_Mjj_PS_Monojet = bookHisto1D("Mjj_PS_Monojet", 50, 250, 5000, "Mjj_PS_Monojet", "Mjj, GeV", "Probability");
        _hist_Mjj_PS_Monojet_HighPt = bookHisto1D("Mjj_PS_Monojet_HighPt", 50, 250, 5000, "Mjj_PS_Monojet_HighPt", "Mjj, GeV", "Probability");
	_hist_Mjj_PS_VBFDM_OR_Monojet = bookHisto1D("Mjj_PS_VBFDM_OR_Monojet", 50, 250, 5000, "Mjj_PS_VBFDM_OR_Monojet", "Mjj, GeV", "Probability");
	_hist_Mjj_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Mjj_PS_VBFDM_OR_Monojet_HighPt", 50, 250, 5000, "Mjj_PS_VBFDM_OR_Monojet_HighPt", "Mjj, GeV", "Probability");
	//PT of Jet 1
        _hist_Jet1PT_PS_VBFZ_Baseline = bookHisto1D("Jet1PT_PS_VBFZ_Baseline", 50, 55, 600, "Jet1PT_PS_VBFZ_Baseline", "Jet1PT, GeV", "Probability");
	_hist_Jet1PT_PS_VBFZ_HighMass = bookHisto1D("Jet1PT_PS_VBFZ_HighMass", 50, 55, 600, "Jet1PT_PS_VBFZ_HighMass", "Jet1PT, GeV", "Probability");
        _hist_Jet1PT_PS_VBFZ_Search = bookHisto1D("Jet1PT_PS_VBFZ_Search", 50, 55, 600, "Jet1PT_PS_VBFZ_Search", "Jet1PT, GeV", "Probability");
        _hist_Jet1PT_PS_VBFDM = bookHisto1D("Jet1PT_PS_VBFDM", 50, 55, 600, "Jet1PT_PS_VBFDM", "Jet1PT, GeV", "Probability");
        _hist_Jet1PT_PS_VBFDM_100 = bookHisto1D("Jet1PT_PS_VBFDM_100", 50, 55, 600, "Jet1PT_PS_VBFDM_100", "Jet1PT, GeV", "Probability");
	_hist_Jet1PT_PS_Monojet = bookHisto1D("Jet1PT_PS_Monojet", 50, 55, 600, "Jet1PT_PS_Monojet", "Jet1PT, GeV", "Probability");
	_hist_Jet1PT_PS_Monojet_HighPt = bookHisto1D("Jet1PT_PS_Monojet_HighPt", 50, 55, 600, "Jet1PT_PS_Monojet_HighPt", "Jet1PT, GeV", "Probability");
	_hist_Jet1PT_PS_VBFDM_OR_Monojet = bookHisto1D("Jet1PT_PS_VBFDM_OR_Monojet", 50, 55, 600, "Jet1PT_PS_VBFDM_OR_Monojet", "Jet1PT, GeV", "Probability");
	_hist_Jet1PT_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Jet1PT_PS_VBFDM_OR_Monojet_HighPt", 50, 55, 600, "Jet1PT_PS_VBFDM_OR_Monojet_HighPt", "Jet1PT, GeV", "Probability");
	//PT of Jet 2
        _hist_Jet2PT_PS_VBFZ_Baseline = bookHisto1D("Jet2PT_PS_VBFZ_Baseline", 50, 45, 600, "Jet2PT_PS_VBFZ_Baseline", "Jet2PT, GeV", "Probability");
	_hist_Jet2PT_PS_VBFZ_HighMass = bookHisto1D("Jet2PT_PS_VBFZ_HighMass", 50, 45, 600, "Jet2PT_PS_VBFZ_HighMass", "Jet2PT, GeV", "Probability");
	_hist_Jet2PT_PS_VBFZ_Search = bookHisto1D("Jet2PT_PS_VBFZ_Search", 50, 45, 600, "Jet2PT_PS_VBFZ_Search", "Jet2PT, GeV", "Probability");
        _hist_Jet2PT_PS_VBFDM = bookHisto1D("Jet2PT_PS_VBFDM", 50, 45, 600, "Jet2PT_PS_VBFDM", "Jet2PT, GeV", "Probability");
        _hist_Jet2PT_PS_VBFDM_100 = bookHisto1D("Jet2PT_PS_VBFDM_100", 50, 45, 600, "Jet2PT_PS_VBFDM_100", "Jet2PT, GeV", "Probability");
	_hist_Jet2PT_PS_Monojet = bookHisto1D("Jet2PT_PS_Monojet", 50, 45, 600, "Jet2PT_PS_Monojet", "Jet2PT, GeV", "Probability");
	_hist_Jet2PT_PS_Monojet_HighPt = bookHisto1D("Jet2PT_PS_Monojet_HighPt", 50, 45, 600, "Jet2PT_PS_Monojet_HighPt", "Jet2PT, GeV", "Probability");
        _hist_Jet2PT_PS_VBFDM_OR_Monojet = bookHisto1D("Jet2PT_PS_VBFDM_OR_Monojet", 50, 45, 600, "Jet2PT_PS_VBFDM_OR_Monojet", "Jet2PT, GeV", "Probability");
	_hist_Jet2PT_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Jet2PT_PS_VBFDM_OR_Monojet_HighPt", 50, 45, 600, "Jet2PT_PS_VBFDM_OR_Monojet_HighPt", "Jet2PT, GeV", "Probability");
	//Number of Jets
        _hist_NumJets_PS_VBFZ_Baseline = bookHisto1D("NumJets_PS_VBFZ_Baseline", 5, 0, 4, "NumJets_PS_VBFZ_Baseline", "NumJets", "Probability");
	_hist_NumJets_PS_VBFZ_HighMass = bookHisto1D("NumJets_PS_VBFZ_HighMass", 5, 0, 4, "NumJets_PS_VBFZ_HighMass", "NumJets", "Probability");
	_hist_NumJets_PS_VBFZ_Search = bookHisto1D("NumJets_PS_VBFZ_Search", 5, 0, 4, "NumJets_PS_VBFZ_Search", "NumJets ", "Probability");
        _hist_NumJets_PS_VBFDM = bookHisto1D("NumJets_PS_VBFDM", 5, 0, 4, "NumJets_PS_VBFDM", "NumJets ", "Probability");
        _hist_NumJets_PS_VBFDM_100 = bookHisto1D("NumJets_PS_VBFDM_100", 5, 0, 4, "NumJets_PS_VBFDM_100", "NumJets ", "Probability");
        _hist_NumJets_PS_Monojet = bookHisto1D("NumJets_PS_Monojet", 5, 0, 4, "NumJets_PS_Monojet", "NumJets ", "Probability");
	_hist_NumJets_PS_Monojet_HighPt = bookHisto1D("NumJets_PS_Monojet_HighPt", 5, 0, 4, "NumJets_PS_Monojet_HighPt", "NumJets ", "Probability");
        _hist_NumJets_PS_VBFDM_OR_Monojet = bookHisto1D("NumJets_PS_VBFDM_OR_Monojet", 5, 0, 4,"NumJets_PS_VBFDM_OR_Monojet", "NumJets ", "Probability");
	_hist_NumJets_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("NumJets_PS_VBFDM_OR_Monojet_HighPt", 5, 0, 4,"NumJets_PS_VBFDM_OR_Monojet_HighPt", "NumJets ", "Probability");
	//Eta Jet 1
        _hist_Jet1Eta_PS_VBFZ_Baseline = bookHisto1D("Jet1Eta_PS_VBFZ_Baseline", 50, -6, 6, "Jet1Eta_PS_VBFZ_Baseline", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_VBFZ_HighMass = bookHisto1D("Jet1Eta_PS_VBFZ_HighMass", 50, -6, 6, "Jet1Eta_PS_VBFZ_HighMass", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_VBFZ_Search = bookHisto1D("Jet1Eta_PS_VBFZ_Search", 50, -6, 6, "Jet1Eta_PS_VBFZ_Search", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_VBFDM = bookHisto1D("Jet1Eta_PS_VBFDM", 50, -6, 6, "Jet1Eta_PS_VBFDM", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_VBFDM_100 = bookHisto1D("Jet1Eta_PS_VBFDM_100", 50, -6, 6, "Jet1Eta_PS_VBFDM_100", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_Monojet = bookHisto1D("Jet1Eta_PS_Monojet", 50, -6, 6, "Jet1Eta_PS_Monojet", "Jet1Eta ", "Probability");
	_hist_Jet1Eta_PS_Monojet_HighPt = bookHisto1D("Jet1Eta_PS_Monojet_HighPt", 50, -6, 6, "Jet1Eta_PS_Monojet_HighPt", "Jet1Eta ", "Probability");
        _hist_Jet1Eta_PS_VBFDM_OR_Monojet = bookHisto1D("Jet1Eta_PS_VBFDM_OR_Monojet", 50, -6, 6, "Jet1Eta_PS_VBFDM_OR_Monojet", "Jet1Eta ", "Probability");
	_hist_Jet1Eta_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Jet1Eta_PS_VBFDM_OR_Monojet_HighPt", 50, -6, 6, "Jet1Eta_PS_VBFDM_OR_Monojet_HighPt", "Jet1Eta ", "Probability");
	//Eta Jet 2                                                                                                               
        _hist_Jet2Eta_PS_VBFZ_Baseline = bookHisto1D("Jet2Eta_PS_VBFZ_Baseline", 50, -6, 6, "Jet2Eta_PS_VBFZ_Baseline", "Jet2Eta ", "Probability");
        _hist_Jet2Eta_PS_VBFZ_HighMass = bookHisto1D("Jet2Eta_PS_VBFZ_HighMass", 50, -6, 6, "Jet2Eta_PS_VBFZ_HighMass", "Jet2Eta ", "Probability");
        _hist_Jet2Eta_PS_VBFZ_Search = bookHisto1D("Jet2Eta_PS_VBFZ_Search", 50, -6, 6, "Jet2Eta_PS_VBFZ_Search", "Jet2Eta ", "Probability");
        _hist_Jet2Eta_PS_VBFDM = bookHisto1D("Jet2Eta_PS_VBFDM", 50, -6, 6, "Jet2Eta_PS_VBFDM", "Jet2Eta ", "Probability");
        _hist_Jet2Eta_PS_VBFDM_100 = bookHisto1D("Jet2Eta_PS_VBFDM_100", 50, -6, 6, "Jet2Eta_PS_VBFDM_100", "Jet2Eta ", "Probability");
	_hist_Jet2Eta_PS_Monojet = bookHisto1D("Jet2Eta_PS_Monojet", 50, -6, 6, "Jet2Eta_PS_Monojet", "Jet2Eta ", "Probability");
	_hist_Jet2Eta_PS_Monojet_HighPt = bookHisto1D("Jet2Eta_PS_Monojet_HighPt", 50, -6, 6, "Jet2Eta_PS_Monojet_HighPt", "Jet2Eta ", "Probability");
        _hist_Jet2Eta_PS_VBFDM_OR_Monojet = bookHisto1D("Jet2Eta_PS_VBFDM_OR_Monojet", 50, -6, 6, "Jet2Eta_PS_VBFDM_OR_Monojet", "Jet2Eta ", "Probability");
	_hist_Jet2Eta_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Jet2Eta_PS_VBFDM_OR_Monojet_HighPt", 50, -6, 6, "Jet2Eta_PS_VBFDM_OR_Monojet_HighPt", "Jet2Eta ", "Probability");
	//Delta Eta
        _hist_DeltaEta_PS_VBFZ_Baseline = bookHisto1D("DeltaEta_PS_VBFZ_Baseline", 50, 0, 9, "DeltaEta_PS_VBFZ_Baseline", "DeltaEta", "Probability");
	_hist_DeltaEta_PS_VBFZ_HighMass = bookHisto1D("DeltaEta_PS_VBFZ_HighMass", 50, 0, 9, "DeltaEta_PS_VBFZ_HighMass", "DeltaEta ", "Probability");
        _hist_DeltaEta_PS_VBFZ_Search = bookHisto1D("DeltaEta_PS_VBFZ_Search", 50, 0, 9, "DeltaEta_PS_VBFZ_Search", "DeltaEta ", "Probability");
        _hist_DeltaEta_PS_VBFDM = bookHisto1D("DeltaEta_PS_VBFDM", 50, 0, 9, "DeltaEta_PS_VBFDM", "DeltaEta ", "Probability");
        _hist_DeltaEta_PS_VBFDM_100 = bookHisto1D("DeltaEta_PS_VBFDM_100", 50, 0, 9, "DeltaEta_PS_VBFDM_100", "DeltaEta ", "Probability");
	_hist_DeltaEta_PS_Monojet = bookHisto1D("DeltaEta_PS_Monojet", 50, 0, 9, "DeltaEta_PS_Monojet", "DeltaEta ", "Probability");
	_hist_DeltaEta_PS_Monojet_HighPt = bookHisto1D("DeltaEta_PS_Monojet_HighPt", 50, 0, 9, "DeltaEta_PS_Monojet_HighPt", "DeltaEta ", "Probability");
        _hist_DeltaEta_PS_VBFDM_OR_Monojet = bookHisto1D("DeltaEta_PS_VBFDM_OR_Monojet", 50, 0, 9, "DeltaEta_PS_VBFDM_OR_Monojet", "DeltaEta ", "Probability");
	_hist_DeltaEta_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("DeltaEta_PS_VBFDM_OR_Monojet_HighPt", 50, 0, 9, "DeltaEta_PS_VBFDM_OR_Monojet_HighPt", "DeltaEta ", "Probability");
	//Delta Phi
        _hist_DeltaPhi_PS_VBFZ_Baseline = bookHisto1D("DeltaPhi_PS_VBFZ_Baseline", 50, 0, 1, "Mjj_PS_VBFZ_Baseline", "DeltaPhi ", "Probability");
        _hist_DeltaPhi_PS_VBFZ_HighMass = bookHisto1D("DeltaPhi_PS_VBFZ_HighMass", 50, 0, 1, "DeltaPhi_PS_VBFZ_HighMass", "DeltaPhi ", "Probability");
        _hist_DeltaPhi_PS_VBFZ_Search = bookHisto1D("DeltaPhi_PS_VBFZ_Search", 50, 0, 1, "DeltaPhi_PS_VBFZ_Search", "DeltaPhi ", "Probability");
        _hist_DeltaPhi_PS_VBFDM = bookHisto1D("DeltaPhi_PS_VBFDM", 50, 0, 1, "DeltaPhi_PS_VBFDM", "DeltaPhi ", "Probability");
        _hist_DeltaPhi_PS_VBFDM_100 = bookHisto1D("DeltaPhi_PS_VBFDM_100", 50, 0, 1, "DeltaPhi_PS_VBFDM_100", "DeltaPhi ", "Probability");
	_hist_DeltaPhi_PS_Monojet = bookHisto1D("DeltaPhi_PS_Monojet", 50, 0, 1, "DeltaPhi_PS_Monojet", "DeltaPhi ", "Probability");
	_hist_DeltaPhi_PS_Monojet_HighPt = bookHisto1D("DeltaPhi_PS_Monojet_HighPt", 50, 0, 1, "DeltaPhi_PS_Monojet_HighPt", "DeltaPhi ", "Probability");
        _hist_DeltaPhi_PS_VBFDM_OR_Monojet = bookHisto1D("DeltaPhi_PS_VBFDM_OR_Monojet", 50, 0, 1, "DeltaPhi_PS_VBFDM_OR_Monojet", "DeltaPhi ", "Probability");
	_hist_DeltaPhi_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("DeltaPhi_PS_VBFDM_OR_Monojet_HighPt", 50, 0, 1, "DeltaPhi_PS_VBFDM_OR_Monojet_HighPt", "DeltaPhi ", "Probability");
	//Missing Energy
	_hist_met_PS_VBFZ_Baseline = bookHisto1D("Etmiss_PS_VBFZ_Baseline", 50, 0.0, 1500, "ETMiss_PS_VBFZ_Baseline", "ETMiss, GeV", "Probability");
	_hist_met_PS_VBFZ_HighMass = bookHisto1D("Etmiss_PS_VBFZ_HighMass", 50, 0.0, 1500, "ETMiss_PS_VBFZ_HighMass", "ETMiss, GeV", "Probability");
	_hist_met_PS_VBFZ_Search = bookHisto1D("Etmiss_PS_VBFZ_Search", 50, 0.0, 1500, "ETMiss_PS_VBFZ_Search", "ETMiss, GeV", "Probability");
	_hist_met_PS_VBFDM = bookHisto1D("Etmiss_PS_VBFDM", 50, 0.0, 1500, "ETMiss_PS_VBFDM", "ETMiss, GeV", "Probability");
        _hist_met_PS_VBFDM_100 = bookHisto1D("Etmiss_PS_VBFDM_100", 50, 0.0, 1500, "ETMiss_PS_VBFDM_100", "ETMiss, GeV", "Probability");
	_hist_met_PS_Monojet = bookHisto1D("Etmiss_PS_Monojet", 50, 0.0, 1500, "ETMiss_PS_Monojet", "ETMiss, GeV", "Probability");
	_hist_met_PS_Monojet_HighPt = bookHisto1D("Etmiss_PS_Monojet_HighPt", 50, 0.0, 1500, "ETMiss_PS_Monojet_HighPt", "ETMiss, GeV", "Probability");
	_hist_met_PS_VBFDM_OR_Monojet = bookHisto1D("Etmiss_PS_VBFDM_OR_Monojet", 50, 0.0, 1500, "ETMiss_PS_VBFDM_OR_Monojet", "ETMiss, GeV", "Probability");
	_hist_met_PS_VBFDM_OR_Monojet_HighPt = bookHisto1D("Etmiss_PS_VBFDM_OR_Monojet_HighPt", 50, 0.0, 1500, "ETMiss_PS_VBFDM_OR_Monojet_HighPt", "ETMiss, GeV", "Probability");
	//Cross Section
        _hist_Count_for_All_PS_CrossSection = bookHisto1D("Count_for_All_PS_CrossSection", 9, 0, 9, "Count_for_All_PS_CrossSection", "Phase Space", "CrossSectionxAcceptedCounts/TotalCounts");

      }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      ++Count_Total;
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");       
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");

      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(45.0*GeV)) {
	if ( fabs(jet.momentum().rapidity()) > 4.4 ) continue;
	//	if ( fabs(deltaR(jet, lepton)) < 0.3 ) continue;
	jets.push_back(jet.momentum());
      }

      if (jets.size() < 1) { vetoEvent; }

      double Mjj, Jet1PT, Jet2PT, Jet1Eta, Jet2Eta, DeltaEta, dphi, DeltaPhi, Countforhist, JC, LC;
      
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
	DeltaPhi = ( dphi<=pi ) ? dphi/pi : (2.*pi-dphi)/pi;
	LC = ( met.rapidity() - ( (jets[0].rapidity() + jets[1].rapidity()) /2 ) ) / (jets[0].rapidity() - jets[1].rapidity());
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

      if(Jet2PT!=0*GeV && Jet2PT<=45*GeV) { vetoEvent; }

      size_t NumJets = jets.size();
      double MET = met.vectorEt().mod();

      bool PS_VBFZ_Baseline = ( Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4  );
      bool PS_VBFZ_HighMass = ( Mjj>1000*GeV && Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 );
      bool PS_VBFZ_Search = ( Mjj>250*GeV && Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 );
      bool PS_VBFDM = ( Mjj>250*GeV && Jet1PT>55*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 && MET>150*GeV );
      bool PS_VBFDM_100 = ( Mjj>250*GeV && Jet1PT>100*GeV && Jet2PT>45*GeV && NumJets>=2 && abseta<4.4 && MET>150*GeV );
      bool PS_Monojet = ( Jet1PT>100*GeV && NumJets>=1 && abseta<4.4 && MET>150*GeV );
      bool PS_Monojet_HighPt = ( Jet1PT>250*GeV && NumJets>=1 && MET>250*GeV );
      bool PS_VBFDM_OR_Monojet = ( PS_VBFDM || PS_Monojet );
      bool PS_VBFDM_OR_Monojet_HighPt = ( PS_VBFDM || PS_Monojet_HighPt );

      if (PS_VBFZ_Baseline) {
	Countforhist = 0;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
	_hist_Mjj_PS_VBFZ_Baseline->fill( Mjj , weight );
	_hist_Jet1PT_PS_VBFZ_Baseline->fill( Jet1PT, weight );
	_hist_Jet2PT_PS_VBFZ_Baseline->fill( Jet2PT, weight );
	_hist_NumJets_PS_VBFZ_Baseline->fill( NumJets, weight );
	_hist_Jet1Eta_PS_VBFZ_Baseline->fill( Jet1Eta, weight );
	_hist_Jet2Eta_PS_VBFZ_Baseline->fill( Jet2Eta, weight );
	_hist_DeltaEta_PS_VBFZ_Baseline->fill( DeltaEta, weight );
	_hist_DeltaPhi_PS_VBFZ_Baseline->fill( DeltaPhi, weight );
	_hist_met_PS_VBFZ_Baseline->fill( MET, weight );

	++Count_VBFZ_Baseline;
      }
      if (PS_VBFZ_HighMass) {
	Countforhist = 1;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Mjj_PS_VBFZ_HighMass->fill( Mjj, weight );
        _hist_Jet1PT_PS_VBFZ_HighMass->fill( Jet1PT, weight );
        _hist_Jet2PT_PS_VBFZ_HighMass->fill( Jet2PT, weight );
        _hist_NumJets_PS_VBFZ_HighMass->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFZ_HighMass->fill( Jet1Eta, weight );
	_hist_Jet2Eta_PS_VBFZ_HighMass->fill( Jet2Eta, weight );       
	_hist_DeltaEta_PS_VBFZ_HighMass->fill( DeltaEta, weight );
        _hist_DeltaPhi_PS_VBFZ_HighMass->fill( DeltaPhi, weight );
        _hist_met_PS_VBFZ_HighMass->fill( MET, weight );

	++Count_VBFZ_HighMass;
      }
      if (PS_VBFZ_Search) {
	Countforhist = 2;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Mjj_PS_VBFZ_Search->fill( Mjj, weight );
        _hist_Jet1PT_PS_VBFZ_Search->fill( Jet1PT, weight );
        _hist_Jet2PT_PS_VBFZ_Search->fill( Jet2PT, weight );
        _hist_NumJets_PS_VBFZ_Search->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFZ_Search->fill( Jet1Eta, weight );
        _hist_Jet2Eta_PS_VBFZ_Search->fill( Jet2Eta, weight );
	_hist_DeltaEta_PS_VBFZ_Search->fill( DeltaEta, weight );
        _hist_DeltaPhi_PS_VBFZ_Search->fill( DeltaPhi, weight );
        _hist_met_PS_VBFZ_Search->fill( MET, weight );

	++Count_VBFZ_Search;
      }
      if (PS_VBFDM) {
	Countforhist = 3;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Mjj_PS_VBFDM->fill( Mjj, weight );
        _hist_Jet1PT_PS_VBFDM->fill( Jet1PT, weight );
        _hist_Jet2PT_PS_VBFDM->fill( Jet2PT, weight );
        _hist_NumJets_PS_VBFDM->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFDM->fill( Jet1Eta, weight );
	_hist_Jet2Eta_PS_VBFDM->fill( Jet2Eta, weight );        
	_hist_DeltaEta_PS_VBFDM->fill( DeltaEta, weight );
        _hist_DeltaPhi_PS_VBFDM->fill( DeltaPhi, weight );
        _hist_met_PS_VBFDM->fill( MET, weight );

	++Count_VBFDM;
      }
      if (PS_VBFDM_100) {
        Countforhist = 4;
        _hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Mjj_PS_VBFDM_100->fill( Mjj, weight );
        _hist_Jet1PT_PS_VBFDM_100->fill( Jet1PT, weight );
        _hist_Jet2PT_PS_VBFDM_100->fill( Jet2PT, weight );
        _hist_NumJets_PS_VBFDM_100->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFDM_100->fill( Jet1Eta, weight );
        _hist_Jet2Eta_PS_VBFDM_100->fill( Jet2Eta, weight );
        _hist_DeltaEta_PS_VBFDM_100->fill( DeltaEta, weight );
        _hist_DeltaPhi_PS_VBFDM_100->fill( DeltaPhi, weight );
        _hist_met_PS_VBFDM_100->fill( MET, weight );

        ++Count_VBFDM_100;
      }

      if (PS_Monojet) {
	Countforhist = 5;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Jet1PT_PS_Monojet->fill( Jet1PT, weight );
        _hist_NumJets_PS_Monojet->fill( NumJets, weight );
        _hist_Jet1Eta_PS_Monojet->fill( Jet1Eta, weight );
        _hist_met_PS_Monojet->fill( MET, weight );
	
	if(NumJets>1){
	  _hist_Mjj_PS_Monojet->fill( Mjj, weight );
	  _hist_Jet2PT_PS_Monojet->fill( Jet2PT, weight );
	  _hist_Jet2Eta_PS_Monojet->fill( Jet2Eta, weight );
	  _hist_DeltaEta_PS_Monojet->fill( DeltaEta, weight );
	  _hist_DeltaPhi_PS_Monojet->fill( DeltaPhi, weight );
	}
       
	++Count_Monojet;
      }
      if (PS_Monojet_HighPt) {
        Countforhist = 6;
        _hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Jet1PT_PS_Monojet_HighPt->fill( Jet1PT, weight );
        _hist_NumJets_PS_Monojet_HighPt->fill( NumJets, weight );
        _hist_Jet1Eta_PS_Monojet_HighPt->fill( Jet1Eta, weight );
        _hist_met_PS_Monojet_HighPt->fill( MET, weight );

        if(NumJets>1){
          _hist_Mjj_PS_Monojet_HighPt->fill( Mjj, weight );
          _hist_Jet2PT_PS_Monojet_HighPt->fill( Jet2PT, weight );
          _hist_Jet2Eta_PS_Monojet_HighPt->fill( Jet2Eta, weight );
          _hist_DeltaEta_PS_Monojet_HighPt->fill( DeltaEta, weight );
          _hist_DeltaPhi_PS_Monojet_HighPt->fill( DeltaPhi, weight );
        }

        ++Count_Monojet_HighPt;
      }
      if (PS_VBFDM_OR_Monojet) {
	Countforhist = 7;
	_hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Jet1PT_PS_VBFDM_OR_Monojet->fill( Jet1PT, weight );
        _hist_NumJets_PS_VBFDM_OR_Monojet->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFDM_OR_Monojet->fill( Jet1Eta, weight );
        _hist_met_PS_VBFDM_OR_Monojet->fill( MET, weight );

	if(NumJets>1){
	  _hist_Mjj_PS_VBFDM_OR_Monojet->fill( Mjj, weight );
	  _hist_Jet2PT_PS_VBFDM_OR_Monojet->fill( Jet2PT, weight );
	  _hist_Jet2Eta_PS_VBFDM_OR_Monojet->fill( Jet2Eta, weight );
	  _hist_DeltaEta_PS_VBFDM_OR_Monojet->fill( DeltaEta, weight );
	  _hist_DeltaPhi_PS_VBFDM_OR_Monojet->fill( DeltaPhi, weight );
	}
	++Count_VBFDM_OR_Monojet;
      }
      if (PS_VBFDM_OR_Monojet_HighPt) {
        Countforhist = 8;
        _hist_Count_for_All_PS_CrossSection->fill( Countforhist );
        _hist_Jet1PT_PS_VBFDM_OR_Monojet_HighPt->fill( Jet1PT, weight );
        _hist_NumJets_PS_VBFDM_OR_Monojet_HighPt->fill( NumJets, weight );
        _hist_Jet1Eta_PS_VBFDM_OR_Monojet_HighPt->fill( Jet1Eta, weight );
        _hist_met_PS_VBFDM_OR_Monojet_HighPt->fill( MET, weight );

        if(NumJets>1){
          _hist_Mjj_PS_VBFDM_OR_Monojet_HighPt->fill( Mjj, weight );
          _hist_Jet2PT_PS_VBFDM_OR_Monojet_HighPt->fill( Jet2PT, weight );
          _hist_Jet2Eta_PS_VBFDM_OR_Monojet_HighPt->fill( Jet2Eta, weight );
          _hist_DeltaEta_PS_VBFDM_OR_Monojet_HighPt->fill( DeltaEta, weight );
          _hist_DeltaPhi_PS_VBFDM_OR_Monojet_HighPt->fill( DeltaPhi, weight );
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

      //      double factor = crossSection()/sumOfWeights();
      double Count_factor = crossSection()/Count_Total;

      // scale(_h_YYYY, crossSection()/sumOfWeights());
      scale(_hist_Count_for_All_PS_CrossSection, Count_factor);

      scale(_hist_Mjj_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_Jet1PT_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_Jet2PT_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_NumJets_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_DeltaEta_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFZ_Baseline, Count_factor);
      scale(_hist_met_PS_VBFZ_Baseline, Count_factor);

      scale(_hist_Mjj_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_Jet1PT_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_Jet2PT_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_NumJets_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_DeltaEta_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFZ_HighMass, Count_factor);
      scale(_hist_met_PS_VBFZ_HighMass, Count_factor);

      scale(_hist_Mjj_PS_VBFZ_Search, Count_factor);
      scale(_hist_Jet1PT_PS_VBFZ_Search, Count_factor);
      scale(_hist_Jet2PT_PS_VBFZ_Search, Count_factor);
      scale(_hist_NumJets_PS_VBFZ_Search, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFZ_Search, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFZ_Search, Count_factor);
      scale(_hist_DeltaEta_PS_VBFZ_Search, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFZ_Search, Count_factor);
      scale(_hist_met_PS_VBFZ_Search, Count_factor);

      scale(_hist_Mjj_PS_VBFDM, Count_factor);
      scale(_hist_Jet1PT_PS_VBFDM, Count_factor);
      scale(_hist_Jet2PT_PS_VBFDM, Count_factor);
      scale(_hist_NumJets_PS_VBFDM, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFDM, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFDM, Count_factor);
      scale(_hist_DeltaEta_PS_VBFDM, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFDM, Count_factor);
      scale(_hist_met_PS_VBFDM, Count_factor);

      scale(_hist_Mjj_PS_VBFDM_100, Count_factor);
      scale(_hist_Jet1PT_PS_VBFDM_100, Count_factor);
      scale(_hist_Jet2PT_PS_VBFDM_100, Count_factor);
      scale(_hist_NumJets_PS_VBFDM_100, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFDM_100, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFDM_100, Count_factor);
      scale(_hist_DeltaEta_PS_VBFDM_100, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFDM_100, Count_factor);
      scale(_hist_met_PS_VBFDM_100, Count_factor);

      scale(_hist_Mjj_PS_Monojet, Count_factor);
      scale(_hist_Jet1PT_PS_Monojet, Count_factor);
      scale(_hist_Jet2PT_PS_Monojet, Count_factor);
      scale(_hist_NumJets_PS_Monojet, Count_factor);
      scale(_hist_Jet1Eta_PS_Monojet, Count_factor);
      scale(_hist_Jet2Eta_PS_Monojet, Count_factor);
      scale(_hist_DeltaEta_PS_Monojet, Count_factor);
      scale(_hist_DeltaPhi_PS_Monojet, Count_factor);
      scale(_hist_met_PS_Monojet, Count_factor);

      scale(_hist_Mjj_PS_Monojet_HighPt, Count_factor);
      scale(_hist_Jet1PT_PS_Monojet_HighPt, Count_factor);
      scale(_hist_Jet2PT_PS_Monojet_HighPt, Count_factor);
      scale(_hist_NumJets_PS_Monojet_HighPt, Count_factor);
      scale(_hist_Jet1Eta_PS_Monojet_HighPt, Count_factor);
      scale(_hist_Jet2Eta_PS_Monojet_HighPt, Count_factor);
      scale(_hist_DeltaEta_PS_Monojet_HighPt, Count_factor);
      scale(_hist_DeltaPhi_PS_Monojet_HighPt, Count_factor);
      scale(_hist_met_PS_Monojet_HighPt, Count_factor);

      scale(_hist_Mjj_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_Jet1PT_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_Jet2PT_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_NumJets_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_DeltaEta_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFDM_OR_Monojet, Count_factor);
      scale(_hist_met_PS_VBFDM_OR_Monojet, Count_factor);

      scale(_hist_Mjj_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_Jet1PT_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_Jet2PT_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_NumJets_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_Jet1Eta_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_Jet2Eta_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_DeltaEta_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_DeltaPhi_PS_VBFDM_OR_Monojet_HighPt, Count_factor);
      scale(_hist_met_PS_VBFDM_OR_Monojet_HighPt, Count_factor);

      /*      
      // normalize(_h_YYYY); // normalize to unity
      normalize(_hist_Count_for_All_PS_CrossSection);
  
      normalize(_hist_Mjj_PS_VBFZ_Baseline);
      normalize(_hist_Jet1PT_PS_VBFZ_Baseline);
      normalize(_hist_Jet2PT_PS_VBFZ_Baseline);
      normalize(_hist_NumJets_PS_VBFZ_Baseline);
      normalize(_hist_Jet1Eta_PS_VBFZ_Baseline);
      normalize(_hist_Jet2Eta_PS_VBFZ_Baseline);
      normalize(_hist_DeltaEta_PS_VBFZ_Baseline);
      normalize(_hist_DeltaPhi_PS_VBFZ_Baseline);
      normalize(_hist_met_PS_VBFZ_Baseline);

      normalize(_hist_Mjj_PS_VBFZ_HighMass);
      normalize(_hist_Jet1PT_PS_VBFZ_HighMass);
      normalize(_hist_Jet2PT_PS_VBFZ_HighMass);
      normalize(_hist_NumJets_PS_VBFZ_HighMass);
      normalize(_hist_Jet1Eta_PS_VBFZ_HighMass);
      normalize(_hist_Jet2Eta_PS_VBFZ_HighMass);
      normalize(_hist_DeltaEta_PS_VBFZ_HighMass);
      normalize(_hist_DeltaPhi_PS_VBFZ_HighMass);
      normalize(_hist_met_PS_VBFZ_HighMass);

      normalize(_hist_Mjj_PS_VBFZ_Search);
      normalize(_hist_Jet1PT_PS_VBFZ_Search);
      normalize(_hist_Jet2PT_PS_VBFZ_Search);
      normalize(_hist_NumJets_PS_VBFZ_Search);
      normalize(_hist_Jet1Eta_PS_VBFZ_Search);
      normalize(_hist_Jet2Eta_PS_VBFZ_Search);
      normalize(_hist_DeltaEta_PS_VBFZ_Search);
      normalize(_hist_DeltaPhi_PS_VBFZ_Search);
      normalize(_hist_met_PS_VBFZ_Search);

      normalize(_hist_Mjj_PS_VBFDM);
      normalize(_hist_Jet1PT_PS_VBFDM);
      normalize(_hist_Jet2PT_PS_VBFDM);
      normalize(_hist_NumJets_PS_VBFDM);
      normalize(_hist_Jet1Eta_PS_VBFDM);
      normalize(_hist_Jet2Eta_PS_VBFDM);
      normalize(_hist_DeltaEta_PS_VBFDM);
      normalize(_hist_DeltaPhi_PS_VBFDM);
      normalize(_hist_met_PS_VBFDM);

      normalize(_hist_Mjj_PS_Monojet);
      normalize(_hist_Jet1PT_PS_Monojet);
      normalize(_hist_Jet2PT_PS_Monojet);
      normalize(_hist_NumJets_PS_Monojet);
      normalize(_hist_Jet1Eta_PS_Monojet);
      normalize(_hist_Jet2Eta_PS_Monojet);
      normalize(_hist_DeltaEta_PS_Monojet);
      normalize(_hist_DeltaPhi_PS_Monojet);
      normalize(_hist_met_PS_Monojet);

      normalize(_hist_Mjj_PS_Monojet_HighPt);
      normalize(_hist_Jet1PT_PS_Monojet_HighPt);
      normalize(_hist_Jet2PT_PS_Monojet_HighPt);
      normalize(_hist_NumJets_PS_Monojet_HighPt);
      normalize(_hist_Jet1Eta_PS_Monojet_HighPt);
      normalize(_hist_Jet2Eta_PS_Monojet_HighPt);
      normalize(_hist_DeltaEta_PS_Monojet_HighPt);
      normalize(_hist_DeltaPhi_PS_Monojet_HighPt);
      normalize(_hist_met_PS_Monojet_HighPt);

      normalize(_hist_Mjj_PS_VBFDM_OR_Monojet);
      normalize(_hist_Jet1PT_PS_VBFDM_OR_Monojet);
      normalize(_hist_Jet2PT_PS_VBFDM_OR_Monojet);
      normalize(_hist_NumJets_PS_VBFDM_OR_Monojet);
      normalize(_hist_Jet1Eta_PS_VBFDM_OR_Monojet);
      normalize(_hist_Jet2Eta_PS_VBFDM_OR_Monojet);
      normalize(_hist_DeltaEta_PS_VBFDM_OR_Monojet);
      normalize(_hist_DeltaPhi_PS_VBFDM_OR_Monojet);
      normalize(_hist_met_PS_VBFDM_OR_Monojet);

      normalize(_hist_Mjj_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_Jet1PT_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_Jet2PT_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_NumJets_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_Jet1Eta_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_Jet2Eta_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_DeltaEta_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_DeltaPhi_PS_VBFDM_OR_Monojet_HighPt);
      normalize(_hist_met_PS_VBFDM_OR_Monojet_HighPt);
      */
    }

  private:
    // Data members like post-cuts event weight counters go here
    /// @name Histograms
    //@{    
    Histo1DPtr _hist_Count_for_All_PS_CrossSection;

    Histo1DPtr _hist_Mjj_PS_VBFZ_Baseline;
    Histo1DPtr _hist_Jet1PT_PS_VBFZ_Baseline;
    Histo1DPtr _hist_Jet2PT_PS_VBFZ_Baseline;
    Histo1DPtr _hist_NumJets_PS_VBFZ_Baseline;
    Histo1DPtr _hist_Jet1Eta_PS_VBFZ_Baseline;
    Histo1DPtr _hist_Jet2Eta_PS_VBFZ_Baseline;
    Histo1DPtr _hist_DeltaEta_PS_VBFZ_Baseline;
    Histo1DPtr _hist_DeltaPhi_PS_VBFZ_Baseline;
    Histo1DPtr _hist_met_PS_VBFZ_Baseline;

    Histo1DPtr _hist_Mjj_PS_VBFZ_HighMass;
    Histo1DPtr _hist_Jet1PT_PS_VBFZ_HighMass;
    Histo1DPtr _hist_Jet2PT_PS_VBFZ_HighMass;
    Histo1DPtr _hist_NumJets_PS_VBFZ_HighMass;
    Histo1DPtr _hist_Jet1Eta_PS_VBFZ_HighMass;
    Histo1DPtr _hist_Jet2Eta_PS_VBFZ_HighMass;
    Histo1DPtr _hist_DeltaEta_PS_VBFZ_HighMass;
    Histo1DPtr _hist_DeltaPhi_PS_VBFZ_HighMass;
    Histo1DPtr _hist_met_PS_VBFZ_HighMass;

    Histo1DPtr _hist_Mjj_PS_VBFZ_Search;
    Histo1DPtr _hist_Jet1PT_PS_VBFZ_Search;
    Histo1DPtr _hist_Jet2PT_PS_VBFZ_Search;
    Histo1DPtr _hist_NumJets_PS_VBFZ_Search;
    Histo1DPtr _hist_Jet1Eta_PS_VBFZ_Search;
    Histo1DPtr _hist_Jet2Eta_PS_VBFZ_Search;
    Histo1DPtr _hist_DeltaEta_PS_VBFZ_Search;
    Histo1DPtr _hist_DeltaPhi_PS_VBFZ_Search;
    Histo1DPtr _hist_met_PS_VBFZ_Search;

    Histo1DPtr _hist_Mjj_PS_VBFDM;
    Histo1DPtr _hist_Jet1PT_PS_VBFDM;
    Histo1DPtr _hist_Jet2PT_PS_VBFDM;
    Histo1DPtr _hist_NumJets_PS_VBFDM;
    Histo1DPtr _hist_Jet1Eta_PS_VBFDM;
    Histo1DPtr _hist_Jet2Eta_PS_VBFDM;
    Histo1DPtr _hist_DeltaEta_PS_VBFDM;
    Histo1DPtr _hist_DeltaPhi_PS_VBFDM;
    Histo1DPtr _hist_met_PS_VBFDM;

    Histo1DPtr _hist_Mjj_PS_VBFDM_100;
    Histo1DPtr _hist_Jet1PT_PS_VBFDM_100;
    Histo1DPtr _hist_Jet2PT_PS_VBFDM_100;
    Histo1DPtr _hist_NumJets_PS_VBFDM_100;
    Histo1DPtr _hist_Jet1Eta_PS_VBFDM_100;
    Histo1DPtr _hist_Jet2Eta_PS_VBFDM_100;
    Histo1DPtr _hist_DeltaEta_PS_VBFDM_100;
    Histo1DPtr _hist_DeltaPhi_PS_VBFDM_100;
    Histo1DPtr _hist_met_PS_VBFDM_100;

    Histo1DPtr _hist_Mjj_PS_Monojet;
    Histo1DPtr _hist_Jet1PT_PS_Monojet;
    Histo1DPtr _hist_Jet2PT_PS_Monojet;
    Histo1DPtr _hist_NumJets_PS_Monojet;
    Histo1DPtr _hist_Jet1Eta_PS_Monojet;
    Histo1DPtr _hist_Jet2Eta_PS_Monojet;
    Histo1DPtr _hist_DeltaEta_PS_Monojet;
    Histo1DPtr _hist_DeltaPhi_PS_Monojet;
    Histo1DPtr _hist_met_PS_Monojet;

    Histo1DPtr _hist_Mjj_PS_Monojet_HighPt;
    Histo1DPtr _hist_Jet1PT_PS_Monojet_HighPt;
    Histo1DPtr _hist_Jet2PT_PS_Monojet_HighPt;
    Histo1DPtr _hist_NumJets_PS_Monojet_HighPt;
    Histo1DPtr _hist_Jet1Eta_PS_Monojet_HighPt;
    Histo1DPtr _hist_Jet2Eta_PS_Monojet_HighPt;
    Histo1DPtr _hist_DeltaEta_PS_Monojet_HighPt;
    Histo1DPtr _hist_DeltaPhi_PS_Monojet_HighPt;
    Histo1DPtr _hist_met_PS_Monojet_HighPt;

    Histo1DPtr _hist_Mjj_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_Jet1PT_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_Jet2PT_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_NumJets_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_Jet1Eta_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_Jet2Eta_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_DeltaEta_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_DeltaPhi_PS_VBFDM_OR_Monojet;
    Histo1DPtr _hist_met_PS_VBFDM_OR_Monojet;

    Histo1DPtr _hist_Mjj_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_Jet1PT_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_Jet2PT_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_NumJets_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_Jet1Eta_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_Jet2Eta_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_DeltaEta_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_DeltaPhi_PS_VBFDM_OR_Monojet_HighPt;
    Histo1DPtr _hist_met_PS_VBFDM_OR_Monojet_HighPt;

    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_VBFDM_Absolute);
 

}
