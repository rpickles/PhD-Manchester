#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include "stylefile.h"
#include <TStyle.h>
#include <TROOT.h>
#include "stdio.h"
#include <TSystem.h>
#include "TLegend.h"
#include "TColor.h"

using namespace std;

Double_t Lumi = 3.;
Double_t CrossSectionDM = 3.;
Double_t CrossSectionRatio = 3.;
Double_t LumiUncertainty = 2.;
Double_t LumiRatioUncertainty = 2.;

// Custom palette for observable phase space plot
void set_plot_style() {
  const int NRGBs = 5, NCont = 80;
  Double_t stops[NRGBs] = { 0.00, 0.3, 0.6, 0.99999, 1.0};
  Double_t red[NRGBs]   = {0.15, 0.20, 0.60, 1.0, 1.0};
  Double_t green[NRGBs] = {0.30, 0.62, 0.93, 1.0, 1.0};
  Double_t blue[NRGBs]  = {0.42, 0.78, 0.98, 1.0, 1.0};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}


// Statistical test with number of options 
Double_t statTest(TH1F* modelRatio, TH1F* dataRatio, TString statOption, TString BinErrorString, Double_t CrossSectionDM, Double_t CrossSectionRatio) {

  if ( statOption == "chi2" ) {
    if ( BinErrorString == "standard" || BinErrorString == "Standard" ){
      return dataRatio->Chi2Test(modelRatio,"WW_P",0);
    }
    else if ( BinErrorString.Contains("Lumi") ){
      if ( BinErrorString.Contains("2016") ){
	Lumi = 30.;
      } 
      else if ( BinErrorString.Contains("2018") ){
	Lumi = 300.;
      }
      else if ( BinErrorString.Contains("2023") ){
	Lumi = 3000.;
      }
      else {
	Lumi = 3.;
      }

      LumiUncertainty = sqrt(Lumi*CrossSectionDM);
      LumiRatioUncertainty = sqrt(Lumi*CrossSectionRatio);

      int nBins = modelRatio->GetXaxis()->GetNbins();
      for ( int i = 0 ; i < nBins ; i++ ){
        modelRatio->SetBinError(i+1,LumiUncertainty);
	dataRatio->SetBinError(i+1,LumiRatioUncertainty);
      }
      return dataRatio->Chi2Test(modelRatio,"WW_P",0);
    }
    else{
      double binError = atof(BinErrorString);
      int nBins = modelRatio->GetXaxis()->GetNbins();
      for ( int i = 0 ; i < nBins ; i++ ){
	modelRatio->SetBinError(i+1,(modelRatio->GetBinContent(i+1))*binError);
      }
      return dataRatio->Chi2Test(modelRatio,"WW_P",0);
    }
  }
  else if ( statOption == "KS" ) {
    return 2;
  }
  else { return 3; }

};

int statisticsTest(TString variable="", TString PS="", TString Dimension="", TString StatOption="", TString BinErrorString=""){

  gStyle->SetTitleBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetHistLineWidth(2);
  gStyle->SetTextFont(42);
  //  gStyle->SetPalette(1);
  set_plot_style();

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  //  TString BinErrorString = Form("BinError%.2f",binError);

  gSystem->Exec("mkdir -p Figures/StatPlots/"+Dimension+"/Absolute/"+BinErrorString+"");

  TCanvas *statCanv;
  statCanv = new TCanvas("statCanv","",3000,1500);

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFuncColor(kRed);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetTextAlign(22);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);
  statCanv->cd();
    
  /*  // Reading in root file containing data histograms                                                                                                                           
  TFile* dataFile = new TFile("~/Documents/Rivet_Analyses/MC_VBFDM/PlotCombinationTool/corrected.root", "READ");
  dataFile->cd();
  TH1F* dataCurve;

  // Checking which variable to compare
  if ( variable == "Mjj" ){
    dataCurve = (TH1F*)gDirectory->Get("Mjj_search");
  }
  else if ( variable == "Etmiss" ){
    dataCurve = (TH1F*)gDirectory->Get("MET_mono");
  }
  else { dataCurve = (TH1F*)gDirectory->Get("MET_mono");}

  cout<<"Number of data bins = "<<dataCurve->GetXaxis()->GetNbins()<<endl;
  */
  // Scale of number of leptons for ratio
  //double ratioScale(6.03);
  int EFTScaleSize(10);
  int DMmassSize(4);

  // Defining low limit on the EFT Scale for each Dimension.
  Double_t EFTScaleMin;
  
  if ( Dimension == "D5c" ) { EFTScaleMin = 3300.; }
  else if ( Dimension == "D5d" ) { EFTScaleMin = 6600.; }
  else if ( Dimension == "D6a" ) { EFTScaleMin = 230.; }
  else if ( Dimension == "D6b" ) { EFTScaleMin = 330.; }
  else { EFTScaleMin = 100. ; }

  Double_t EFTScale[10];
  TString EFTScaleString[10];
  //  EFTScale[0] = EFTScaleMin;
  //EFTScaleString[0] = "#Lambda_Min";

  for ( int i = 0 ; i < EFTScaleSize ; ++i ){
    Double_t j = 1+(i*0.25);
    EFTScale[i] = j*EFTScaleMin;
    EFTScaleString[i] = Form("%.2f #Lambda_Min", j);
  }

  // Defining different mass options
  Double_t DMmass[4] = { 1.0, 10.0, 100.0, 1000.0 }; //0.1

  TString DMmassString[4];
  //  DMmassString[0] = Form("%.1f",DMmass[0]);
  for ( int i = 0 ; i < DMmassSize ; ++i ){
    DMmassString[i] = Form("%.0f",DMmass[i]);
  }

  // Initialising statistics result 2d histogram
  TH2F* statHist = new TH2F("statHist", "statHist", DMmassSize, DMmass[0]-(DMmass[0]/2), DMmass[2]+(DMmass[2]/2), EFTScaleSize, EFTScale[0]-(EFTScale[0]/2), EFTScale[1]+(EFTScale[1]/2));

  TH1F* ratioHist[DMmassSize][EFTScaleSize];
  TH1F* SMRatioHist;
  TH1F* xsec_RatioCurve;
  Double_t statTestResult[DMmassSize][EFTScaleSize];
  double EFTScaleForRatio[DMmassSize][EFTScaleSize];

  for ( int i = 0 ; i < DMmassSize ; ++i ){
    for ( int j = 0 ; j < EFTScaleSize ; ++j ){
    // Reading in root file containing DM model and background histograms
    TFile* modelFile = new TFile("~/Documents/Rivet_Analyses/MC_VBFDM/PlotCombinationTool/Figures/"+Dimension+"/Absolute/"+Dimension+"_"+variable+"_PS_"+PS+".root", "READ");
    modelFile->cd();
    TH1F* DMcurve = (TH1F*)gDirectory->Get("Mass"+DMmassString[i]+"");
    TH1F* EWK_and_fifthQCD = (TH1F*)gDirectory->Get("EWK_and_fifthQCD");
    TH1F* cRatio = (TH1F*)gDirectory->Get("cRatio");
    TH1F* zmumu = (TH1F*)gDirectory->Get("zmumu");

    // Reading in root file containing Cross section histogram
    TFile* crossSectionFile = new TFile("~/Documents/Rivet_Analyses/MC_VBFDM/PlotCombinationTool/Figures/Mass"+DMmassString[i]+"/Absolute/Mass"+DMmassString[i]+"_Count_for_All_PS_CrossSection.root", "READ");
    crossSectionFile->cd();
    TH1F* crossSectionDMCurve = (TH1F*)gDirectory->Get(""+Dimension+"");
    TH1F* crossSectionEWK_and_fifthQCD = (TH1F*)gDirectory->Get("EWK_and_fifthQCD");
    TH1F* crossSectioncRatio = (TH1F*)gDirectory->Get("cRatio");
    TH1F* crossSectionzmumu = (TH1F*)gDirectory->Get("zmumu");

    cout<<"Number of DM bins = "<<DMcurve->GetXaxis()->GetNbins()<<endl;

    // Defining the scaling power for each dimension                                                                                                                             
    double power = 0.;
    if ( Dimension == "D5a" || Dimension == "D5b" || Dimension == "D5c" || Dimension == "D5d" ) { power = 4; }
    if ( Dimension == "D6a" || Dimension == "D6b" ) { power = 5; }
    if ( Dimension == "D7a" || Dimension == "D7b" || Dimension == "D7c" || Dimension == "D7d" ) { power = 6; }

    EFTScaleForRatio[i][j] = pow(((EFTScaleMin)/(EFTScale[j])),power);

    DMcurve->Scale(EFTScaleForRatio[i][j]);
    crossSectionDMCurve->Scale(EFTScaleForRatio[i][j]);

    // Producing the ratio of DM+Backgrounds/Backgrounds cross-sections and scaling due to number of leptons
    TH1F xsec_up_NotP = (*crossSectionEWK_and_fifthQCD)+(*crossSectionDMCurve);
    TH1F* xsec_EWK_and_fifthQCD_and_DMcurve = &xsec_up_NotP;
    TH1F* xsec_modelRatio = (TH1F*)xsec_EWK_and_fifthQCD_and_DMcurve->Clone("xsec_modelRatio");
    xsec_modelRatio->Sumw2();
    xsec_modelRatio->Divide(crossSectionzmumu);
    
    // Producing the ratio of DM+Backgrounds/Backgrounds and scaling due to number of leptons 
    TH1F EWK_and_fifthQCD_and_DMcurve_NotP = (*EWK_and_fifthQCD)+(*DMcurve);
    TH1F* EWK_and_fifthQCD_and_DMcurve = &EWK_and_fifthQCD_and_DMcurve_NotP;
    TH1F* modelRatio = (TH1F*)EWK_and_fifthQCD_and_DMcurve->Clone("modelRatio");
    modelRatio->Sumw2();
    modelRatio->Divide(zmumu);

    ratioHist[i][j] = (TH1F*)modelRatio->Clone("ratioHist");
    //    ratioHist[i][j]->GetXaxis()->SetNdivisions(5);
    SMRatioHist = (TH1F*)cRatio->Clone("SMRatioHist");
    xsec_RatioCurve = (TH1F*)crossSectioncRatio->Clone("xsec_RatioCurve");

    CrossSectionDM = xsec_modelRatio->GetBinContent(3);
    CrossSectionRatio = xsec_RatioCurve->GetBinContent(3);

    // Running the histograms through a statistical test and outputting the p-value
    statTestResult[i][j] = statTest(ratioHist[i][j], SMRatioHist, StatOption, BinErrorString, CrossSectionDM, CrossSectionRatio);
    if ( statTestResult[i][j] < 1e-43 ){ statTestResult[i][j] = 1e-45; }
    cout<<"EFTScaleForRatio = "<<EFTScaleForRatio[i][j]<<endl;
    statHist->Fill(DMmassString[i], EFTScaleString[j], statTestResult[i][j]);
  
    }
  }

  // Writing output to a root file.
  TFile * outfile = new TFile("Figures/StatPlots/"+Dimension+"/Absolute/"+BinErrorString+"/Stats_"+Dimension+"_"+variable+"_"+BinErrorString+"_PS_"+PS+".root", "RECREATE");
  outfile->cd();

  for ( int i = 0 ; i < DMmassSize ; ++i ){
    for ( int j = 0 ; j < EFTScaleSize ; ++j ){
      ratioHist[i][j]->Write();
     
    }
  }
  //  dataCurve->Write();
  SMRatioHist->Write();
  statHist->Write();
  outfile->Close();
  
  // Two pads are manipulated for appearence

  TPad* mainPad = new TPad("mainPad","", 0.52, 0.15, 0.99, 0.99); mainPad->Draw();
  TPad* ratioPad = new TPad("ratioPad","", 0.01, 0.15, 0.48, 0.99); ratioPad->Draw();
  TLegend* ratioLegend = new TLegend(0.07,0.04,0.99,0.15); ratioLegend->Draw();

  mainPad->cd(); statHist->GetZaxis()->SetRangeUser(0, 1); statHist->Draw("COLZ"); statHist->Draw("Text same"); statHist->GetXaxis()->SetTitle("DM mass (GeV)"); /*statHist->GetYaxis()->SetTitle("EFT Scale");*/ statHist->GetYaxis()->SetTitleOffset(1.5); gPad->SetRightMargin(1); statHist->GetZaxis()->SetLabelOffset(0); gPad->Update();// gPad->SetLogz(1);

  ratioPad->cd(); ratioHist[0][0]->GetYaxis()->SetTitle("#frac{#sigma((Z#rightarrow #nu #bar{#nu})jj) + #sigma((Z#rightarrow #chi #bar{#chi})jj)}{#sigma((Z#rightarrow #mu+ #mu-)jj)}"); ratioHist[0][0]->GetYaxis()->SetTitleSize(50); ratioHist[0][0]->GetYaxis()->SetTitleFont(43); ratioHist[0][0]->GetYaxis()->SetTitleOffset(2); ratioHist[0][0]->GetXaxis()->SetTitle(""+variable+""); ratioHist[0][0]->GetXaxis()->SetTitleSize(50); ratioHist[0][0]->GetXaxis()->SetTitleFont(43); ratioHist[0][0]->GetXaxis()->SetTitleOffset(2); ratioHist[0][0]->SetTickLength(0.08); ratioHist[0][0]->SetMinimum(0); ratioHist[0][0]->SetMaximum(25); 

  //   gPad->SetLogy(1);
  
  int k(2);
  int n(0);
  for ( int i = 0 ; i < 4 ; ++i ){
    int m(0);
    for ( int j = 0 ; j < 6 ; ++j ){
      if ( m < EFTScaleSize && n < DMmassSize ){ 
	ratioHist[n][m]->SetMarkerColor(k);
	ratioHist[n][m]->SetLineColor(k);
	ratioHist[n][m]->Draw("ep same");
	m = m+2;
	n = n+1;
	k++;
      }
    }
  }
  /*
  dataCurve->SetMarkerStyle(1);
  dataCurve->SetMarkerColor(1);
  dataCurve->SetLineColor(1);
  dataCurve->Draw("ep same");
  */
  SMRatioHist->SetMarkerStyle(1);
  SMRatioHist->SetMarkerColor(1);
  SMRatioHist->SetLineColor(1);
  SMRatioHist->Draw("ep same");

  ratioLegend->SetNColumns(3);
  ratioLegend->SetFillColor(0);
  ratioLegend->SetTextSize(0.02);
  ratioLegend->SetBorderSize(0);
  ratioLegend->SetTextFont(42);
  //  ratioLegend->AddEntry(dataCurve,"data","PL");
  ratioLegend->AddEntry(SMRatioHist, "#frac{#sigma ( (Z#rightarrow #nu #bar{#nu}) jj )}{#sigma ( (Z#rightarrow #mu+ #mu-) jj )}","PL");

  int p(0);
  for ( int i = 0 ; i < 5 ; ++i ){
    int m(0);
    for ( int j = 0 ; j < 6 ; ++j ){
      if ( m < EFTScaleSize && p < DMmassSize ){
	ratioLegend->AddEntry(ratioHist[p][m], Form(""+Dimension+", DM mass = %.1f GeV, #Lambda = %.0f GeV",DMmass[p],EFTScale[m]), "PL");
	m = m+2;
	p = p+1;
      }
  }
}

  statCanv->cd();

  // THIS ACTUALLY WRITES OUT THE PLOTS (PDF/PNG/EPS FORMAT)
  gPad->RedrawAxis();statCanv->Print("Figures/StatPlots/"+Dimension+"/Absolute/"+BinErrorString+"/Stats_"+Dimension+"_"+variable+"_"+BinErrorString+"_PS_"+PS+".pdf");
  gPad->RedrawAxis();statCanv->Print("Figures/StatPlots/"+Dimension+"/Absolute/"+BinErrorString+"/Stats_"+Dimension+"_"+variable+"_"+BinErrorString+"_PS_"+PS+".png");
  gPad->RedrawAxis();statCanv->Print("Figures/StatPlots/"+Dimension+"/Absolute/"+BinErrorString+"/Stats_"+Dimension+"_"+variable+"_"+BinErrorString+"_PS_"+PS+".eps");

  return 0;

}


