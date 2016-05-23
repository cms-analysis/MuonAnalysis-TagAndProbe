#include <iostream>
#include <vector>
#include <sstream>
#include "plotUtil_Tracking.cxx"
#include "TGraphAsymmErrors.h"

void plotSingleDelta(TString plotName, TString binName, TString xAxisName, TString fake = "");

void plotCorrection(TString earlyTag = ""){

  gROOT->SetBatch();

  gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
  gStyle->SetOptStat(0);
  //gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php plots/comparisonMatches/");

  datalbl = "Raw efficiency";
  reflbl  = "Corrected efficiency";
  cmslabel = "CMS Preliminary";
  cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";

  TString match = earlyTag + "dr030e030";
  std::vector<TString> binning{"eff_aeta"};
//abseta", "eff_eta", "eff_eta2", "eff_eta3", "eff", "eff_vtx", "eff_phi"};
  
  TString fileName;
  if (match.Contains("tk0") )  fileName = "plots/Run2015D_16Dec2015-v1_IsoMu20_vs_DY76_IsoMu20/" + match + "/fits.root";
  else fileName = "plots/Run2015D_16Dec2015-v1_IsoMu20_vs_DY76_IsoMu20/fits.root";
  std::cout << fileName << std::endl;
  TFile *file = new TFile(fileName);

  TString nameEffData = "fit_" + binning.at(0) + "_" + match ;
  TString nameEffCorrData = binning.at(0) + "_" + match + "_corr";
  //TString nameEffMC = "ref_" + binning.at(0) + "_" + match ;
  std::cout << nameEffData << std::endl;
  std::cout << nameEffCorrData << std::endl;
  TGraphAsymmErrors *graphEff  = (TGraphAsymmErrors*)file->Get(nameEffData);
  TGraphAsymmErrors *graphEffCorr  = (TGraphAsymmErrors*)file->Get(nameEffCorrData);
  if (graphEff==0 || graphEffCorr==0)  {
    std::cerr << graphEff << " not found. " << std::endl;
    return;
  }
  TCanvas *c1 = new TCanvas();
  c1->cd();
  //TLegend *leg;
  //leg = new TLegend(.42, .75, .72, .93 );
  
  double yMin = 0.0;
  if (match.Contains("tk0") ) yMin = 0.958;
  else  yMin = 0.978;
  double yMax = 1.002;

  if (match.Contains("tk0") ){
    graphEff->SetLineColor(kViolet-2);
    graphEff->SetMarkerColor(kViolet-2);
  } else {
    graphEff->SetLineColor(65);
    graphEff->SetMarkerColor(65);
  } 				
  graphEffCorr->SetLineColor(kBlack);
  graphEffCorr->SetMarkerColor(kBlack);
  graphEff->SetLineWidth(2);
  graphEff->SetMarkerSize(1.6);
  graphEffCorr->SetLineWidth(2);
  graphEffCorr->SetMarkerSize(1.6);

  //leg->AddEntry(graphEff, "Raw efficiancy", "LP");
  //leg->AddEntry(graphEffCorr, "Corrected efficiency", "LP");

  graphEff->GetYaxis()->SetTitle("Efficiency");
  graphEff->GetYaxis()->SetRangeUser(yMin, yMax);
  if (match.Contains("tk0") )  graphEff->GetYaxis()->SetNdivisions(505);
  else  graphEff->GetYaxis()->SetNdivisions(503);
  graphEff->GetYaxis()->SetDecimals(true);
  graphEff->GetXaxis()->SetTitle("muon |#eta|");

  graphEff->Draw("AP");
  graphEffCorr->Draw("PSAME");
  if (datalbl) doLegend(graphEff,graphEffCorr,datalbl, reflbl);
  //leg->Draw();
  c1->Print("effcomparison.png", ".png");
  c1->Print("effcomparison.pdf", ".pdf");

}
/*
  std::vector<TString> dataOrMC{"fit_", "ref_", "ratio_"};
  std::vector<TString> binning{"eff_abseta", "eff_eta", "eff_eta2", "eff_eta3", "eff", "eff_vtx", "eff_phi"};
  std::vector<TString> xvars{ "muon |#eta|", "muon #eta","muon #eta","muon #eta", "muon #eta", "N(primary vertices)", "muon #phi"};

  for(unsigned int n = 0; n < dataOrMC.size() ; ++n){
    for(unsigned int bin = 0; bin < binning.size() ; ++bin){
    
      TString plotName = dataOrMC.at(n) + binning.at(bin) + "_";
      //std::cout << plotName << std::endl;
      TString binName = "";
      plotSingleDelta(plotName, binName, xvars.at(bin));

      if(dataOrMC.at(n) != "ratio_"){
        plotSingleDelta(plotName, binName, xvars.at(bin), "_fake1");
        plotSingleDelta(plotName, binName, xvars.at(bin), "_fake2");
      }

    }
  }

}

void plotSingleDelta(TString plotName, TString binName, TString xAxisName, TString fake = ""){

  std::vector<TString> matches;
  matches.push_back("dr010e010");
  matches.push_back("dr030e015");
  matches.push_back("dr050e020");
  matches.push_back("dr070e030");
  matches.push_back("dr100e040");

  std::vector<std::string> matchesLeg;
  matchesLeg.push_back("#DeltaR, #Delta#eta = 0.10, 0.10");
  matchesLeg.push_back("#DeltaR, #Delta#eta = 0.30, 0.15");
  matchesLeg.push_back("#DeltaR, #Delta#eta = 0.50, 0.20");
  matchesLeg.push_back("#DeltaR, #Delta#eta = 0.70, 0.30");
  matchesLeg.push_back("#DeltaR, #Delta#eta = 1.00, 0.40");

  std::vector<TColor*> colors;
//  TColor *col26 = gROOT->GetColor(26);
  colors.push_back(gROOT->GetColor(16));
  colors.push_back(gROOT->GetColor(46));
  colors.push_back(gROOT->GetColor(66));
  colors.push_back(gROOT->GetColor(22));
  colors.push_back(gROOT->GetColor(30));

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg;
  if (!fake.Contains("fake")) leg = new TLegend(.42, .15, .72, .35 );
  if (fake.Contains("fake")) leg = new TLegend(.42, .75, .72, .93 );

  for(unsigned int n = 0; n < matches.size() ; ++n){

    TString fileName = "plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + matches.at(n) + "/fits.root";
    //std::cout << fileName << std::endl;
    TFile *file = new TFile(fileName);

    TString graphName = plotName + matches.at(n) ;
    if (fake != 0) graphName += fake;
    if (graphName.Contains("ratio")) graphName += "_corr";
    //std::cout << graphName << std::endl;
    TGraphAsymmErrors *g1  = (TGraphAsymmErrors*)file->Get(graphName);
    if (g1==0)  {
      std::cerr << graphName << " not found. " << std::endl;
      return;
    }

    g1->SetLineColor(kRed+n);
    g1->SetLineWidth(2);
    g1->SetMarkerColor(kRed+n);
    g1->SetMarkerSize(1.6);
    leg->AddEntry(g1, matchesLeg.at(n).c_str(), "LP");

    if( n == 0 ){
      double yMin = 0.8;
      double yMax = 1.009;
      g1->GetXaxis()->SetTitle(xAxisName);
      g1->GetYaxis()->SetTitleOffset(1.2);

      if (graphName.Contains("fit_")) {
        g1->GetYaxis()->SetTitle("Data efficiency");
        if (graphName.Contains("_fake")) {
          yMin = 0.0; yMax = 1.2;
          g1->GetYaxis()->SetTitle("Data Fake match rate");
        }
        g1->GetYaxis()->SetRangeUser(yMin, yMax);
      }
      if (graphName.Contains("ref_")) { 
        g1->GetYaxis()->SetTitle("MC efficiency");
        if (graphName.Contains("_fake")) {
          yMin = 0.0; yMax = 1.2;
          g1->GetYaxis()->SetTitle("MC Fake match rate");
        }
        g1->GetYaxis()->SetRangeUser(yMin, yMax);
      }

      if (graphName.Contains("ratio_")) {
        yMin = 0.6; yMax = 1.4;
        g1->GetYaxis()->SetTitle("Data / Sim. ratio");
        g1->GetYaxis()->SetRangeUser(yMin, yMax);
      }
      g1->Draw("AP");

    }
    if( n > 0 ) g1->Draw("Psame");


  }

  leg->Draw();
  if (fake != 0) plotName += fake;
  c1->Print("plots/comparisonMatches/"+plotName+binName+".png", ".png");

  return;
}
*/
