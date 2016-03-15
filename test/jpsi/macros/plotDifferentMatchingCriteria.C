#include <iostream>
#include <vector>
#include "../plotUtil_Tracking.cxx"

void plotSingleDelta(TString plotName, TString binName, TString xAxisName, TString fake = "");

void plotDifferentMatchingCriteria(){
  cmslabel = "CMS Preliminary";
  cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";

  gROOT->SetBatch();

  gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
  gStyle->SetOptStat(0);
//  gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php plots/comparisonMatches/");

  std::vector<TString> dataOrMC{"fit_", "ref_", "ratio_", ""};
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

  std::vector<int> colors;
  colors.push_back(99);
  colors.push_back(95);
  colors.push_back(78);
  colors.push_back(66);
  colors.push_back(57);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg;
  if (!fake.Contains("fake")) leg = new TLegend(.30, .20, .75, .45 );
  if (fake.Contains("fake")) leg = new TLegend(.30, .65, .75, .90 );
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);

  for(unsigned int n = 0; n < matches.size() ; ++n){
    TString fileName = "../plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + matches.at(n) + "/fits.root";
    //std::cout << fileName << std::endl;
    TFile *file = new TFile(fileName);

    TString graphName = plotName + matches.at(n) ;
    if (fake != 0) graphName += fake;
    if (graphName.Contains("ratio")) graphName += "_corr";
    if (!graphName.Contains("_fake") && !graphName.Contains("fit_") && !graphName.Contains("ref_") && !graphName.Contains("ratio")) graphName += "_corr";
    //std::cout << graphName << std::endl;
    TGraphAsymmErrors *g1  = (TGraphAsymmErrors*)file->Get(graphName);
    if (g1==0)  {
      std::cerr << graphName << " not found. " << std::endl;
      return;
    }

    g1->SetLineColor(colors.at(n));
    g1->SetLineWidth(2);
    g1->SetMarkerColor(colors.at(n));
    g1->SetMarkerSize(1.6);
    leg->AddEntry(g1, matchesLeg.at(n).c_str(), "LP");

    if( n == 0 ){
      double yMin = 0.90;
      double yMax = 1.01;//009;
      g1->GetXaxis()->SetTitle(xAxisName);
      g1->GetYaxis()->SetTitleOffset(1.2);
      g1->GetYaxis()->SetTitle("Data efficiency");
      g1->GetYaxis()->SetRangeUser(yMin, yMax);

      if (graphName.Contains("fit_")) {
        g1->GetYaxis()->SetTitle("Raw data efficiency");
        if (graphName.Contains("_fake")) {
          yMin = 0.0; yMax = 1.5;
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
        yMin = 0.9; yMax = 1.1;
        g1->GetYaxis()->SetTitle("Data / Sim. ratio");
        g1->GetYaxis()->SetRangeUser(yMin, yMax);
      }
      g1->Draw("AP");

    }
    if( n > 0 ) g1->Draw("Psame");


  }

  //printcmsprelim();
  //printcmsconditionslabel();
  leg->Draw();
  if (fake != 0) plotName += fake;
  c1->Print(plotName+binName+".pdf", ".pdf");
  c1->Print(plotName+binName+".png", ".png");

  return;
}
