#include <iostream>
#include <vector>
#include <sstream>
#include "../plotUtil_Tracking.cxx"
#include "TGraphAsymmErrors.h"

void plotPt(){

  gROOT->SetBatch();

  gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
  gStyle->SetOptStat(0);
  //gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php plots/comparisonMatches/");

  datalbl = "Raw efficiency";
  reflbl  = "Corrected efficiency";
  cmslabel = "CMS Preliminary";
  cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";

  TString match = "dr030e015";
  TString tk0match = "tk0_dr030e015";
  std::vector<TString> binning{"eff_pt"};
//abseta", "eff_eta", "eff_eta2", "eff_eta3", "eff", "eff_vtx", "eff_phi"};
  
  TString fileNameGeneral;
  TString fileNameEarly;
  fileNameGeneral = "../plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + match + "/fits_pt.root";
  fileNameEarly =   "../plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + match + "/tk0/fits_pt.root";
  std::cout << fileNameGeneral << std::endl;
  std::cout << fileNameEarly   << std::endl;
  TFile *fileGeneral = new TFile(fileNameGeneral);
  TFile *fileEarly = new TFile(fileNameEarly);

  TString nameEffDataCorrGeneral = binning.at(0) + "_" + match + "_corr" ;
  TString nameEffMCCorrGeneral = "ref_" + binning.at(0) + "_" + match + "_corr";
  TString nameEffDataCorrEarly = binning.at(0) + "_" + tk0match + "_corr" ;
  TString nameEffMCCorrEarly = "ref_" + binning.at(0) + "_" + tk0match + "_corr";
  std::cout << nameEffDataCorrGeneral << std::endl;
  std::cout << nameEffMCCorrEarly << std::endl;
  TGraphAsymmErrors *graphEffCorrDataGeneral = (TGraphAsymmErrors*)fileGeneral->Get(nameEffDataCorrGeneral);
  TGraphAsymmErrors *graphEffCorrMCGeneral = (TGraphAsymmErrors*)fileGeneral->Get(nameEffMCCorrGeneral);
  TGraphAsymmErrors *graphEffCorrDataEarly = (TGraphAsymmErrors*)fileEarly->Get(nameEffDataCorrEarly);
  TGraphAsymmErrors *graphEffCorrMCEarly = (TGraphAsymmErrors*)fileEarly->Get(nameEffMCCorrEarly);
  if (graphEffCorrDataGeneral==0 || graphEffCorrMCGeneral==0 || graphEffCorrDataEarly==0 || graphEffCorrMCEarly==0 )  {
    std::cerr << "At least one graph not found. " << std::endl;
    return;
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  double yMin = 0.96;
  double yMax = 1.002;

  graphEffCorrDataGeneral->SetLineColor(kBlack);
  graphEffCorrDataGeneral->SetMarkerColor(kBlack);
  graphEffCorrMCGeneral->SetLineColor(65);
  graphEffCorrMCGeneral->SetMarkerColor(65);
  graphEffCorrMCGeneral->SetFillColorAlpha(65, 0.3);

  graphEffCorrDataEarly->SetLineColor(kViolet-2);
  graphEffCorrDataEarly->SetMarkerColor(kViolet-2);
  graphEffCorrMCEarly->SetLineColor(kViolet-2);
  graphEffCorrMCEarly->SetMarkerColor(kViolet-2);
  graphEffCorrMCEarly->SetFillColorAlpha(kViolet-2,0.5);

  graphEffCorrDataGeneral->SetLineWidth(2);
  graphEffCorrDataGeneral->SetMarkerSize(1.6);
  graphEffCorrMCGeneral->SetLineWidth(2);
  graphEffCorrMCGeneral->SetMarkerSize(1.6);

  graphEffCorrDataEarly->SetLineWidth(2);
  graphEffCorrDataEarly->SetMarkerSize(1.6);
  graphEffCorrMCEarly->SetLineWidth(2);
  graphEffCorrMCEarly->SetMarkerSize(1.6);

  graphEffCorrMCGeneral->GetYaxis()->SetTitle("Efficiency");
  graphEffCorrMCGeneral->GetYaxis()->SetRangeUser(yMin, yMax);
  graphEffCorrMCGeneral->GetYaxis()->SetNdivisions(505);
  graphEffCorrMCGeneral->GetYaxis()->SetDecimals(true);
  graphEffCorrMCGeneral->GetXaxis()->SetTitle("p_{T}");

  graphEffCorrMCGeneral->Draw("A2");
  graphEffCorrMCEarly->Draw("2SAME");
  graphEffCorrDataEarly->Draw("PSAME");
  graphEffCorrDataGeneral->Draw("PSAME");

  TLegend *leg;
  leg = new TLegend(.20, .20, .90, .40 );
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->AddEntry(graphEffCorrDataGeneral, "Data - generalTracks", "LP");
  leg->AddEntry(graphEffCorrMCGeneral, "Simulation - generalTracks", "F");
  leg->AddEntry(graphEffCorrDataEarly, "Data - earlyGeneralTracks", "LP");
  leg->AddEntry(graphEffCorrMCEarly, "Simulation - earlyGeneralTracks", "F");
  leg->Draw("same");

  TPaveText *cmsprel = new TPaveText(0.46 ,.96,.94,.99,"NDC");
  cmsprel->SetTextSize(0.05);
  cmsprel->SetFillColor(0);
  cmsprel->SetFillStyle(0);
  cmsprel->SetLineStyle(2);
  cmsprel->SetLineColor(0);
  cmsprel->SetTextAlign(12);
  cmsprel->SetTextFont(42);
  cmsprel->AddText(cmsconditionslabel);
  cmsprel->Draw("same");

    TPaveText *cmsprel2 = new TPaveText(0.20 ,.41,.27,.45,"NDC");
    cmsprel2->SetTextSize(0.04);
    cmsprel2->SetFillColor(0);
    cmsprel2->SetFillStyle(0);
    cmsprel2->SetLineStyle(2);
    cmsprel2->SetLineColor(0);
    cmsprel2->SetTextAlign(12);
    cmsprel2->SetTextFont(42);
    cmsprel2->AddText(cmslabel);
    cmsprel2->Draw("same");


  c1->Print("eff_pt_comparison.png", ".png");
  c1->Print("eff_pt_comparison.pdf", ".pdf");

}



void plotPt_ratio(){

  gROOT->SetBatch();

  gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
  gStyle->SetOptStat(0);
  //gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php plots/comparisonMatches/");

  datalbl = "Raw efficiency";
  reflbl  = "Corrected efficiency";
  cmslabel = "CMS Preliminary";
  cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";

  TString match = "dr030e015";
  TString tk0match = "tk0_dr030e015";
  std::vector<TString> binning{"eff_pt"};
//abseta", "eff_eta", "eff_eta2", "eff_eta3", "eff", "eff_vtx", "eff_phi"};
  
  TString fileNameGeneral;
  TString fileNameEarly;
  fileNameGeneral = "../plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + match + "/fits_pt.root";
  fileNameEarly =   "../plots/Run2015D_16Dec2015-v1_Json_vs_JPsiPt8_" + match + "/tk0/fits_pt.root";
  std::cout << fileNameGeneral << std::endl;
  std::cout << fileNameEarly   << std::endl;
  TFile *fileGeneral = new TFile(fileNameGeneral);
  TFile *fileEarly = new TFile(fileNameEarly);

  TString nameEffRatioCorrGeneral = "ratio_" + binning.at(0) + "_" + match + "_corr" ;
  TString nameEffRatioCorrEarly = "ratio_" + binning.at(0) + "_" + tk0match + "_corr" ;

  std::cout << nameEffRatioCorrGeneral << std::endl;
  std::cout << nameEffRatioCorrEarly << std::endl;
  TGraphAsymmErrors *graphEffCorrRatioGeneral = (TGraphAsymmErrors*)fileGeneral->Get(nameEffRatioCorrGeneral);
  TGraphAsymmErrors *graphEffCorrRatioEarly = (TGraphAsymmErrors*)fileEarly->Get(nameEffRatioCorrEarly);
  if (graphEffCorrRatioGeneral==0 || graphEffCorrRatioEarly==0 )  {
    std::cerr << "At least one graph not found. " << std::endl;
    return;
  }


  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  double yMin = 0.950;
  double yMax = 1.020;

  graphEffCorrRatioGeneral->SetLineColor(kBlack);
  graphEffCorrRatioGeneral->SetMarkerColor(kBlack);

  graphEffCorrRatioEarly->SetLineColor(kViolet-2);
  graphEffCorrRatioEarly->SetMarkerColor(kViolet-2);

  graphEffCorrRatioGeneral->SetLineWidth(2);
  graphEffCorrRatioGeneral->SetMarkerSize(1.6);

  graphEffCorrRatioEarly->SetLineWidth(2);
  graphEffCorrRatioEarly->SetMarkerSize(1.6);

  graphEffCorrRatioEarly->GetYaxis()->SetTitle("Data / Sim. ratio");
  graphEffCorrRatioEarly->GetYaxis()->SetRangeUser(yMin, yMax);
  graphEffCorrRatioEarly->GetYaxis()->SetNdivisions(508);
  graphEffCorrRatioEarly->GetYaxis()->SetDecimals(true);
  graphEffCorrRatioEarly->GetXaxis()->SetRangeUser(3., 20.);
  graphEffCorrRatioEarly->GetXaxis()->SetTitle("p_{T}");

  graphEffCorrRatioEarly->Draw("AP");

  TLine line(graphEffCorrRatioGeneral->GetX()[0]-graphEffCorrRatioGeneral->GetErrorXlow(0), 1, graphEffCorrRatioGeneral->GetX()[graphEffCorrRatioGeneral->GetN()-1]+graphEffCorrRatioGeneral->GetErrorXhigh(graphEffCorrRatioGeneral->GetN()-1), 1);
  line.SetLineWidth(4);
  line.SetLineColor(65);
  line.DrawClone("SAME");

  graphEffCorrRatioGeneral->Draw("PSAME");

  TLegend *leg;
  leg = new TLegend(.20, .20, .65, .30 );
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->AddEntry(graphEffCorrRatioGeneral, "generalTracks", "LP");
  leg->AddEntry(graphEffCorrRatioEarly, "earlyGeneralTracks", "LP");
  leg->Draw("same");

  TPaveText *cmsprel = new TPaveText(0.46 ,.96,.94,.99,"NDC");
  cmsprel->SetTextSize(0.05);
  cmsprel->SetFillColor(0);
  cmsprel->SetFillStyle(0);
  cmsprel->SetLineStyle(2);
  cmsprel->SetLineColor(0);
  cmsprel->SetTextAlign(12);
  cmsprel->SetTextFont(42);
  cmsprel->AddText(cmsconditionslabel);
  cmsprel->Draw("same");

  TPaveText *cmsprel2 = new TPaveText(0.20 ,.31,.27,.35,"NDC");
  cmsprel2->SetTextSize(0.04);
  cmsprel2->SetFillColor(0);
  cmsprel2->SetFillStyle(0);
  cmsprel2->SetLineStyle(2);
  cmsprel2->SetLineColor(0);
  cmsprel2->SetTextAlign(12);
  cmsprel2->SetTextFont(42);
  cmsprel2->AddText(cmslabel);
  cmsprel2->Draw("same");


  c1->Print("eff_pt_comparison_ratio.png", ".png");
  c1->Print("eff_pt_comparison_ratio.pdf", ".pdf");

}
