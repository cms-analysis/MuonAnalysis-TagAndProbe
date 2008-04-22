// -*- C++ -*-
//
// Package:    Histogrammer
// Class:      Histogrammer
// 
/**\class Histogrammer Histogrammer.cc MuonAnalysis/Histogrammer/src/Histogrammer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Adam Hunt"
//         Created:  Sun Apr 20 10:35:25 CDT 2008
// $Id: Histogrammer.cc,v 1.1 2008/04/22 15:50:24 ahunt Exp $
//
//

#include "MuonAnalysis/TagAndProbe/interface/Histogrammer.h"

Histogrammer::Histogrammer(const edm::ParameterSet& iConfig)
{
  using namespace std;

  fileNames       = iConfig.getUntrackedParameter< std::vector<std::string> >("inputFileNames");
  quantities      = iConfig.getUntrackedParameter< std::vector<std::string> >("quantities"); 
  conditions      = iConfig.getUntrackedParameter< std::vector<std::string> >("conditions"); 
  outputFileNames = iConfig.getUntrackedParameter< std::vector<std::string> >("outputFileNames");
  XBins           = iConfig.getUntrackedParameter< std::vector<unsigned int> >("XBins");
  XMin            = iConfig.getUntrackedParameter< std::vector<double> >("XMin");
  XMax            = iConfig.getUntrackedParameter< std::vector<double> >("XMax");
  logY            = iConfig.getUntrackedParameter< std::vector<unsigned int> >("logY");

  // Chain files together
  std::string tempString;
  fChain = new TChain("evttree");

  vector<string>::iterator it;
  for( it=fileNames.begin(); it < fileNames.end(); it++){
    tempString = *it;
    fChain->Add(tempString.c_str());
  }

  vectorSize = quantities.size();

  Histograms = new TH1F[vectorSize];

}

Histogrammer::~Histogrammer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// ------------ method called to for each event  ------------
void
Histogrammer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  
  if(outputFileNames.size() != vectorSize){
    std::cout << "outputFileNames is not the same size as quantities" << std::endl;    
  }else if(conditions.size() != vectorSize){
    std::cout << "conditions is not the same size as quantities" << std::endl;    
  }else if(XBins.size() != vectorSize){
    std::cout << "XBins is not the same size as quantities" << std::endl;    
  }else if(XMax.size() != vectorSize){
    std::cout << "XMax is not the same size as quantities" << std::endl;    
  }else if(XMin.size() != vectorSize){
    std::cout << "XMin is not the same size as quantities" << std::endl;    
  }else if(logY.size() != vectorSize){
    std::cout << "logY is not the same size as quantities" << std::endl;    
  }else{
    
    for(unsigned int i = 0; i < vectorSize; i++){
      CreateHistogram(Histograms[i], i);
      SaveHistogram(Histograms[i], i);
    }
  }
}

int Histogrammer::CreateHistogram(TH1F& Histo, int i){
  
  std::stringstream tempsstream;
  std::stringstream HistoName;

  // initialize
  HistoName.str(std::string());
  HistoName << "Histo" << i;

  tempsstream.str(std::string());    
  tempsstream << quantities[i] << " >>  " << HistoName.str();  

  if(fChain->FindBranch(quantities[i].substr(0,quantities[i].find("[")).c_str())){ 
    Histo =  TH1F(HistoName.str().c_str(), "", XBins[i], XMin[i], XMax[i]);
    fChain->Draw(tempsstream.str().c_str());
  }else{
    std::cout << "Branch does not exist: " << quantities[i] << std::endl;
  }

  return 0;
}

int Histogrammer::SaveHistogram(TH1F& Histo, int i){
  
  TCanvas* c1 = new TCanvas("c1","c1",700,500);
  c1->GetPad(0)->SetTicks(1,1);
  c1->SetLogy(logY[i]);
  
  Histo.Draw();
  
  c1->SaveAs(outputFileNames[i].c_str());
  
  delete c1;

  return 0;
}


void Histogrammer::SideBandSubtraction(const TH1F* Total, TH1F* &Result, Double_t Peak, Double_t SD){

  const Double_t BinWidth  = Total->GetXaxis()->GetBinWidth(1);
  const Int_t nbins = Total->GetNbinsX();
  const Double_t xmin = Total->GetXaxis()->GetXmin();

  const Int_t PeakBin = (Int_t)(Peak - xmin)/BinWidth + 1; // Peak
  const Double_t SDBin = SD/BinWidth; // Standard deviation
  const Double_t I = 3*SDBin; // Interval
  const Double_t D = 5*SDBin;  // Distance from peak

  const Double_t IntegralRight = Total->Integral(PeakBin + D, PeakBin + D + I);
  const Double_t IntegralLeft = Total->Integral(PeakBin - D - I, PeakBin - D);

  std::cout << "Peak: " << Peak << " " << PeakBin << std::endl;
  std::cout << "SD: " << SD << " " << SDBin << std::endl;
  std::cout << IntegralRight << std::endl;
  std::cout << IntegralLeft << std::endl;

  Double_t SubValue = 0.0;
  Double_t NewValue = 0.0;

  for(Int_t bin = 1; bin < (nbins + 1); bin++){
    SubValue = ((IntegralRight - IntegralLeft)/(2*D+I)*(bin - PeakBin - D - I/2.0) + IntegralRight)/I;
    NewValue = Total->GetBinContent(bin)-SubValue;
    std::cout << Total->GetBinContent(bin)  << " - " << SubValue <<  " = " << NewValue << std::endl;
    if(NewValue > 0){
      Result->SetBinContent(bin, NewValue);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
Histogrammer::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Histogrammer::endJob() {
}

//define this as a plug-in

