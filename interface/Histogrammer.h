// system include files
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>

class Histogrammer : public edm::EDAnalyzer{
 public:
  explicit Histogrammer(const edm::ParameterSet&);
  ~Histogrammer();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int CreateHistogram(TH1F& Histo, int i);
  int SaveHistogram(TH1F& Histo, std::string outFileName, Int_t LogY);

  void SideBandSubtraction(const TH1F& Total, TH1F& Result, Double_t Peak, Double_t SD);
  void SBSExample();


  std::vector<std::string> fileNames;
  std::vector<std::string> quantities;  
  std::vector<std::string> conditions;  
  std::vector<std::string> outputFileNames;

  std::vector<unsigned int> XBins;
  std::vector<double> XMax;
  std::vector<double> XMin;
  std::vector<unsigned int> logY;

  TChain* fChain;
  TH1F* Histograms;
  int* NumEvents;

  
  unsigned int vectorSize;
  
};
