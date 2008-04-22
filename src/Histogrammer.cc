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
  
  unsigned int vectorSize = quantities.size();

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
    std::stringstream tempsstream;
    TCanvas* c1 = new TCanvas("c1","c1",700,500);
    c1->GetPad(0)->SetTicks(1,1);
    TH1F* Histo1;
    

    for(unsigned int i = 0; i < vectorSize; i++){
      if(fChain->FindBranch(quantities[i].substr(0,quantities[i].find("[")).c_str())){ 
	c1->SetLogy(logY[i]);
	tempsstream.str(std::string());  
	tempsstream << quantities[i] << " >> Histo1";  
	Histo1 = new TH1F("Histo1","",XBins[i], XMin[i], XMax[i]);
	fChain->Draw(tempsstream.str().c_str());
	c1->SaveAs(outputFileNames[i].c_str());
	delete Histo1;
      }else{
	std::cout << "Branch does not exist: " << quantities[i] << std::endl;
      }
    }
    delete c1;
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

