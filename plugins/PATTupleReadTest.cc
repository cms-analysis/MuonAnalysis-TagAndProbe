// -*- C++ -*-
//
// Package:    PATTupleReadTest
// Class:      PATTupleReadTest
// 
/**\class PATTupleReadTest PATTupleReadTest.cc MuonAnalysis/TagAndProbe/src/PATTupleReadTest.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tommaso Boccali
// Modified for muons: Jonathan Hollar
//         Created:  Tue Nov 25 15:50:50 CET 2008
// $Id: PATTupleReadTest.cc,v 1.5 2010/06/02 11:37:28 jjhollar Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//
// class decleration
//


#include "DataFormats/MuonReco/interface/Muon.h"   
#include "DataFormats/MuonReco/interface/MuonFwd.h"   

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TFile.h> 
#include <TH1D.h> 
#include <TTree.h> 
#include <TRandom3.h>

class PATTupleReadTest : public edm::EDAnalyzer {
public:
  explicit PATTupleReadTest(const edm::ParameterSet&);
  ~PATTupleReadTest();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
};

using namespace edm; 
using namespace std; 
using namespace reco; 

//
// constructors and destructor
//
PATTupleReadTest::PATTupleReadTest(const edm::ParameterSet& iConfig)

{
}


PATTupleReadTest::~PATTupleReadTest()
{
}

// ------------ method called to for each event  ------------
void
PATTupleReadTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<pat::Muon> > patmuons;
  iEvent.getByLabel("selectedPatMuonsWithEff",patmuons);
  edm::View<pat::Muon>::const_iterator patmuon;
  for (edm::View<pat::Muon>::const_iterator patmuon = patmuons->begin(), end = patmuons->end(); patmuon != end; ++patmuon) 
    {
      std::string effname = "GlobalMuon_Data_CaloMuonProbe_JPsi";
      pat::Muon myMuon = *patmuon; // copy

      cout << "\tRead PAT efficiency " 
	   << (myMuon.efficiency("GlobalMuon_Data_CaloMuonProbe_JPsi")).value() << " + " 
	   << (myMuon.efficiency("GlobalMuon_Data_CaloMuonProbe_JPsi_UpperError")).value() << " - " 
	   << (myMuon.efficiency("GlobalMuon_Data_CaloMuonProbe_JPsi_LowerError")).value() 
	   << " (" <<  effname << ")" << endl;
    }
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATTupleReadTest::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATTupleReadTest::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATTupleReadTest);
