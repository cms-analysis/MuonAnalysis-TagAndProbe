// -*- C++ -*-
//
// Package:    MuTestPAT
// Class:      MuTestPAT
// 
/**\class MuTestPAT MuTestPAT.cc MuonAnalysis/TagAndProbe/src/MuTestPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tommaso Boccali
// Modified for muons: Jonathan Hollar
//         Created:  Tue Nov 25 15:50:50 CET 2008
// $Id: MuTestPAT.cc,v 1.6 2010/06/02 11:37:28 jjhollar Exp $
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

#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h"
#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"  

#include "DataFormats/MuonReco/interface/Muon.h"   
#include "DataFormats/MuonReco/interface/MuonFwd.h"   

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TFile.h> 
#include <TH1D.h> 
#include <TTree.h> 
#include <TRandom3.h>

class MuTestPAT : public edm::EDProducer {
public:
  explicit MuTestPAT(const edm::ParameterSet&);
  ~MuTestPAT();
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
private:
  std::vector<std::string> algonames;
  std::string rootfilename;
  bool useAbsEtaVals;
  virtual void beginJob() ;
  //  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  MuonPerformanceReadback *effreader;
};

using namespace edm; 
using namespace std; 

using namespace reco; 

//
// constructors and destructor
//
MuTestPAT::MuTestPAT(const edm::ParameterSet& iConfig)

{
  algonames =  iConfig.getParameter< std::vector<std::string> >("AlgoNames");
  useAbsEtaVals = iConfig.getParameter< bool >("UseAbsEtaVals"); 

  produces<std::vector<pat::Muon> >();

}


MuTestPAT::~MuTestPAT()
{
}

// ------------ method called to for each event  ------------
void
MuTestPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::string algoname;
  double pateta, patpt, patphi;
  int patchg = 0;

  /* Fill PAT efficiencies from the DB */
  std::vector<pat::Muon> * patCorrMuons = new std::vector<pat::Muon>();

  edm::Handle<edm::View<pat::Muon> > patmuons;
  iEvent.getByLabel("selectedPatMuons",patmuons);
  edm::View<pat::Muon>::const_iterator patmuon;

  for (edm::View<pat::Muon>::const_iterator patmuon = patmuons->begin(), end = patmuons->end(); patmuon != end; ++patmuon) 
    {
      patpt = patmuon->pt();
      pateta = patmuon->eta();
      patphi = patmuon->phi();
      patchg = patmuon->charge();
      pat::Muon myMuon = *patmuon; // copy

      if(useAbsEtaVals == true)
	pateta = fabs(pateta);

      sort(algonames.begin(), algonames.end());

      for(unsigned int i = 0; i < algonames.size(); ++i)
	{
	  std::string effname = algonames[i]; 
	  std::string effnameuppererror = effname + "_UpperError";
	  std::string effnamelowererror = effname + "_LowerError";

	  const MuonPerformance &muonefficiency = effreader->getPerformanceRecord(effname, iSetup);
	  const MuonPerformance &muonefficiencyuppererror = effreader->getPerformanceRecord(effnameuppererror, iSetup); 
          const MuonPerformance &muonefficiencylowererror = effreader->getPerformanceRecord(effnamelowererror, iSetup); 

	  double myeff = effreader->getEff(patpt, pateta, patphi, patchg, muonefficiency);
	  double myefflowererr = effreader->getEff(patpt, pateta, patphi, patchg, muonefficiencylowererror); 
          double myeffuppererr = effreader->getEff(patpt, pateta, patphi, patchg, muonefficiencyuppererror);  
	  double myefferr = 0.0;

	  cout << "\tDB efficiency " 
	       << myeff << " + "
	       << myeffuppererr << " - " 
	       << myefflowererr << " (" 
	       << effname << ")" << endl;

	  myMuon.setEfficiency(effname, pat::LookupTableRecord(myeff, myefferr, 0)); 
	  myMuon.setEfficiency(effnamelowererror, pat::LookupTableRecord(myefflowererr, myefferr, 0)); 
          myMuon.setEfficiency(effnameuppererror, pat::LookupTableRecord(myeffuppererr, myefferr, 0)); 
	  
	  // Fill asymmetric errors
	  //      myMuon.setEfficiency(effname, pat::LookupTableRecord(myeff, myefferr, myefferr, myefferr, 0));
	}
      patCorrMuons->push_back(myMuon);
    }

  std::auto_ptr<std::vector<pat::Muon> > ptr(patCorrMuons);
  iEvent.put(ptr);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuTestPAT::beginJob()
{
  effreader = new MuonPerformanceReadback();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuTestPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuTestPAT);
