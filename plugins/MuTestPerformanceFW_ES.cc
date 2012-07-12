// -*- C++ -*-
//
// Package:    MuTestPerformanceFW_ES
// Class:      MuTestPerformanceFW_ES
// 
/**\class MuTestPerformanceFW_ES MuTestPerformanceFW_ES.cc MuonAnalysis/TagAndProbe/src/MuTestPerformanceFW_ES.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tommaso Boccali
// Modified for muons: Jonathan Hollar
//         Created:  Tue Nov 25 15:50:50 CET 2008
// $Id: MuTestPerformanceFW_ES.cc,v 1.8 2011/01/30 10:54:00 gpetrucc Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

#include <TFile.h> 
#include <TH1D.h> 
#include <TTree.h> 
#include <TRandom3.h>

class MuTestPerformanceFW_ES : public edm::EDAnalyzer {
public:
  explicit MuTestPerformanceFW_ES(const edm::ParameterSet&);
  ~MuTestPerformanceFW_ES();
  
  
private:
  std::vector<std::string> algonames;
  std::string rootfilename;
  bool useAbsEtaVals;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  MuonPerformanceReadback *effreader;

  TFile *thefile;
  TH1D *hMCtruth;
  TH1D *hMCCorrected;
  TH1D *hReconstructed;
  TH1D *hMCtruthEta; 
  TH1D *hMCCorrectedEta; 
  TH1D *hReconstructedEta; 
  TRandom3 *rnd;
};

using namespace edm; 
using namespace std; 
using namespace HepMC; 
using namespace reco; 

//
// constructors and destructor
//
MuTestPerformanceFW_ES::MuTestPerformanceFW_ES(const edm::ParameterSet& iConfig)

{
  algonames =  iConfig.getParameter< std::vector<std::string> >("AlgoNames");
  rootfilename = iConfig.getUntrackedParameter<std::string>("outfilename","test.root");
  useAbsEtaVals = iConfig.getParameter< bool >("UseAbsEtaVals");  
}


MuTestPerformanceFW_ES::~MuTestPerformanceFW_ES()
{
}

// ------------ method called to for each event  ------------
void
MuTestPerformanceFW_ES::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* Retrieve the DB "MuonPerformance" record through the EventSetup. */

  double eta, pt, phi;
  int pdgid = 0;
  int chg = 0;
  std::string algoname;

  Handle<GenParticleCollection> genParticles;  
  iEvent.getByLabel( "genParticles", genParticles );  

  for ( size_t i = 0; i < genParticles->size(); ++ i )  
    { 
      const Candidate & p = (*genParticles)[ i ]; 

      if(p.status() == 3) 
	continue; 

      pt = p.pt();
      eta = p.eta(); 
      phi = p.phi();
      pdgid=p.pdgId(); 

      if(useAbsEtaVals == true)
	eta = fabs(eta);
      
      /* Look at all true muons. */
      if(pdgid == 13 || pdgid == -13)
	{
	  if(pdgid == 13) 
	    chg = 1;
	  if(pdgid == -13) 
	    chg = -1;

	  bool passescorrection = true;

	  /* Now lookup the efficiency and error from the DB based on the muon pT, eta, phi, and charge */
	  for(unsigned int i = 0; i < algonames.size(); ++i) 
	    {
	      std::string effname = algonames[i];  
	      std::string effnameuppererror = effname + "_UpperError"; 
	      std::string effnamelowererror = effname + "_LowerError"; 
 
	      const MuonPerformance &muonefficiency = effreader->getPerformanceRecord(effname, iSetup); 
	      //const MuonPerformance &muonefficiencyuppererror = effreader->getPerformanceRecord(effnameuppererror, iSetup);  
	      //const MuonPerformance &muonefficiencylowererror = effreader->getPerformanceRecord(effnamelowererror, iSetup);  

	      //double myeff = effreader->getEff(pt, eta, phi, chg, muonefficiency); 
	      //double myefflowererr = effreader->getEff(pt, eta, phi, chg, muonefficiencylowererror);  
	      //double myeffuppererr = effreader->getEff(pt, eta, phi, chg, muonefficiencyuppererror);   
 
	      /* Fill a histogram of the pT spectrum for all true muons */
	      hMCtruth->Fill(pt); 
	      hMCtruthEta->Fill(eta);  

	      /* 
	       * Now apply the actual efficiency correction. The example here uses "PID-Killing"; just do a muon-by-muon 
	       * accept-reject and keep all muons that pass. This is just an example. Alternativly "PID-weighting"; when we have 
	       * tables for the T/P efficiency from both data & MC, apply the ratio as a weight, etc, etc.
	       */
	      if(effreader->passesPIDKilling(pt, eta, phi, chg, muonefficiency) == false) // Apply the muon efficiency correction
		{
		  passescorrection = false;
		}
	    }
	  if(passescorrection == true)
	    {
	      hMCCorrected->Fill(pt); 
	      hMCCorrectedEta->Fill(eta);  
	    }
	}
    }

  /* For comparison, also look at the FullSIM/RECO muons */
  Handle<reco::MuonCollection> muons; 
  iEvent.getByLabel("muons", muons); 
  reco::MuonCollection::const_iterator muon1; 
  for (muon1 = muons->begin(); muon1 != muons->end(); ++muon1) 
    {
      hReconstructed->Fill(muon1->pt());
      hReconstructedEta->Fill(muon1->eta()); 
    }
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuTestPerformanceFW_ES::beginJob()
{
  effreader = new MuonPerformanceReadback();

  thefile = new TFile(rootfilename.c_str(),"recreate");  
  thefile->cd();  
  hMCtruth= new TH1D("hMCtruth","hMCtruth",10,0,20); 
  hMCCorrected= new TH1D("hMCCorrected","hMCCorrected",10,0,20);  
  hReconstructed= new TH1D("hReconstructed","hReconstructed",10,0,20);   
  hMCtruthEta= new TH1D("hMCtruthEta","hMCtruthEta",10,-3,3);  
  hMCCorrectedEta= new TH1D("hMCCorrectedEta","hMCCorrectedEta",10,-3,3);   
  hReconstructedEta= new TH1D("hReconstructedEta","hReconstructedEta",10,-3,3);    

  rnd = new TRandom3(); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuTestPerformanceFW_ES::endJob() {
  hMCtruth->Write(); 
  hMCCorrected->Write(); 
  hReconstructed->Write();
  hMCtruthEta->Write();  
  hMCCorrectedEta->Write();  
  hReconstructedEta->Write(); 
  thefile->Close(); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuTestPerformanceFW_ES);
