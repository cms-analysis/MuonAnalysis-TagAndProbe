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
// $Id: MuTestPerformanceFW_ES.cc,v 1.1 2009/09/30 09:47:26 jjhollar Exp $
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

#include "RecoMuon/Records/interface/MuonPerformanceRecord.h"

#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h" 
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
  std::string name1;
  std::string name2;
  std::string rootfilename;
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  double getEff(double, double, double, int, const MuonPerformance &);
  double getEffError(double, double, double, int, const MuonPerformance &); 
  bool passesPIDKilling(double, double, double, int, const MuonPerformance &);

  // ----------member data ---------------------------

  TFile *thefile;
  TH1D *hMCtruth;
  TH1D *hMCCorrected;
  TH1D *hReconstructed;
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
  name1 =  iConfig.getParameter<std::string>("AlgoName1");
  name2 =  iConfig.getParameter<std::string>("AlgoName2");
  rootfilename = iConfig.getUntrackedParameter<std::string>("outfilename","test.root");
}


MuTestPerformanceFW_ES::~MuTestPerformanceFW_ES()
{
}


//
// member functions
//
double
MuTestPerformanceFW_ES::getEff(double pt, double eta, double phi, int charge, const MuonPerformance & theperf)
{
  BinningPointByMap p; 
  p.insert(BinningVariables::MuonEta,eta); 
  p.insert(BinningVariables::MuonPt,pt); 
  p.insert(BinningVariables::MuonPhi,phi); 
  p.insert(BinningVariables::MuonCharge,charge);  

  double theeff = theperf.getResult(PerformanceResult::MUEFF,p);  

  /* Debugging printout. */ 
  std::cout <<" test eta=" << eta << ", pt=" << pt << ", phi = " << phi << ", " << "charge = " << charge << std::endl;   
  std::cout <<" Discriminant is "<<theperf.workingPoint().discriminantName()<<std::endl;  
  std::cout <<" mueff/muerr ="<<theperf.getResult(PerformanceResult::MUEFF,p)<<"/"<<theperf.getResult(PerformanceResult::MUERR,p)<<std::endl;  

 
  return theeff;
}

double 
MuTestPerformanceFW_ES::getEffError(double pt, double eta, double phi, int charge, const MuonPerformance & theperf) 
{ 
  BinningPointByMap p;  
  p.insert(BinningVariables::MuonEta,eta);  
  p.insert(BinningVariables::MuonPt,pt);  
  p.insert(BinningVariables::MuonPhi,phi);  
  p.insert(BinningVariables::MuonCharge,charge);   
 
  double theefferr = theperf.getResult(PerformanceResult::MUERR,p);  
  
  return theefferr; 
} 

bool
MuTestPerformanceFW_ES::passesPIDKilling(double pt, double eta, double phi, int charge, const MuonPerformance & theperf)
{
  bool passes = false;

  BinningPointByMap p;   
  p.insert(BinningVariables::MuonEta,eta);   
  p.insert(BinningVariables::MuonPt,pt);   
  p.insert(BinningVariables::MuonPhi,phi);   
  p.insert(BinningVariables::MuonCharge,charge);    
  
  double theeff = theperf.getResult(PerformanceResult::MUEFF,p);   

  double randthrow = rnd->Rndm(); 
 
  if(randthrow < theeff) 
    { 
      passes = true;
    }

  return passes;
}

// ------------ method called to for each event  ------------
void
MuTestPerformanceFW_ES::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* Retrieve the DB "MuonPerformance" record through the EventSetup. */
  edm::ESHandle<MuonPerformance> perfSA;
  iSetup.get<MuonPerformanceRecord>().get(name1,perfSA);
  const MuonPerformance & standaloneMuon = *(perfSA.product());

  edm::ESHandle<MuonPerformance> perfTRK; 
  iSetup.get<MuonPerformanceRecord>().get(name2,perfTRK); 
  const MuonPerformance & tracking = *(perfTRK.product()); 

  double eta, pt, phi;
  int pdgid = 0;
  int status = 0;
  int chg = 0;

  Handle<GenParticleCollection> genParticles;  
  iEvent.getByLabel( "genParticles", genParticles );  

  for ( size_t i = 0; i < genParticles->size(); ++ i )  
    { 
      const Candidate & p = (*genParticles)[ i ]; 

      status = p.status(); 
      if(p.status() == 3) 
	continue; 

      pt = p.pt();
      eta = p.eta(); 
      phi = p.phi();
      pdgid=p.pdgId(); 
      
      /* Look at all true muons. */
      if(pdgid == 13 || pdgid == -13)
	{
	  if(pdgid == 13) 
	    chg = 1;
	  if(pdgid == -13) 
	    chg = -1;

	  /* Now lookup the efficiency and error from the DB based on the muon pT, eta, phi, and charge */

	  getEff(pt, eta, phi, chg, standaloneMuon);
	  getEffError(pt, eta, phi, chg, standaloneMuon); 
          getEff(pt, eta, phi, chg, tracking); 
          getEffError(pt, eta, phi, chg, tracking);  


	  /* Fill a histogram of the pT spectrum for all true muons */
          hMCtruth->Fill(pt); 

	  /* 
	   * Now apply the actual efficiency correction. The example here uses "PID-Killing"; just do a muon-by-muon 
	   * accept-reject and keep all muons that pass. Alternativly "PID-weighting"; when we have tables for the T/P 
	   * efficiency from both data & MC, apply the ratio as a weight.
	   */
	  if(passesPIDKilling(pt, eta, phi, chg, standaloneMuon)) // Apply the SA-muon efficiency correction
	    {
	      if(passesPIDKilling(pt, eta, phi, chg, tracking)) // Apply the tracking efficiency correction
		{
		  hMCCorrected->Fill(pt);
		}
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
    }

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuTestPerformanceFW_ES::beginJob(const edm::EventSetup&)
{
  thefile = new TFile(rootfilename.c_str(),"recreate");  
  thefile->cd();  
  hMCtruth= new TH1D("hMCtruth","hMCtruth",100,0,100); 
  hMCCorrected= new TH1D("hMCCorrected","hMCCorrected",100,0,100);  
  hReconstructed= new TH1D("hReconstructed","hReconstructed",100,0,100);   
  rnd = new TRandom3(); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuTestPerformanceFW_ES::endJob() {
  hMCtruth->Write(); 
  hMCCorrected->Write(); 
  hReconstructed->Write();
  thefile->Close(); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuTestPerformanceFW_ES);
