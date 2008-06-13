#ifndef MuonAnalysis_TagAndProbe_eTriggerCandProducer_h
#define MuonAnalysis_TagAndProbe_eTriggerCandProducer_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
// forward declarations

class eTriggerCandProducer : public edm::EDProducer 
{
 public:
  explicit eTriggerCandProducer(const edm::ParameterSet&);
  ~eTriggerCandProducer();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual bool TriggerDecision(edm::Event&, reco::GsfElectron&);
  virtual bool MatchObjects( const reco::Candidate* , 
			     const reco::GsfElectron&);

  // ----------member data ---------------------------
      
  std::string _inputProducer;
  edm::InputTag triggerEventTag_;
  edm::InputTag hltL1Tag_;
  edm::InputTag hltTag_;
  double delRMatchingCut_;

};

#endif
