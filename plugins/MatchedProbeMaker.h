#ifndef MuonAnalysis_TagAndProbe_MatchedProbeMaker_H
#define MuonAnalysis_TagAndProbe_MatchedProbeMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// class decleration
//

class MatchedProbeMaker : public edm::EDProducer {
public:
  explicit MatchedProbeMaker(const edm::ParameterSet&);
  ~MatchedProbeMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag m_candidateSource;
  edm::InputTag m_referenceSource;
  edm::InputTag m_resMatchMapSource;

};
#endif
