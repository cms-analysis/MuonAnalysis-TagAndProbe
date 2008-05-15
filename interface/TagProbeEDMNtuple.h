#ifndef MuonAnalysis_TagAndProbe_TagProbeEDMNtuple_h
#define MuonAnalysis_TagAndProbe_TagProbeEDMNtuple_h
// -*- C++ -*-
//
// Package:     TagAndProbe
// Class  :     TagProbeEDMNtuples
// 
/**\class TagProbeEDMNtuples TagProbeEDMNtuple.h MuonAnalysis/TagAndProbe/interface/TagProbeEDMNtuple.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Mon May  5 09:05:35 CDT 2008
// $Id: TagProbeEDMNtuple.h,v 1.2 2008/05/11 11:56:29 neadam Exp $
//

// system include files
#include <memory>
#include <string>

// user include files
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class TagProbeEDMNtuple : public edm::EDProducer 
{
   public:
      typedef edm::Ref< reco::CandidateView > CandViewRef;

      explicit TagProbeEDMNtuple(const edm::ParameterSet&);
      ~TagProbeEDMNtuple();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // Functions to fill the varios event data sets
      void fillRunEventInfo();
      void fillTriggerInfo();
      void fillMCInfo();
      void fillTrackInfo();
      void fillTagProbeInfo();
      void fillTrueEffInfo();

      bool CandFromZ( reco::GenParticleRef mcRef );

      int ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, 
				 edm::Handle<reco::CandidateView>& passprobes );

      bool MatchObjects( const reco::Candidate *hltObj, 
			 const reco::CandidateBaseRef& tagObj,
			 bool exact = true );
      
      // ----------member data ---------------------------
      edm::Event* m_event;
      const edm::EventSetup* m_setup;

      // Type of Cands (used for matching and PDGId)
      std::string candType_;

      // PDG id of Cands
      int candPDGId_;

      // MC particles to store
      std::vector<int> mcParticles_;
      std::vector<int> mcParents_;
  
      // Track Collection Tags
      std::vector<edm::InputTag> trackTags_;

      // Tag probe map tags
      std::vector<edm::InputTag> tagProbeMapTags_;

      // Candidate collection tags
      edm::InputTag genParticlesTag_;
      std::vector<edm::InputTag> tagCandTags_;
      std::vector<edm::InputTag> allProbeCandTags_;
      std::vector<edm::InputTag> passProbeCandTags_;

      // Truth matching
      std::vector<edm::InputTag> tagTruthMatchMapTags_;
      std::vector<edm::InputTag> allProbeTruthMatchMapTags_;
      std::vector<edm::InputTag> passProbeTruthMatchMapTags_;

      // Trigger parameters
      edm::InputTag triggerEventTag_;
      edm::InputTag hltL1Tag_;
      edm::InputTag hltTag_;

      // Matching parameters
      double delRMatchingCut_;
      double delPtRelMatchingCut_;
      
      // MC parameter
      bool isMC_;      
};

#endif
