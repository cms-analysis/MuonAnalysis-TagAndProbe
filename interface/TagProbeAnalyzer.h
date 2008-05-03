#ifndef MuonAnalysis_TagAndProbe_TagProbeAnalyzer_h
#define MuonAnalysis_TagAndProbe_TagProbeAnalyzer_h
// -*- C++ -*-
//
// Package:     TagAndProbe
// Class  :     TagProbeAnalyzer
// 
/**\class TagAndProbe TagProbeAnalyzer.h MuonAnalysis/TagAndProbe/interface/TagProbeAnalyzer.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Nadia Adam
//         Created:  Fri Dec 28 13:22:23 CST 2007
// $Id: TagProbeAnalyzer.h,v 1.3 2008/04/14 19:16:08 neadam Exp $
//

// system include files
#include <memory>
#include <string>

// user include files
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/TagAndProbe/interface/FillEvtTree.h"

// forward declarations

class TagProbeAnalyzer : public edm::EDAnalyzer 
{
   public:
      explicit TagProbeAnalyzer(const edm::ParameterSet&);
      ~TagProbeAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // Functions to fill the varios event data sets
      void fillRunEventInfo();
      void fillTriggerInfo();
      void fillMCInfo();
      void fillTrackInfo();
      void fillTagProbeInfo();
      void fillTrueEffInfo();
      void fillVertexInfo();

      bool CandFromZ( const reco::CandidateRef& cand, 
		      edm::Handle<reco::CandMatchMap>& matchMap );

      int ProbePassProbeOverlap( const reco::CandidateRef& probe, 
				 edm::Handle<reco::CandidateCollection>& passprobes );

      bool MatchObjects( const reco::Candidate *hltObj, 
			 const reco::CandidateRef& tagObj );

      // ----------member data ---------------------------
      const edm::Event* m_event;
      const edm::EventSetup* m_setup;

      EvtTree     m_evtTree;
      FillEvtTree m_fillEvt;
  
      std::string rootFile_;

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

      // Vertexing
      edm::InputTag PvxtTag_;

      // Trigger parameters
      edm::InputTag triggerEventTag_;
      edm::InputTag hltL1Tag_;
      edm::InputTag hltTag_;
      double delRMatchingCut_;
      
      // MC parameter
      bool isMC_;

      math::XYZPoint* m_vtx;
};


#endif
