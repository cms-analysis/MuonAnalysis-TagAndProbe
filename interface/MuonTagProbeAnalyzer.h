#ifndef MuonAnalysis_TagAndProbe_MuonTagProbeAnalyzer_h
#define MuonAnalysis_TagAndProbe_MuonTagProbeAnalyzer_h
// -*- C++ -*-
//
// Package:     TagAndProbe
// Class  :     MuonTagProbeAnalyzer
// 
/**\class TagAndProbe MuonTagProbeAnalyzer.h MuonAnalysis/TagAndProbe/interface/MuonTagProbeAnalyzer.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Nadia Adam
//         Created:  Fri Dec 28 13:22:23 CST 2007
// $Id: MuonTagProbeAnalyzer.h,v 1.1 2008/04/02 19:59:02 neadam Exp $
//

// system include files
#include <memory>
#include <string>

// user include files
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/TagAndProbe/interface/FillEvtTree.h"

// forward declarations

class MuonTagProbeAnalyzer : public edm::EDAnalyzer 
{
   public:
      explicit MuonTagProbeAnalyzer(const edm::ParameterSet&);
      ~MuonTagProbeAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // Functions to fill the varios event data sets
      void fillRunEventInfo();
      void fillTriggerInfo();
      void fillMCInfo();
      void fillTrackInfo();
      void fillMuonInfo();
      void fillTagProbeInfo();
      void fillMuonEffInfo();
      void fillZInfo();
      void fillVertexInfo();

      bool MatchObjects( reco::CandidateBaseRef hltObj, const reco::MuonRef& tagObj );

      // ----------member data ---------------------------
      const edm::Event* m_event;
      const edm::EventSetup* m_setup;

      EvtTree     m_evtTree;
      FillEvtTree m_fillEvt;
  
      std::string rootFile_;

      std::vector<int> mcParticles_;
      std::vector<int> mcParents_;
  
      edm::InputTag tracksTag_;
      edm::InputTag satracksTag_;
      edm::InputTag ctracksTag_;

      // Muon collection tags
      edm::InputTag muonsTag_;
      std::vector<edm::InputTag> tagMuonTags_;
      std::vector<edm::InputTag> allProbeMuonTags_;
      std::vector<edm::InputTag> passProbeMuonTags_;

      // Tag probe map tags
      std::vector<edm::InputTag> tagProbeMapTags_;

      // Candidate collection tags
      edm::InputTag genParticlesTag_;
      std::vector<edm::InputTag> tagCandTags_;
      std::vector<edm::InputTag> allProbeCandTags_;
      std::vector<edm::InputTag> passProbeCandTags_;

      // Z tags
      std::vector<edm::InputTag> zCandTags_;

      // Truth matching
      std::vector<edm::InputTag> tagTruthMatchMapTags_;
      std::vector<edm::InputTag> allProbeTruthMatchMapTags_;
      std::vector<edm::InputTag> passProbeTruthMatchMapTags_;

      // Vertexing
      edm::InputTag PvxtTag_;

      // Trigger parameters
      edm::InputTag hltL1Tag_;
      edm::InputTag hltTag_;
      double delRMatchingCut_;
      
      // MC parameter
      bool isMC_;


      math::XYZPoint* m_vtx;
};


#endif
