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
// $Id$
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

      bool CandFromZ( const reco::CandidateRef& cand, 
		      edm::Handle<reco::CandMatchMap>& matchMap );

      int ProbePassProbeOverlap( const reco::CandidateRef& probe, 
				 edm::Handle<reco::CandidateCollection>& passprobes );

      bool MatchObjects( const reco::Candidate *hltObj, 
			 const reco::CandidateRef& tagObj );
      
      // ----------member data ---------------------------
      const edm::Event* m_event;
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
      double delRMatchingCut_;
      
      // MC parameter
      bool isMC_;

      // Private EDM Ntuple variables
      std::auto_ptr< int > run_;
      std::auto_ptr< int > event_;

      std::auto_ptr< int > nl1_;
      std::auto_ptr< int > nhlt_;

      std::auto_ptr< int > nrtp_;
      std::auto_ptr< std::vector<int> > tp_type_;
      std::auto_ptr< std::vector<int> > tp_true_;
      std::auto_ptr< std::vector<int> > tp_ppass_;
      std::auto_ptr< std::vector<int> > tp_l1_;  
      std::auto_ptr< std::vector<int> > tp_hlt_; 

      std::auto_ptr< std::vector<float> > tp_mass_;
      std::auto_ptr< std::vector<float> > tp_p_;  
      std::auto_ptr< std::vector<float> > tp_pt_; 
      std::auto_ptr< std::vector<float> > tp_px_; 
      std::auto_ptr< std::vector<float> > tp_py_; 
      std::auto_ptr< std::vector<float> > tp_pz_; 
      std::auto_ptr< std::vector<float> > tp_e_;  
      std::auto_ptr< std::vector<float> > tp_et_; 

      std::auto_ptr< std::vector<float> > tp_tag_p_;   
      std::auto_ptr< std::vector<float> > tp_tag_px_;  
      std::auto_ptr< std::vector<float> > tp_tag_py_;  
      std::auto_ptr< std::vector<float> > tp_tag_pz_;  
      std::auto_ptr< std::vector<float> > tp_tag_pt_;  
      std::auto_ptr< std::vector<float> > tp_tag_e_;   
      std::auto_ptr< std::vector<float> > tp_tag_et_;  
      std::auto_ptr< std::vector<float> > tp_tag_q_;   
      std::auto_ptr< std::vector<float> > tp_tag_eta_; 
      std::auto_ptr< std::vector<float> > tp_tag_phi_; 

      std::auto_ptr< std::vector<float> > tp_probe_p_;   
      std::auto_ptr< std::vector<float> > tp_probe_px_;  
      std::auto_ptr< std::vector<float> > tp_probe_py_;  
      std::auto_ptr< std::vector<float> > tp_probe_pz_;  
      std::auto_ptr< std::vector<float> > tp_probe_pt_;  
      std::auto_ptr< std::vector<float> > tp_probe_e_;   
      std::auto_ptr< std::vector<float> > tp_probe_et_;  
      std::auto_ptr< std::vector<float> > tp_probe_q_;   
      std::auto_ptr< std::vector<float> > tp_probe_eta_; 
      std::auto_ptr< std::vector<float> > tp_probe_phi_; 

      std::auto_ptr< int > ncnd_;
      std::auto_ptr< std::vector<int> > cnd_type_;   
      std::auto_ptr< std::vector<int> > cnd_tag_;    
      std::auto_ptr< std::vector<int> > cnd_aprobe_; 
      std::auto_ptr< std::vector<int> > cnd_pprobe_; 
      std::auto_ptr< std::vector<int> > cnd_moid_;   
      std::auto_ptr< std::vector<int> > cnd_gmid_;   

      std::auto_ptr< std::vector<float> > cnd_p_;   
      std::auto_ptr< std::vector<float> > cnd_pt_;  
      std::auto_ptr< std::vector<float> > cnd_px_;  
      std::auto_ptr< std::vector<float> > cnd_py_;  
      std::auto_ptr< std::vector<float> > cnd_pz_;  
      std::auto_ptr< std::vector<float> > cnd_e_;   
      std::auto_ptr< std::vector<float> > cnd_et_;  
      std::auto_ptr< std::vector<float> > cnd_q_;   
      std::auto_ptr< std::vector<float> > cnd_eta_; 
      std::auto_ptr< std::vector<float> > cnd_phi_; 

      std::auto_ptr< std::vector<float> > cnd_rp_;   
      std::auto_ptr< std::vector<float> > cnd_rpt_;  
      std::auto_ptr< std::vector<float> > cnd_rpx_;  
      std::auto_ptr< std::vector<float> > cnd_rpy_;  
      std::auto_ptr< std::vector<float> > cnd_rpz_;  
      std::auto_ptr< std::vector<float> > cnd_re_;   
      std::auto_ptr< std::vector<float> > cnd_ret_;  
      std::auto_ptr< std::vector<float> > cnd_rq_;   
      std::auto_ptr< std::vector<float> > cnd_reta_; 
      std::auto_ptr< std::vector<float> > cnd_rphi_; 
      
};

#endif
