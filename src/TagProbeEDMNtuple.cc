// -*- C++ -*-
//
// Package:    TagProbeEDMNtuple
// Class:      TagProbeEDMNtuple
// 
/**\class TagProbeEDMNtuple TagProbeEDMNtuple.cc MuonAnalysis/TagProbeEDMNtuple/src/TagProbeEDMNtuple.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nadia Adam
//         Created:  Mon May  5 08:47:29 CDT 2008
// $Id$
//
//


// system include files
#include <cmath>

// user include files
#include "MuonAnalysis/TagAndProbe/interface/TagProbeEDMNtuple.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "AnalysisDataFormats/TagAndProbe/interface/CandidateAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/FitResult.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

//
// constants, enums and typedefs
//
using namespace std; 
using namespace reco; 
using namespace edm;
using namespace l1extra;
using namespace trigger;

//
// static data member definitions
//

//
// constructors and destructor
//
TagProbeEDMNtuple::TagProbeEDMNtuple(const edm::ParameterSet& iConfig)
{
   candType_ = iConfig.getUntrackedParameter<string>("tagProbeType","Muon");
   candPDGId_ = 13;
   if( candType_ != "Muon" ) candPDGId_ = 11;

   // Get the id's of any MC particles to store.
   vector<int> defaultPIDs;
   mcParticles_ = iConfig.getUntrackedParameter< vector<int> >(
      "mcParticles",defaultPIDs);

   vector<int> defaultPPIDs;
   mcParents_ = iConfig.getUntrackedParameter< vector<int> >(
      "mcParents",defaultPPIDs);

   if( mcParents_.size() != mcParticles_.size() )
   {
      mcParents_.clear();
      for(unsigned int i=0; i<mcParticles_.size(); ++i )
      {
	 mcParents_.push_back(0);
      }
   }

   // ********** Reco Tracks ********** //
   vector< edm::InputTag > dTrackTags;
   dTrackTags.push_back( edm::InputTag("ctfWithMaterialTracks"));
   trackTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "trackTags",dTrackTags);
   // ********************************* //

   // ********** Tag-Probes ********** //
   vector< edm::InputTag > defaultTagProbeMapTags;
   tagProbeMapTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "tagProbeMapTags",defaultTagProbeMapTags);
   for( int i=0; i<(int)tagProbeMapTags_.size(); ++i )
      cout << tagProbeMapTags_[i] << endl;
   // ******************************** //

   // Make sure vector sizes are correct!
   int map_size = (int)tagProbeMapTags_.size();
   const edm::InputTag dBlankTag("Blank");

   // ********** Candidate collections ********** //
   const edm::InputTag dGenParticlesTag("genParticleCandidates");
   genParticlesTag_ = iConfig.getUntrackedParameter<edm::InputTag>(
      "genParticlesTag",dGenParticlesTag);

   vector< edm::InputTag > defaultTagCandTags;
   tagCandTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "tagCandTags",defaultTagCandTags);
   for( int i=0; i<(int)tagCandTags_.size(); ++i ) 
      cout << tagCandTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)tagCandTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)tagCandTags_.size()); ++i ) 
	 tagCandTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultAllProbeCandTags;
   allProbeCandTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "allProbeCandTags",defaultAllProbeCandTags);
   for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
      cout << allProbeCandTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)allProbeCandTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)allProbeCandTags_.size()); ++i ) 
	 allProbeCandTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultPassProbeCandTags;
   passProbeCandTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "passProbeCandTags",defaultPassProbeCandTags);
   for( int i=0; i<(int)passProbeCandTags_.size(); ++i ) 
      cout << passProbeCandTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)passProbeCandTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)passProbeCandTags_.size()); ++i ) 
	 passProbeCandTags_.push_back(dBlankTag);
   } 
   // ******************************************* //

   // ********** Truth matching ********** //
   vector< edm::InputTag > defaultTagTruthMatchMapTags;
   tagTruthMatchMapTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "tagTruthMatchMapTags",defaultTagTruthMatchMapTags);
   for( int i=0; i<(int)tagTruthMatchMapTags_.size(); ++i ) 
      cout << tagTruthMatchMapTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)tagTruthMatchMapTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)tagTruthMatchMapTags_.size()); ++i ) 
	 tagTruthMatchMapTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultAllProbeTruthMatchMapTags;
   allProbeTruthMatchMapTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "allProbeTruthMatchMapTags",defaultAllProbeTruthMatchMapTags);
   for( int i=0; i<(int)allProbeTruthMatchMapTags_.size(); ++i ) 
      cout << allProbeTruthMatchMapTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)allProbeTruthMatchMapTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)allProbeTruthMatchMapTags_.size()); ++i ) 
	 allProbeTruthMatchMapTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultPassProbeTruthMatchMapTags;
   passProbeTruthMatchMapTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "passProbeTruthMatchMapTags",defaultPassProbeTruthMatchMapTags);
   for( int i=0; i<(int)passProbeTruthMatchMapTags_.size(); ++i )
      cout << passProbeTruthMatchMapTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)passProbeTruthMatchMapTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)passProbeTruthMatchMapTags_.size()); ++i ) 
	 passProbeTruthMatchMapTags_.push_back(dBlankTag);
   } 
   // ************************************ //

   // Trigger
   const edm::InputTag dTriggerEventTag("triggerSummaryRAW");
   triggerEventTag_ = 
      iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);

   const edm::InputTag dHLTL1Tag("SingleMuIsoLevel1Seed");
   hltL1Tag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltL1Tag",dHLTL1Tag);

   const edm::InputTag dHLTTag("SingleMuIsoL3IsoFiltered");
   hltTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltTag",dHLTTag);

   delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch",0.15);

   //register your products
   produces< int >( "run" ).setBranchAlias( "Run" );
   produces< int >( "event" ).setBranchAlias( "Event" );

   produces< int >( "nl1" ).setBranchAlias( "nL1" );
   produces< int >( "nhlt" ).setBranchAlias( "nHLT" );

   produces< int >( "nrTP" ).setBranchAlias( "nrTP" );
   produces< std::vector<int> >( "TPtype" ).setBranchAlias( "TPtype" );
   produces< std::vector<int> >( "TPtrue" ).setBranchAlias( "TPtrue" );
   produces< std::vector<int> >( "TPppass" ).setBranchAlias( "TPppass" );
   produces< std::vector<int> >( "TPl1" ).setBranchAlias( "TPl1" );
   produces< std::vector<int> >( "TPhlt" ).setBranchAlias( "TPhlt" );

   produces<std::vector<float> >( "TPmass" ).setBranchAlias( "TPmass" );
   produces<std::vector<float> >( "TPp"    ).setBranchAlias( "TPp"    );
   produces<std::vector<float> >( "TPpt"   ).setBranchAlias( "TPpt"   );
   produces<std::vector<float> >( "TPpx"   ).setBranchAlias( "TPpx"   );
   produces<std::vector<float> >( "TPpy"   ).setBranchAlias( "TPpy"   );
   produces<std::vector<float> >( "TPpz"   ).setBranchAlias( "TPpz"   );
   produces<std::vector<float> >( "TPe"    ).setBranchAlias( "TPe"    );
   produces<std::vector<float> >( "TPet"   ).setBranchAlias( "TPet"   );

   produces<std::vector<float> >( "TPTagp"   ).setBranchAlias( "TPTagp"   );
   produces<std::vector<float> >( "TPTagpx"  ).setBranchAlias( "TPTagpx"  );
   produces<std::vector<float> >( "TPTagpy"  ).setBranchAlias( "TPTagpy"  );
   produces<std::vector<float> >( "TPTagpz"  ).setBranchAlias( "TPTagpz"  );
   produces<std::vector<float> >( "TPTagpt"  ).setBranchAlias( "TPTagpt"  );
   produces<std::vector<float> >( "TPTage"   ).setBranchAlias( "TPTage"   );
   produces<std::vector<float> >( "TPTaget"  ).setBranchAlias( "TPTaget"  );
   produces<std::vector<float> >( "TPTagq"   ).setBranchAlias( "TPTagq"   );
   produces<std::vector<float> >( "TPTageta" ).setBranchAlias( "TPTageta" );
   produces<std::vector<float> >( "TPTagphi" ).setBranchAlias( "TPTagphi" );

   produces<std::vector<float> >( "TPProbep"   ).setBranchAlias( "TPProbep"   );
   produces<std::vector<float> >( "TPProbepx"  ).setBranchAlias( "TPProbepx"  );
   produces<std::vector<float> >( "TPProbepy"  ).setBranchAlias( "TPProbepy"  );
   produces<std::vector<float> >( "TPProbepz"  ).setBranchAlias( "TPProbepz"  );
   produces<std::vector<float> >( "TPProbept"  ).setBranchAlias( "TPProbept"  );
   produces<std::vector<float> >( "TPProbee"   ).setBranchAlias( "TPProbee"   );
   produces<std::vector<float> >( "TPProbeet"  ).setBranchAlias( "TPProbeet"  );
   produces<std::vector<float> >( "TPProbeq"   ).setBranchAlias( "TPProbeq"   );
   produces<std::vector<float> >( "TPProbeeta" ).setBranchAlias( "TPProbeeta" );
   produces<std::vector<float> >( "TPProbephi" ).setBranchAlias( "TPProbephi" );

   produces< int >( "nCnd" ).setBranchAlias( "nCnd" );
   produces< std::vector<int> >( "Cndtype" ).setBranchAlias( "Cndtype" );
   produces< std::vector<int> >( "Cndtag" ).setBranchAlias( "Cndtag" );
   produces< std::vector<int> >( "Cndaprobe" ).setBranchAlias( "Cndaprobe" );
   produces< std::vector<int> >( "Cndpprobe" ).setBranchAlias( "Cndpprobe" );
   produces< std::vector<int> >( "Cndmoid" ).setBranchAlias( "Cndmoid" );
   produces< std::vector<int> >( "Cndgmid" ).setBranchAlias( "Cndgmid" );

   produces<std::vector<float> >( "Cndp" ).setBranchAlias( "Cndp" );
   produces<std::vector<float> >( "Cndpt" ).setBranchAlias( "Cndpt" );
   produces<std::vector<float> >( "Cndpx" ).setBranchAlias( "Cndpx" );
   produces<std::vector<float> >( "Cndpy" ).setBranchAlias( "Cndpy" );
   produces<std::vector<float> >( "Cndpz" ).setBranchAlias( "Cndpz" );
   produces<std::vector<float> >( "Cnde" ).setBranchAlias( "Cnde" );
   produces<std::vector<float> >( "Cndet" ).setBranchAlias( "Cndet" );
   produces<std::vector<float> >( "Cndq" ).setBranchAlias( "Cndq" );
   produces<std::vector<float> >( "Cndeta" ).setBranchAlias( "Cndeta" );
   produces<std::vector<float> >( "Cndphi" ).setBranchAlias( "Cndphi" );

   produces<std::vector<float> >( "Cndrp" ).setBranchAlias( "Cndrp" );
   produces<std::vector<float> >( "Cndrpt" ).setBranchAlias( "Cndrpt" );
   produces<std::vector<float> >( "Cndrpx" ).setBranchAlias( "Cndrpx" );
   produces<std::vector<float> >( "Cndrpy" ).setBranchAlias( "Cndrpy" );
   produces<std::vector<float> >( "Cndrpz" ).setBranchAlias( "Cndrpz" );
   produces<std::vector<float> >( "Cndre" ).setBranchAlias( "Cndre" );
   produces<std::vector<float> >( "Cndret" ).setBranchAlias( "Cndret" );
   produces<std::vector<float> >( "Cndrq" ).setBranchAlias( "Cndrq" );
   produces<std::vector<float> >( "Cndreta" ).setBranchAlias( "Cndreta" );
   produces<std::vector<float> >( "Cndrphi" ).setBranchAlias( "Cndrphi" );


}


TagProbeEDMNtuple::~TagProbeEDMNtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TagProbeEDMNtuple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   cout << "Here in muon efficiency producer ..." << endl;

   // Set the private event & setup pointers
   m_event = &iEvent;
   m_setup = &iSetup;

   run_.reset( new int );
   event_.reset( new int );

   nl1_.reset( new int );
   nhlt_.reset( new int );

   nrtp_.reset( new int );
   tp_type_.reset( new vector<int> );
   tp_true_.reset( new vector<int> );
   tp_ppass_.reset( new vector<int> );
   tp_l1_.reset( new vector<int> );  
   tp_hlt_.reset( new vector<int> ); 

   tp_mass_.reset( new vector<float> );
   tp_p_.reset( new vector<float> );  
   tp_pt_.reset( new vector<float> ); 
   tp_px_.reset( new vector<float> ); 
   tp_py_.reset( new vector<float> ); 
   tp_pz_.reset( new vector<float> ); 
   tp_e_.reset( new vector<float> );  
   tp_et_.reset( new vector<float> ); 

   tp_tag_p_.reset( new vector<float> );   
   tp_tag_px_.reset( new vector<float> );  
   tp_tag_py_.reset( new vector<float> );  
   tp_tag_pz_.reset( new vector<float> );  
   tp_tag_pt_.reset( new vector<float> );  
   tp_tag_e_.reset( new vector<float> );   
   tp_tag_et_.reset( new vector<float> );  
   tp_tag_q_.reset( new vector<float> );   
   tp_tag_eta_.reset( new vector<float> ); 
   tp_tag_phi_.reset( new vector<float> ); 

   tp_probe_p_.reset( new vector<float> );   
   tp_probe_px_.reset( new vector<float> );  
   tp_probe_py_.reset( new vector<float> );  
   tp_probe_pz_.reset( new vector<float> );  
   tp_probe_pt_.reset( new vector<float> );  
   tp_probe_e_.reset( new vector<float> );   
   tp_probe_et_.reset( new vector<float> );  
   tp_probe_q_.reset( new vector<float> );   
   tp_probe_eta_.reset( new vector<float> ); 
   tp_probe_phi_.reset( new vector<float> ); 

   ncnd_.reset( new int );
   cnd_type_.reset( new vector<int> );   
   cnd_tag_.reset( new vector<int> );    
   cnd_aprobe_.reset( new vector<int> ); 
   cnd_pprobe_.reset( new vector<int> ); 
   cnd_moid_.reset( new vector<int> );   
   cnd_gmid_.reset( new vector<int> );   

   cnd_p_.reset( new vector<float> );   
   cnd_pt_.reset( new vector<float> );  
   cnd_px_.reset( new vector<float> );  
   cnd_py_.reset( new vector<float> );  
   cnd_pz_.reset( new vector<float> );  
   cnd_e_.reset( new vector<float> );   
   cnd_et_.reset( new vector<float> );  
   cnd_q_.reset( new vector<float> );   
   cnd_eta_.reset( new vector<float> ); 
   cnd_phi_.reset( new vector<float> ); 

   cnd_rp_.reset( new vector<float> );   
   cnd_rpt_.reset( new vector<float> );  
   cnd_rpx_.reset( new vector<float> );  
   cnd_rpy_.reset( new vector<float> );  
   cnd_rpz_.reset( new vector<float> );  
   cnd_re_.reset( new vector<float> );   
   cnd_ret_.reset( new vector<float> );  
   cnd_rq_.reset( new vector<float> );   
   cnd_reta_.reset( new vector<float> ); 
   cnd_rphi_.reset( new vector<float> ); 

   // Fill the run and event number information
   fillRunEventInfo();

   // Fill event level trigger info
   fillTriggerInfo();

   // Fill Tag-Probe Info
   fillTagProbeInfo();

   // Fill Efficiency Info for true muons
   fillTrueEffInfo();

   // Put the objects into the event
   iEvent.put( run_, "run" );
   iEvent.put( event_, "event" );

   iEvent.put( nl1_, "nl1" );
   iEvent.put( nhlt_, "nhlt" );

   iEvent.put( nrtp_, "nrTP" );

   iEvent.put( tp_type_, "TPtype" );
   iEvent.put( tp_true_, "TPtrue" );
   iEvent.put( tp_ppass_, "TPppass" );
   iEvent.put( tp_l1_, "TPl1" );  
   iEvent.put( tp_hlt_, "TPhlt" ); 

   iEvent.put( tp_mass_, "TPmass" );
   iEvent.put( tp_p_, "TPp" );  
   iEvent.put( tp_pt_, "TPpt" ); 
   iEvent.put( tp_px_, "TPpx" ); 
   iEvent.put( tp_py_, "TPpy" ); 
   iEvent.put( tp_pz_, "TPpz" ); 
   iEvent.put( tp_e_, "TPe" );  
   iEvent.put( tp_et_, "TPet" ); 

   iEvent.put( tp_tag_p_, "TPTagp" );   
   iEvent.put( tp_tag_px_, "TPTagpx" );  
   iEvent.put( tp_tag_py_, "TPTagpy" );  
   iEvent.put( tp_tag_pz_, "TPTagpz" );  
   iEvent.put( tp_tag_pt_, "TPTagpt" );  
   iEvent.put( tp_tag_e_, "TPTage" );   
   iEvent.put( tp_tag_et_, "TPTaget" );  
   iEvent.put( tp_tag_q_, "TPTagq" );   
   iEvent.put( tp_tag_eta_, "TPTageta" ); 
   iEvent.put( tp_tag_phi_, "TPTagphi" ); 

   iEvent.put( tp_probe_p_, "TPProbep" );   
   iEvent.put( tp_probe_px_, "TPProbepx" );  
   iEvent.put( tp_probe_py_, "TPProbepy" );  
   iEvent.put( tp_probe_pz_, "TPProbepz" );  
   iEvent.put( tp_probe_pt_, "TPProbept" );  
   iEvent.put( tp_probe_e_, "TPProbee" );   
   iEvent.put( tp_probe_et_, "TPProbeet" );  
   iEvent.put( tp_probe_q_, "TPProbeq" );   
   iEvent.put( tp_probe_eta_, "TPProbeeta" ); 
   iEvent.put( tp_probe_phi_, "TPProbephi" ); 

   iEvent.put( ncnd_, "nCnd" );

   iEvent.put( cnd_type_, "Cndtype" );
   iEvent.put( cnd_tag_, "Cndtag" );
   iEvent.put( cnd_aprobe_, "Cndaprobe" );
   iEvent.put( cnd_pprobe_, "Cndpprobe" );  
   iEvent.put( cnd_moid_, "Cndmoid" ); 
   iEvent.put( cnd_gmid_, "Cndgmid" ); 

   iEvent.put( cnd_p_, "Cndp" );  
   iEvent.put( cnd_pt_, "Cndpt" ); 
   iEvent.put( cnd_px_, "Cndpx" ); 
   iEvent.put( cnd_py_, "Cndpy" ); 
   iEvent.put( cnd_pz_, "Cndpz" ); 
   iEvent.put( cnd_e_, "Cnde" );  
   iEvent.put( cnd_et_, "Cndet" ); 
   iEvent.put( cnd_q_, "Cndq" );  
   iEvent.put( cnd_eta_, "Cndeta" ); 
   iEvent.put( cnd_phi_, "Cndphi" ); 

   iEvent.put( cnd_rp_, "Cndrp" );  
   iEvent.put( cnd_rpt_, "Cndrpt" ); 
   iEvent.put( cnd_rpx_, "Cndrpx" ); 
   iEvent.put( cnd_rpy_, "Cndrpy" ); 
   iEvent.put( cnd_rpz_, "Cndrpz" ); 
   iEvent.put( cnd_re_, "Cndre" );  
   iEvent.put( cnd_ret_, "Cndret" ); 
   iEvent.put( cnd_rq_, "Cndrq" );  
   iEvent.put( cnd_reta_, "Cndreta" ); 
   iEvent.put( cnd_rphi_, "Cndrphi" ); 


   return;
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
TagProbeEDMNtuple::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TagProbeEDMNtuple::endJob() {
}

// ****************************************************************** //
// ********************* Tree Fill Functions! *********************** //
// ****************************************************************** //


// ********************* Fill Run & Event Info *********************** //
void
TagProbeEDMNtuple::fillRunEventInfo()
{
   int run = m_event->id().run();
   int event = m_event->id().event();

   *run_ = run;
   *event_ = event;
}
// ****************************************************************** //

// ********************* Fill Trigger Info *********************** //
void
TagProbeEDMNtuple::fillTriggerInfo()
{
   
   // Trigger Info
   Handle<TriggerEventWithRefs> trgEvent;
   try{ m_event->getByLabel(triggerEventTag_,trgEvent); } 
   catch (...) 
   { 
      LogWarning("TagAndProbe") << "Could not extract trigger event summary "
				<< "with tag " << triggerEventTag_;
   }
   
   *nl1_ = 0;
   *nhlt_ = 0;
   if( trgEvent.isValid() )
   {
      vector< L1MuonParticleRef > muonL1CandRefs;
      size_type l1index = trgEvent->filterIndex(hltL1Tag_.label());
      if( l1index < trgEvent->size() )
      {
	 trgEvent->getObjects( l1index, trigger::TriggerL1Mu, muonL1CandRefs );
	 *nl1_ = muonL1CandRefs.size(); 
      }

      vector< RecoChargedCandidateRef > muonCandRefs;
      size_type index = trgEvent->filterIndex(hltTag_.label());
      if( index < trgEvent->size() )
      {
	 trgEvent->getObjects( index, trigger::TriggerMuon, muonCandRefs );
	 *nhlt_ = muonCandRefs.size(); 
      }
   }
}
// ****************************************************************** //



// ********************* Fill Tag-Probe Info *********************** //
void
TagProbeEDMNtuple::fillTagProbeInfo()
{

   // Trigger Info
   Handle<TriggerEventWithRefs> trgEvent;
   try{ m_event->getByLabel(triggerEventTag_,trgEvent); } 
   catch (...) 
   { 
      LogWarning("Z") << "Could not extract trigger event summary "
		      << "with tag " << triggerEventTag_;
   }

   // Fill some information about the tag & probe collections
   int nrtp = 0;
   for( int itype=0; itype<(int)tagProbeMapTags_.size(); ++itype )
   {
      // Tag-Probes
      Handle<CandCandAssociationCollection> tagprobes;
      try{ m_event->getByLabel(tagProbeMapTags_[itype],tagprobes); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract tag-probe map with input tag " 
			 << tagProbeMapTags_[itype];
      }

      // Tag MC Truth Match Maps
      Handle<CandMatchMap> tagmatch;
      try{ m_event->getByLabel(tagTruthMatchMapTags_[itype],tagmatch); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract muon tag match map "
			 << "with input tag " << tagTruthMatchMapTags_[itype];
      }

      // Truth match for probe
      Handle<CandMatchMap> allprobematch;
      try{ m_event->getByLabel(allProbeTruthMatchMapTags_[itype],allprobematch); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract muon allprobe match map "
			 << "with input tag " << allProbeTruthMatchMapTags_[itype];
      }

      // Passing Probe Candidates
      Handle<CandidateCollection> passprobemuons;
      try{ m_event->getByLabel(passProbeCandTags_[itype],passprobemuons); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract tag muon cands with input tag " 
			 << passProbeCandTags_[itype];
      }

      if( tagprobes.isValid() )
      {
	 CandCandAssociationCollection::const_iterator tpItr = tagprobes->begin();
	 for( ; tpItr != tagprobes->end(); ++tpItr )
	 {
	    const CandidateRef &tag = tpItr->key;
	    const CandidateRefVector &probes = tpItr->val;

	    // If there are two probes with the tag continue
	    if( probes.size() > 1 ) continue;

	    math::XYZTLorentzVector tpP4 = tag->p4() + probes[0]->p4();

	    // Is this Tag-Probe pair from a true Z?
	    // See if both the daughters are matched.
	    int tptrue = 0;
	    bool tagFromZ   = CandFromZ(tag,tagmatch);
	    bool probeFromZ = CandFromZ(probes[0],allprobematch);

	    // If both tag and probe are from Z .. set to true
	    if( tagFromZ && probeFromZ ) tptrue = 1;
	    cout << "Is true Z? " << tptrue << endl;

	    // Is this probe in the set that pass the efficiency criteria?
	    int ppass = ProbePassProbeOverlap(probes[0],passprobemuons);

	    // Did this tag cause a L1 and/or HLT trigger?
	    bool l1Trigger = false;
	    bool hltTrigger = false;
	    if( trgEvent.isValid() )
	    {
	       vector< L1MuonParticleRef > muonL1CandRefs;
	       size_type l1index = trgEvent->filterIndex(hltL1Tag_.label());
	       if( l1index < trgEvent->size() )
	       {
		  trgEvent->getObjects( l1index, trigger::TriggerL1Mu, muonL1CandRefs );
		  int npart = muonL1CandRefs.size();
		  for(int ipart = 0; ipart != npart; ++ipart)
		  { 
		     const L1MuonParticle* l1mu = (muonL1CandRefs[ipart]).get();
		     l1Trigger = MatchObjects( (const Candidate*)l1mu, tag );
		     if( l1Trigger ) break;
		  }
	       }

	       vector< RecoChargedCandidateRef > muonCandRefs;
	       size_type index = trgEvent->filterIndex(hltTag_.label());
	       if( index < trgEvent->size() )
	       {
		  trgEvent->getObjects( index, trigger::TriggerMuon, muonCandRefs );
		  int npart = muonCandRefs.size();
		  for(int ipart = 0; ipart != npart; ++ipart)
		  { 
		     const RecoChargedCandidate* muon = (muonCandRefs[ipart]).get();
		     hltTrigger = MatchObjects( (const Candidate*)muon, tag );
		     if( hltTrigger ) break;
		  }
	       }
	    }

	    double mass = tpP4.M();
	    double px   = tpP4.Px();
	    double py   = tpP4.Py();
	    double pz   = tpP4.Pz();
	    double pt   = tpP4.Pt();
	    double p    = tpP4.P();
	    double e    = tpP4.E();
	    double et   = tpP4.Et();

	    tp_type_->push_back(   itype );
	    tp_true_->push_back(   tptrue );
	    tp_ppass_->push_back(  ppass );
	    tp_l1_->push_back(  l1Trigger );
	    tp_hlt_->push_back(  hltTrigger );
	    tp_p_->push_back(   p );
	    tp_px_->push_back(  px );
	    tp_py_->push_back(  py );
	    tp_pz_->push_back(  pz );
	    tp_pt_->push_back(  pt );
	    tp_e_->push_back(   e );
	    tp_et_->push_back(  et );
	    tp_mass_->push_back(  mass );

	    double dp   = tag->p();
	    double dpx  = tag->px();
	    double dpy  = tag->py();
	    double dpz  = tag->pz();
	    double dpt  = tag->pt();
	    double de   = tag->energy();
	    double det  = tag->et();
	    double dq   = tag->charge();
	    double deta = tag->eta();
	    double dphi = tag->phi();

	    tp_tag_p_->push_back(    dp );
	    tp_tag_px_->push_back(   dpx );
	    tp_tag_py_->push_back(   dpy );
	    tp_tag_pz_->push_back(   dpz );
	    tp_tag_pt_->push_back(   dpt );
	    tp_tag_e_->push_back(    de );
	    tp_tag_et_->push_back(   det );
	    tp_tag_q_->push_back(    dq );
	    tp_tag_eta_->push_back(  deta );
	    tp_tag_phi_->push_back(  dphi );

	    dp   = probes[0]->p();
	    dpx  = probes[0]->px();
	    dpy  = probes[0]->py();
	    dpz  = probes[0]->pz();
	    dpt  = probes[0]->pt();
	    de   = probes[0]->energy();
	    det  = probes[0]->et();
	    dq   = probes[0]->charge();
	    deta = probes[0]->eta();
	    dphi = probes[0]->phi();

	    tp_probe_p_->push_back(    dp );
	    tp_probe_px_->push_back(   dpx );
	    tp_probe_py_->push_back(   dpy );
	    tp_probe_pz_->push_back(   dpz );
	    tp_probe_pt_->push_back(   dpt );
	    tp_probe_e_->push_back(    de );
	    tp_probe_et_->push_back(   det );
	    tp_probe_q_->push_back(    dq );
	    tp_probe_eta_->push_back(  deta );
	    tp_probe_phi_->push_back(  dphi );

	    ++nrtp;
	 }
      }
   }
   nrtp_.reset( new int(nrtp) );
}
// ******************************************************************* //


// ********************* Fill True Efficiency Info *********************** //
void
TagProbeEDMNtuple::fillTrueEffInfo()
{

   // Should change this to get the eff info for all types of tag-probe!!
   Handle<CandidateCollection> genparticles;
   try{ m_event->getByLabel(genParticlesTag_,genparticles); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract gen particles with input tag " 
			<< genParticlesTag_;
   }

   int ncnd = 0;
   // Fill some information about the muon efficiency
   if( genparticles.isValid() )
   {
      for( int itype=0; itype<(int)tagProbeMapTags_.size(); ++itype )
      {
	 Handle<CandMatchMap> tagmatch;
	 try{ m_event->getByLabel(tagTruthMatchMapTags_[itype],tagmatch); }
	 catch(...)
	 { 
	    LogWarning("Z") << "Could not extract muon match map "
			    << "with input tag " << tagTruthMatchMapTags_[itype];
	 }

	 Handle<CandMatchMap> apmatch;
	 try{ m_event->getByLabel(allProbeTruthMatchMapTags_[itype],apmatch); }
	 catch(...)
	 { 
	    LogWarning("Z") << "Could not extract all probe match map "
			    << "with input tag " << allProbeTruthMatchMapTags_[itype];
	 }

	 Handle<CandMatchMap> ppmatch;
	 try{ m_event->getByLabel(passProbeTruthMatchMapTags_[itype],ppmatch); }
	 catch(...)
	 { 
	    LogWarning("Z") << "Could not extract pass probe match map "
			    << "with input tag " << passProbeTruthMatchMapTags_[itype];
	 }
   
	 for( unsigned int i=0; i<genparticles->size(); i++ )
	 {
	    int pdg_id = (*genparticles)[i].pdgId();

	    // If this is not a muon keep going!
	    if( abs( pdg_id ) != candPDGId_ ) continue;

	    int moid  = -1;
	    int gmoid = -1;
	    const Candidate *mcand = (*genparticles)[i].mother();
	    if( mcand != 0 )
	    {
	       moid = mcand->pdgId();
	       if( mcand->pdgId() == pdg_id )
	       {
		  if( mcand->mother() != 0 )
		  {
		     const Candidate *gcand = mcand->mother();
		     gmoid = gcand->pdgId();
		  }
	       }
	    }

	    int ftag = 0;
	    int fapb = 0;
	    int fppb = 0;

	    double p   = (*genparticles)[i].p();
	    double px = (*genparticles)[i].px();
	    double py = (*genparticles)[i].py();
	    double pz = (*genparticles)[i].pz();
	    double pt = (*genparticles)[i].pt();
	    double e  = (*genparticles)[i].energy();
	    double et  = (*genparticles)[i].et();
	    double q  = (*genparticles)[i].charge();
	    double phi = (*genparticles)[i].phi();
	    double eta = (*genparticles)[i].eta();

	    double rp   = 0;
	    double rpx  = 0;
	    double rpy  = 0;
	    double rpz  = 0; 
	    double rpt  = 0; 
	    double re   = 0;
	    double ret  = 0;
	    double rq   = 0;
	    double rphi = 0;
	    double reta = 0;

	    CandidateRef mcRef( genparticles, (size_t)i );
	    if( tagmatch.isValid() )
	    {
	       CandMatchMap::const_iterator f = tagmatch->begin();
	       for( ; f != tagmatch->end(); ++f )
	       {
		  const Candidate *mcMatchRef  = &*(f->val);

		  if( &(*mcRef)==&(*mcMatchRef) ) 
		  {
		     ftag = 1;
		     const Candidate *cnd = &*(f->key);
		     rp   = cnd->p();
		     rpx  = cnd->px();
		     rpy  = cnd->py();
		     rpz  = cnd->pz(); 
		     rpt  = cnd->pt(); 
		     re   = cnd->energy();
		     ret  = cnd->et();
		     rq   = cnd->charge();
		     rphi = cnd->phi();
		     reta = cnd->eta();
		  }
	       }
	    }
	    if( apmatch.isValid() )
	    {
	       CandMatchMap::const_iterator t = apmatch->begin();
	       for( ; t != apmatch->end(); ++t )
	       {
		  const Candidate *mcMatchRef  = &*(t->val);
		  if( &(*mcRef)==&(*mcMatchRef) )
		  {
		     fapb = 1;
		     if( ftag == 0 )
		     {
			const Candidate *tmu = &*(t->key);
			rp   = tmu->p();
			rpx  = tmu->px();
			rpy  = tmu->py();
			rpz  = tmu->pz(); 
			rpt  = tmu->pt(); 
			re   = tmu->energy();
			ret  = tmu->et();
			rq   = tmu->charge();
			rphi = tmu->phi();
			reta = tmu->eta();
		     }
		  }
	       }
	    }
	    if( ppmatch.isValid() )
	    {
	       CandMatchMap::const_iterator t = ppmatch->begin();
	       for( ; t != ppmatch->end(); ++t )
	       {
		  const Candidate *mcMatchRef  = &*(t->val);
		  if( &(*mcRef)==&(*mcMatchRef) )
		  {
		     fppb = 1;
		     if( ftag == 0 && fapb == 0 )
		     {
			const Candidate *tmu = &*(t->key);
			rp   = tmu->p();
			rpx  = tmu->px();
			rpy  = tmu->py();
			rpz  = tmu->pz(); 
			rpt  = tmu->pt(); 
			re   = tmu->energy();
			ret  = tmu->et();
			rq   = tmu->charge();
			rphi = tmu->phi();
			reta = tmu->eta();
		     }
		  }
	       }
	    }

	    cnd_type_->push_back(     itype );

	    cnd_tag_->push_back(     ftag );
	    cnd_aprobe_->push_back(  fapb );
	    cnd_pprobe_->push_back(  fppb );
	    cnd_moid_->push_back(  moid );
	    cnd_gmid_->push_back(  gmoid );

	    cnd_p_->push_back(  p );
	    cnd_px_->push_back(  px );
	    cnd_py_->push_back(  py );
	    cnd_pz_->push_back(  pz );
	    cnd_pt_->push_back(  pt );
	    cnd_e_->push_back(  e );
	    cnd_et_->push_back(  et );
	    cnd_q_->push_back(  q );
	    cnd_phi_->push_back(  phi );
	    cnd_eta_->push_back(  eta );

	    cnd_rp_->push_back(  rp );
	    cnd_rpx_->push_back(  rpx );
	    cnd_rpy_->push_back(  rpy );
	    cnd_rpz_->push_back(  rpz );
	    cnd_rpt_->push_back(  rpt );
	    cnd_re_->push_back(  re );
	    cnd_ret_->push_back(  ret );
	    cnd_rq_->push_back(  rq );
	    cnd_rphi_->push_back(  rphi );
	    cnd_reta_->push_back(  reta );

	    ncnd++;
	 }
      }
   }
   ncnd_.reset( new int(ncnd) );

}
// *********************************************************************** //

// ***************** Are Candidates from a Z? ******************** //
bool TagProbeEDMNtuple::CandFromZ( const reco::CandidateRef& cand, 
				   edm::Handle<reco::CandMatchMap>& matches )
{
   bool isFromZ = false;

   if( matches.isValid() )
   {
      CandMatchMap::const_iterator f = matches->begin();
      for( ; f != matches->end(); ++f )
      {
	 const Candidate *matchRef = &*(f->key);
	 if( &(*cand)==&(*matchRef) ) 
	 {
	    //const Candidate *mcMatchRef  = &*(f->val);
	    CandidateRef mcMatchRef  = f->val;
		  
	    if( mcMatchRef.isNonnull() )
	    {
	       //cout << "Id " << mcMatchRef->pdgId() << endl;
	       int moid = -99;
	       int gmoid = -99;
	       const Candidate *mcand = mcMatchRef->mother();
	       if( mcand != 0 )
	       {
		  moid = mcand->pdgId();
		  if( mcand->pdgId() == mcMatchRef->pdgId() )
		  {
		     if( mcand->mother() != 0 )
		     {
			const Candidate *gcand = mcand->mother();
			gmoid = gcand->pdgId();
		     }
		  }
	       }
	       if( moid == 23 || gmoid == 23 ) isFromZ = true;
	    }
	 }
      }
   }

   return isFromZ;
}
// *************************************************************** //

// ***************** Does this probe pass? ******************** //
int TagProbeEDMNtuple::ProbePassProbeOverlap( const CandidateRef& probe, 
					      Handle<CandidateCollection>& passprobes )
{
   int ppass = 0;
   if( candType_ == "Muon" )
   {
      MuonRef probeRef;
      if( probe->hasMasterClone() )
      {
	 probeRef = (probe->masterClone()).castTo<MuonRef>();
      }
      if( passprobes.isValid() && probeRef.isNonnull() )
      {
	 for( int ipp=0; ipp<(int)passprobes->size(); ++ipp )
	 {
	    bool isOverlap = false;

	    MuonRef ppRef;
	    if( (*passprobes)[ipp].hasMasterClone() )
	    {
	       ppRef = ((*passprobes)[ipp].masterClone()).castTo<MuonRef>();
	    }
	    if( ppRef.isNonnull() )
	    {
	       if( (ppRef->track()).isNonnull() && 
		   (probeRef->track()).isNonnull() && 
		   ppRef->track() == probeRef->track() )
		  isOverlap = true;
	       if( (ppRef->standAloneMuon()).isNonnull() && 
		   (probeRef->standAloneMuon()).isNonnull() && 
		   ppRef->standAloneMuon() == probeRef->standAloneMuon() )
		  isOverlap = true;
	       if( (ppRef->combinedMuon()).isNonnull() && 
		   (probeRef->combinedMuon()).isNonnull() && 
		   ppRef->combinedMuon() == probeRef->combinedMuon() )
		  isOverlap = true;
	       if( isOverlap ) ppass = 1;
	    }
	 } 
      }
   }
   else
   {
      if( passprobes.isValid() )
      {
	 for( int ipp=0; ipp<(int)passprobes->size(); ++ipp )
	 {
	    bool isOverlap = MatchObjects(&((*passprobes)[ipp]),probe);
	    if( isOverlap ) ppass = 1;
	 } 
      }      
   }

   return ppass;
}
// ************************************************************ //


// ***************** Trigger object matching ******************** //
bool TagProbeEDMNtuple::MatchObjects( const Candidate *hltObj, 
				      const CandidateRef& tagObj )
{
   double tEta = tagObj->eta();
   double tPhi = tagObj->phi();
   double hEta = hltObj->eta();
   double hPhi = hltObj->phi();

   double dRval = deltaR(tEta, tPhi, hEta, hPhi);
   return ( dRval < delRMatchingCut_ );
}
// ************************************************************** //

//define this as a plug-in
//DEFINE_FWK_MODULE(TagProbeEDMNtuple);
