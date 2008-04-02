// -*- C++ -*-
//
// Package:    TagAndProbe
// Class:      MuonTagProbeAnalyzer
// 
/**\class TagAndProbe MuonTagProbeAnalyzer.cc MuonAnalysis/TagAndProbe/src/MuonTagProbeAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nadia Adam
//         Created:  Tue Mar  4 12:27:53 CST 2008
// $Id$
//
//


// system include files
#include <cmath>

// user include files
#include "MuonAnalysis/TagAndProbe/interface/MuonTagProbeAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "AnalysisDataFormats/Muon/interface/MuonTagProbeAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
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
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/FitResult.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
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

#define DEBUGLVL 0

//
// class decleration
//

//
// constants, enums and typedefs
//
using namespace std; 
using namespace reco; 
using namespace edm;

//
// static data member definitions
//

//
// constructors and destructor           
//
MuonTagProbeAnalyzer::MuonTagProbeAnalyzer(const ParameterSet& iConfig)
{

   rootFile_ = iConfig.getUntrackedParameter<string>("rootFile");

   //now do what ever initialization is needed
   m_fillEvt.init(rootFile_, &m_evtTree);

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
   const edm::InputTag dTrackTag("ctfWithMaterialTracks");
   tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tracksTag",dTrackTag);

   const edm::InputTag dsaTrackTag("standAloneMuons");
   satracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>(
      "standAloneMuonTracksTag",dsaTrackTag);

   const edm::InputTag dcTrackTag("globalMuons");
   ctracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>(
      "globalMuonTracksTag",dcTrackTag);
   // ********************************* //

   // ********** Muons ********** //
   const edm::InputTag dMuonTag("muons");
   muonsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("muonsTag",dMuonTag);
   // *************************** //

   // ********** Tag-Probes ********** //
   vector< edm::InputTag > defaultTagProbeMapTags;
   tagProbeMapTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "tagProbeMapTags",defaultTagProbeMapTags);
   for( int i=0; i<(int)tagProbeMapTags_.size(); ++i )
      cout << tagProbeMapTags_[i] << endl;
   // ******************************** //

   int map_size = (int)tagProbeMapTags_.size();
   const edm::InputTag dBlankTag("Blank");

   // ********** Tag & Probe Muons *********** //
   vector< edm::InputTag > defaultTagMuonTags;
   tagMuonTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "tagMuonTags",defaultTagMuonTags);
   for( int i=0; i<(int)tagMuonTags_.size(); ++i ) 
      cout << tagMuonTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)tagMuonTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)tagMuonTags_.size()); ++i ) 
	 tagMuonTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultAllProbeMuonTags;
   allProbeMuonTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "allProbeMuonTags",defaultAllProbeMuonTags);
   for( int i=0; i<(int)allProbeMuonTags_.size(); ++i ) 
      cout << allProbeMuonTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)allProbeMuonTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)allProbeMuonTags_.size()); ++i ) 
	 allProbeMuonTags_.push_back(dBlankTag);
   } 

   vector< edm::InputTag > defaultPassProbeMuonTags;
   passProbeMuonTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "passProbeMuonTags",defaultPassProbeMuonTags);
   for( int i=0; i<(int)passProbeMuonTags_.size(); ++i ) 
      cout << passProbeMuonTags_[i] << endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)passProbeMuonTags_.size() < map_size )
   {
      cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << endl;
      for( int i=0; i<(map_size-(int)passProbeMuonTags_.size()); ++i ) 
	 passProbeMuonTags_.push_back(dBlankTag);
   } 
   // **************************************** //

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

   // ********** Z Cands ********** //
   vector< edm::InputTag > defaultZCandTags;
   zCandTags_ = iConfig.getUntrackedParameter< vector<edm::InputTag> >(
      "zCandTags",defaultZCandTags);
   for( int i=0; i<(int)zCandTags_.size(); ++i ) 
      cout << zCandTags_[i] << endl;
   // ***************************** //


   // HLT
   const edm::InputTag dHLTL1Tag("SingleMuIsoLevel1Seed");
   hltL1Tag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltL1Tag",dHLTL1Tag);

   const edm::InputTag dHLTTag("SingleMuIsoL3IsoFiltered");
   hltTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltTag",dHLTTag);

   delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch",0.15);

   // Primary Vertex
   const edm::InputTag dPrimaryVertexTag("offlinePrimaryVerticesFromCTFTracks");
   PvxtTag_  = iConfig.getUntrackedParameter<edm::InputTag>(
      "verticesTag",dPrimaryVertexTag);

   isMC_ = iConfig.getUntrackedParameter< bool >("isMC",true);

}


MuonTagProbeAnalyzer::~MuonTagProbeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MuonTagProbeAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
   cout << "Here in muon efficiency analyzer ..." << endl;

   // Set the private event & setup pointers
   m_event = &iEvent;
   m_setup = &iSetup;

   // Initialize the Tree values
   m_evtTree.init();

   // Fill the run and event number information
   fillRunEventInfo();

   // Fill event level trigger info
   fillTriggerInfo();

   // Fill Generator Level Information
   fillMCInfo();

   // Fill Track Info
   fillTrackInfo();

   // Fill Muon Info
   fillMuonInfo();

   // Fill Tag-Probe Info
   fillTagProbeInfo();

   // Fill Efficiency Info for true muons
   fillMuonEffInfo();

   // Fill Z Cand Info
   fillZInfo();

   // Fill Vertex Info
   fillVertexInfo();

   // Finally fill the tree for this event!
   m_fillEvt.fill();   

   return;
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonTagProbeAnalyzer::beginJob(const EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonTagProbeAnalyzer::endJob() {

   m_fillEvt.finalize();
}

// ****************************************************************** //
// ********************* Tree Fill Functions! *********************** //
// ****************************************************************** //


// ********************* Fill Run & Event Info *********************** //
void
MuonTagProbeAnalyzer::fillRunEventInfo()
{
   int run = m_event->id().run();
   int event = m_event->id().event();

   m_evtTree.run = run;
   m_evtTree.event = event;
   m_evtTree.xs = 0;
   m_evtTree.num = 0;
}
// ****************************************************************** //

// ********************* Fill Trigger Info *********************** //
void
MuonTagProbeAnalyzer::fillTriggerInfo()
{
   
   // L1 Info
   Handle<HLTFilterObjectWithRefs> l1MuonTrg;
   try{ m_event->getByLabel(hltL1Tag_,l1MuonTrg); } 
   catch (...) 
   { 
      LogWarning("Z") << "Could not extract HLT L1 object with tag "
		      << hltL1Tag_;
   }

   // HLT Info
   Handle<HLTFilterObjectWithRefs> hltMuonTrg;
   try{ m_event->getByLabel(hltTag_,hltMuonTrg); } 
   catch (...) 
   { 
      LogWarning("Z") << "Could not extract HLT object with tag "
		      << hltTag_;
   }

   m_evtTree.nl1 = 0;
   if( l1MuonTrg.isValid() ) m_evtTree.nl1 = l1MuonTrg->size(); 

   m_evtTree.nhlt = 0;
   if( hltMuonTrg.isValid() ) m_evtTree.nhlt = hltMuonTrg->size(); 
}
// ****************************************************************** //

// ********************* Fill MC Info *********************** //
void
MuonTagProbeAnalyzer::fillMCInfo()
{
   // Extract the MC information from the frame, if we are running on MC
   if( isMC_ )
   {      
      // Set the number of different particle types we're storing
      int nmcpart = mcParticles_.size();
      m_evtTree.nmcpart = nmcpart;

      vector< Handle< HepMCProduct > > EvtHandles ;
      m_event->getManyByType( EvtHandles ) ;
   
      // Loop over all MC products
      for ( unsigned int i=0; (i<EvtHandles.size() && nmcpart>0); i++ )
      {
	 if( i>0 ) continue;

	 if ( EvtHandles[i].isValid() )
	 {
	    // Get the event
	    const HepMC::GenEvent* Evt = EvtHandles[i]->GetEvent() ;
   
	    // Loop over particles and extract any that the user has asked to
	    // be stored
	    int *nmc = new int[nmcpart];
	    for( int i=0; i<nmcpart; i++ ) nmc[i] = 0;
	    HepMC::GenEvent::particle_const_iterator Part = Evt->particles_begin();
	    HepMC::GenEvent::particle_const_iterator PartEnd = Evt->particles_end();
	    for(; Part != PartEnd; Part++ )
	    {
	       int pdg = (*Part)->pdg_id();
	       
	       for( int j=0; j<nmcpart; j++ )
	       {
		  if( abs(pdg) == mcParticles_[j] && nmc[j] < EvtTree::MAXNMCT)
		  {
		     HepMC::FourVector p4 = (*Part)->momentum();

		     int    pid = pdg;
		     int    barcode = (*Part)->barcode();
		     double p   = p4.mag();
		     double px = p4.px();
		     double py = p4.py();
		     double pz = p4.pz();
		     double pt = p4.perp();
		     double e  = p4.e();
		     double phi = p4.phi() ;
		     double eta = -log(tan(p4.theta()/2.));
                     double mass = e*e - p*p;
                     if( mass > 0 ) mass = sqrt(mass);
                     else           mass = -sqrt(fabs(mass));

		     // Get the parent barcode (in general there
		     // can be more than one parent, but for now
		     // we will just take the first).
		     int parent_pi = 0;
		     int parent_bc = 0;
		     HepMC::GenVertex *pvtx = (*Part)->production_vertex();

		     // Loop over the particles in this vertex
		     if( pvtx != NULL )
		     {
		        if( (*pvtx).particles_in_size() > 0 )
			{
			   HepMC::GenVertex::particles_in_const_iterator pIt = 
			      (*pvtx).particles_in_const_begin();
			   parent_pi = (*pIt)->pdg_id(); 
			   parent_bc = (*pIt)->barcode();
			}
		     }

		     int num_lep = 0;
		     if( abs(pdg) == 23 )
		     {
			HepMC::GenVertex *evtx = (*Part)->end_vertex();
			
			// Loop over the particles in this vertex
			if( evtx != NULL )
			{
			   if( (*evtx).particles_out_size() > 0 )
			   {
			      HepMC::GenVertex::particles_out_const_iterator pIt = 
				 (*evtx).particles_out_const_begin();
			      for( ; pIt != (*evtx).particles_out_const_end(); pIt++ )
			      {
				 if( abs((*pIt)->pdg_id()) >= 11 && 
				     abs((*pIt)->pdg_id()) <= 14 ) 
				 {
				    ++num_lep;
				 }
			      }
			   }
			}
		     }
		     if( abs(mcParticles_[j]) == 23 && num_lep == 0 ) continue;

		     // Check for the parents if required
		     if( mcParents_[j] != 0 && abs(parent_pi) != mcParents_[j] ) continue;

		     // set the tree values
		     m_evtTree.mc_pbc[j][nmc[j]] = parent_bc;
		     m_evtTree.mc_ppid[j][nmc[j]] = parent_pi;
		     m_evtTree.mc_bc[j][nmc[j]] = barcode;
		     m_evtTree.mc_pid[j][nmc[j]] = pid;
		     m_evtTree.mc_p[j][nmc[j]] = p;
		     m_evtTree.mc_mass[j][nmc[j]] = mass;
		     m_evtTree.mc_px[j][nmc[j]] = px;
		     m_evtTree.mc_py[j][nmc[j]] = py;
		     m_evtTree.mc_pz[j][nmc[j]] = pz;
		     m_evtTree.mc_pt[j][nmc[j]] = pt;
		     m_evtTree.mc_e[j][nmc[j]] = e;
		     m_evtTree.mc_phi[j][nmc[j]] = phi;
		     m_evtTree.mc_eta[j][nmc[j]] = eta;

		     nmc[j]++;
		  }
	       }
	    }

	    for( int k=0; k<nmcpart; k++ ) m_evtTree.nmc[k] = nmc[k];
	 }
      }
   }
}
// ********************************************************** //

// ********************* Fill Track Info *********************** //
void
MuonTagProbeAnalyzer::fillTrackInfo()
{

   // Fill info for all tracks
   Handle<TrackCollection> tracks;
   try{ m_event->getByLabel(tracksTag_,tracks); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract reco tracks with input tag " 
			<< tracksTag_;
   }

   // Fill some information about the reconstructed tracks
   int nrtrk = 0;
   if( tracks.isValid() )
   {
      for( unsigned int i=0; i<tracks->size(); i++ )
      {
	 if( nrtrk >= EvtTree::MAXNTRK ) break;

	 int id = i+1;
	 int type = 0;
	 double chi2 = (*tracks)[i].chi2();
	 double ndof = (*tracks)[i].ndof();
	 double p   = (*tracks)[i].p();
	 double px = (*tracks)[i].px();
	 double py = (*tracks)[i].py();
	 double pz = (*tracks)[i].pz();
	 double pt = (*tracks)[i].pt();
	 double e  = (*tracks)[i].p();
	 double q  = (*tracks)[i].charge();
	 double phi = (*tracks)[i].phi();
	 double eta = (*tracks)[i].eta();
	 double dxy = (*tracks)[i].dxy();
	 double d0  = (*tracks)[i].d0();
	 double dsz = (*tracks)[i].dsz();
	 double dz  = (*tracks)[i].dz();
	 double vx  = (*tracks)[i].vx();
	 double vy  = (*tracks)[i].vy();
	 double vz  = (*tracks)[i].vz();

	 m_evtTree.trk_id[nrtrk] = id;
	 m_evtTree.trk_type[nrtrk] = type;
	 m_evtTree.trk_chi2[nrtrk] = chi2;
	 m_evtTree.trk_ndof[nrtrk] = ndof;
	 m_evtTree.trk_p[nrtrk] = p;
	 m_evtTree.trk_px[nrtrk] = px;
	 m_evtTree.trk_py[nrtrk] = py;
	 m_evtTree.trk_pz[nrtrk] = pz;
	 m_evtTree.trk_pt[nrtrk] = pt;
	 m_evtTree.trk_e[nrtrk] = e;
	 m_evtTree.trk_q[nrtrk] = q;
	 m_evtTree.trk_phi[nrtrk] = phi;
	 m_evtTree.trk_eta[nrtrk] = eta;
	 m_evtTree.trk_dxy[nrtrk] = dxy;
	 m_evtTree.trk_d0[nrtrk] = d0;
	 m_evtTree.trk_dsz[nrtrk] = dsz;
	 m_evtTree.trk_dz[nrtrk] = dz;
	 m_evtTree.trk_vx[nrtrk] = vx;
	 m_evtTree.trk_vy[nrtrk] = vy;
	 m_evtTree.trk_vz[nrtrk] = vz;

	 nrtrk++;
      }
   }


   // Fill info for all stand alone muon tracks
   Handle<TrackCollection> satracks;
   try{ m_event->getByLabel(satracksTag_,satracks); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract reco satracks with input tag " 
			<< satracksTag_;
   }
   

   // Fill some information about the reconstructed satracks
   if( satracks.isValid() )
   {
      for( unsigned int i=0; i<satracks->size(); i++ )
      {
	 if( nrtrk >= EvtTree::MAXNTRK ) break;

	 int id = i+1;
	 int type = 1;
	 double chi2 = (*satracks)[i].chi2();
	 double ndof = (*satracks)[i].ndof();
	 double p   = (*satracks)[i].p();
	 double px = (*satracks)[i].px();
	 double py = (*satracks)[i].py();
	 double pz = (*satracks)[i].pz();
	 double pt = (*satracks)[i].pt();
	 double e  = (*satracks)[i].p();
	 double q  = (*satracks)[i].charge();
	 double phi = (*satracks)[i].phi();
	 double eta = (*satracks)[i].eta();
	 double dxy = (*satracks)[i].dxy();
	 double d0  = (*satracks)[i].d0();
	 double dsz = (*satracks)[i].dsz();
	 double dz  = (*satracks)[i].dz();
	 double vx  = (*satracks)[i].vx();
	 double vy  = (*satracks)[i].vy();
	 double vz  = (*satracks)[i].vz();

	 m_evtTree.trk_id[nrtrk] = id;
	 m_evtTree.trk_type[nrtrk] = type;
	 m_evtTree.trk_chi2[nrtrk] = chi2;
	 m_evtTree.trk_ndof[nrtrk] = ndof;
	 m_evtTree.trk_p[nrtrk] = p;
	 m_evtTree.trk_px[nrtrk] = px;
	 m_evtTree.trk_py[nrtrk] = py;
	 m_evtTree.trk_pz[nrtrk] = pz;
	 m_evtTree.trk_pt[nrtrk] = pt;
	 m_evtTree.trk_e[nrtrk] = e;
	 m_evtTree.trk_q[nrtrk] = q;
	 m_evtTree.trk_phi[nrtrk] = phi;
	 m_evtTree.trk_eta[nrtrk] = eta;
	 m_evtTree.trk_dxy[nrtrk] = dxy;
	 m_evtTree.trk_d0[nrtrk] = d0;
	 m_evtTree.trk_dsz[nrtrk] = dsz;
	 m_evtTree.trk_dz[nrtrk] = dz;
	 m_evtTree.trk_vx[nrtrk] = vx;
	 m_evtTree.trk_vy[nrtrk] = vy;
	 m_evtTree.trk_vz[nrtrk] = vz;

	 nrtrk++;
      }
   }


   // Fill info for global muon tracks
   Handle<TrackCollection> ctracks;
   try{ m_event->getByLabel(ctracksTag_,ctracks); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract reco ctracks with input tag " 
			<< ctracksTag_;
   }

   // Fill some information about the reconstructed ctracks
   if( ctracks.isValid() )
   {

      for( unsigned int i=0; i<ctracks->size(); i++ )
      {
	 if( nrtrk >= EvtTree::MAXNTRK ) break;

	 int id = i+1;
	 int type = 2;
	 double chi2 = (*ctracks)[i].chi2();
	 double ndof = (*ctracks)[i].ndof();
	 double p   = (*ctracks)[i].p();
	 double px = (*ctracks)[i].px();
	 double py = (*ctracks)[i].py();
	 double pz = (*ctracks)[i].pz();
	 double pt = (*ctracks)[i].pt();
	 double e  = (*ctracks)[i].p();
	 double q  = (*ctracks)[i].charge();
	 double phi = (*ctracks)[i].phi();
	 double eta = (*ctracks)[i].eta();
	 double dxy = (*ctracks)[i].dxy();
	 double d0  = (*ctracks)[i].d0();
	 double dsz = (*ctracks)[i].dsz();
	 double dz  = (*ctracks)[i].dz();
	 double vx  = (*ctracks)[i].vx();
	 double vy  = (*ctracks)[i].vy();
	 double vz  = (*ctracks)[i].vz();

	 m_evtTree.trk_id[nrtrk] = id;
	 m_evtTree.trk_type[nrtrk] = type;
	 m_evtTree.trk_chi2[nrtrk] = chi2;
	 m_evtTree.trk_ndof[nrtrk] = ndof;
	 m_evtTree.trk_p[nrtrk] = p;
	 m_evtTree.trk_px[nrtrk] = px;
	 m_evtTree.trk_py[nrtrk] = py;
	 m_evtTree.trk_pz[nrtrk] = pz;
	 m_evtTree.trk_pt[nrtrk] = pt;
	 m_evtTree.trk_e[nrtrk] = e;
	 m_evtTree.trk_q[nrtrk] = q;
	 m_evtTree.trk_phi[nrtrk] = phi;
	 m_evtTree.trk_eta[nrtrk] = eta;
	 m_evtTree.trk_dxy[nrtrk] = dxy;
	 m_evtTree.trk_d0[nrtrk] = d0;
	 m_evtTree.trk_dsz[nrtrk] = dsz;
	 m_evtTree.trk_dz[nrtrk] = dz;
	 m_evtTree.trk_vx[nrtrk] = vx;
	 m_evtTree.trk_vy[nrtrk] = vy;
	 m_evtTree.trk_vz[nrtrk] = vz;

	 nrtrk++;

      }
   }

   m_evtTree.nrtrk = nrtrk;
}
// ************************************************************* //



// ********************* Fill Muon Object Info *********************** //
void
MuonTagProbeAnalyzer::fillMuonInfo()
{

   Handle<MuonCollection> muons;
   try{ m_event->getByLabel(muonsTag_,muons); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract reco muons with input tag " 
			<< muonsTag_;
   }

   // Fill some information about the reconstructed muons
   int nrmu = 0;
   if( muons.isValid() )
   {
      //cout << "There are " << muons->size() << " global muons" << endl;
      for( unsigned int i=0; i<muons->size(); i++ )
      {
	 if( nrmu >= EvtTree::MAXNMU ) break;

	 int id = i+1;
	 int gbl = (int)(*muons)[i].isGlobalMuon();
	 int sta = (int)(*muons)[i].isStandAloneMuon();
	 int trk = (int)(*muons)[i].isTrackerMuon();
	 double p   = (*muons)[i].p();
	 double px = (*muons)[i].px();
	 double py = (*muons)[i].py();
	 double pz = (*muons)[i].pz();
	 double pt = (*muons)[i].pt();
	 double e  = (*muons)[i].energy();
	 double et  = (*muons)[i].et();
	 double q  = (*muons)[i].charge();
	 double phi = (*muons)[i].phi();
	 double eta = (*muons)[i].eta();

	 m_evtTree.mu_id[nrmu] = id;
	 m_evtTree.mu_gbl[nrmu] = gbl;
	 m_evtTree.mu_sta[nrmu] = sta;
	 m_evtTree.mu_trk[nrmu] = trk;
	 m_evtTree.mu_p[nrmu] = p;
	 m_evtTree.mu_px[nrmu] = px;
	 m_evtTree.mu_py[nrmu] = py;
	 m_evtTree.mu_pz[nrmu] = pz;
	 m_evtTree.mu_pt[nrmu] = pt;
	 m_evtTree.mu_e[nrmu] = e;
	 m_evtTree.mu_et[nrmu] = et;
	 m_evtTree.mu_q[nrmu] = q;
	 m_evtTree.mu_phi[nrmu] = phi;
	 m_evtTree.mu_eta[nrmu] = eta;

	 nrmu++;
      }
   }

   m_evtTree.nrmu = nrmu;
}
// ******************************************************************* //

// ********************* Fill Tag-Probe Info *********************** //
void
MuonTagProbeAnalyzer::fillTagProbeInfo()
{

   // L1 Info
   Handle<HLTFilterObjectWithRefs> l1MuonTrg;
   try{ m_event->getByLabel(hltL1Tag_,l1MuonTrg); } 
   catch (...) 
   { 
      LogWarning("Z") << "Could not extract HLT L1 object with tag "
		      << hltL1Tag_;
   }

   // HLT Info
   Handle<HLTFilterObjectWithRefs> hltMuonTrg;
   try{ m_event->getByLabel(hltTag_,hltMuonTrg); } 
   catch (...) 
   { 
      LogWarning("Z") << "Could not extract HLT object with tag "
		      << hltTag_;
   }

   // Fill some information about the tag & probe collections
   int nrtp = 0;
   for( int itype=0; itype<(int)tagProbeMapTags_.size(); ++itype )
   {
      // Tag-Probes
      Handle<MuonTagProbeAssociationCollection> tagprobes;
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

      // Tag Candidate collections
      Handle<CandidateCollection> taggedmuons;
      try{ m_event->getByLabel(tagCandTags_[itype],taggedmuons); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract tag muon cands with input tag " 
			 << tagCandTags_[itype];
      }
      for( int i=0; i<(int)taggedmuons->size(); ++i )
      {
	 const Track * trk = (*taggedmuons)[i].get<const Track *>();
	 if( trk == 0 ) cout << "No track reference!" << endl;
	 cout << "Has master clone? " << (*taggedmuons)[i].hasMasterClone() << endl;
	 if( (*taggedmuons)[i].hasMasterClone() )
	 {
	    MuonRef muRef = ((*taggedmuons)[i].masterClone()).castTo<MuonRef>();
	    cout << "Muon pt " << muRef->pt() << endl;
	 }
      }


      // Truth match for probe
      Handle<CandMatchMap> allprobematch;
      try{ m_event->getByLabel(allProbeTruthMatchMapTags_[itype],allprobematch); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract muon allprobe match map "
			 << "with input tag " << allProbeTruthMatchMapTags_[itype];
      }

      // Probe Candidates
      Handle<CandidateCollection> allprobemuons;
      try{ m_event->getByLabel(allProbeCandTags_[itype],allprobemuons); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract tag muon cands with input tag " 
			 << allProbeCandTags_[itype];
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
	 MuonTagProbeAssociationCollection::const_iterator tpItr = tagprobes->begin();
	 for( ; tpItr != tagprobes->end(); ++tpItr )
	 {
	    if( nrtp >= EvtTree::MAXNPAIR ) break;

	    const MuonRef &tag = tpItr->key;
	    const MuonRefVector &probes = tpItr->val;
	    //cout << "Tag key " << tag.key() << endl;

	    // If there are two probes with the tag continue
	    if( probes.size() > 1 ) continue;

	    math::XYZTLorentzVector tpP4 = tag->p4() + probes[0]->p4();

	    // Is this Tag-Probe pair from a true Z?
	    // See if both the daughters are matched.
	    int tptrue = 0;

	    bool tagFromZ = false;
	    if( taggedmuons.isValid() )
	    {
	       CandidateRef tagRef( taggedmuons, tag.key() );
	       if( tagmatch.isValid() )
	       {
		  CandMatchMap::const_iterator f = tagmatch->begin();
		  for( ; f != tagmatch->end(); ++f )
		  {
		     cout << "Here in tag match!" << endl;
		     const Candidate *tagMatchRef = &*(f->key);
		     if( &(*tagRef)==&(*tagMatchRef) ) 
		     {
			//const Candidate *mcMatchRef  = &*(f->val);
			CandidateRef mcMatchRef  = f->val;
		  
			if( mcMatchRef.isNonnull() )
			{
			   cout << "Id " << mcMatchRef->pdgId() << endl;

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
			   if( moid == 23 || gmoid == 23 ) tagFromZ = true;
			}
		     }
		  }
	       }
	    }

	    bool probeFromZ = false;
	    if( allprobemuons.isValid() )
	    {
	       CandidateRef probeRef( allprobemuons, probes[0].key() );
	       if( allprobematch.isValid() )
	       {
		  CandMatchMap::const_iterator f = allprobematch->begin();
		  for( ; f != allprobematch->end(); ++f )
		  {
		     cout << "Here in probe match!" << endl;
		     const Candidate *probeMatchRef = &*(f->key);
		     if( &(*probeRef)==&(*probeMatchRef) ) 
		     {
			CandidateRef mcMatchRef  = f->val;
		  
			if( mcMatchRef.isNonnull() )
			{
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
			   if( moid == 23 || gmoid == 23 ) probeFromZ = true;
			}
		     }
		  }
	       }
	    }
	    // If both tag and probe are from Z .. set to true
	    if( tagFromZ && probeFromZ ) tptrue = 1;
	    cout << "Is true Z? " << tptrue << endl;

	    // Is this probe in the set that pass the efficiency criteria?
	    int ppass = 0;
	    if( passprobemuons.isValid() )
	    {
	       for( int ipp=0; ipp<(int)passprobemuons->size(); ++ipp )
	       {
		  bool isOverlap = false;
		  MuonRef ppRef;
		  if( (*passprobemuons)[ipp].hasMasterClone() )
		  {
		     cout << "Pass probe has master clone!" << endl;
		     ppRef = ((*passprobemuons)[ipp].masterClone()).castTo<MuonRef>();
		  }
		  if( ppRef.isNonnull() )
		  {
		     if( (ppRef->track()).isNonnull() && 
			 (probes[0]->track()).isNonnull() && 
			 ppRef->track() == probes[0]->track() )
			isOverlap = true;
		     if( (ppRef->standAloneMuon()).isNonnull() && 
			 (probes[0]->standAloneMuon()).isNonnull() && 
			 ppRef->standAloneMuon() == probes[0]->standAloneMuon() )
			isOverlap = true;
		     if( (ppRef->combinedMuon()).isNonnull() && 
			 (probes[0]->combinedMuon()).isNonnull() && 
			 ppRef->combinedMuon() == probes[0]->combinedMuon() )
			isOverlap = true;
		     if( isOverlap ) ppass = 1;
		  }
	       } 
	    }

	    // Did this tag cause a L1 and/or HLT trigger?
	    bool l1Trigger = false;
	    if( l1MuonTrg.isValid() )
	    {
	       int npart = l1MuonTrg->size();
	       for(int ipart = 0; ipart != npart; ++ipart)
	       { 
		  l1Trigger = MatchObjects( l1MuonTrg->getParticleRef(ipart), tag );
		  if( l1Trigger ) break;
	       }
	    }

	    bool hltTrigger = false;
	    if( hltMuonTrg.isValid() )
	    {
	       int npart = hltMuonTrg->size();
	       for(int ipart = 0; ipart != npart; ++ipart)
	       { 
		  hltTrigger = MatchObjects( hltMuonTrg->getParticleRef(ipart), tag );
		  if( hltTrigger ) break;
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

	    m_evtTree.tp_type[nrtp]  = itype;
	    m_evtTree.tp_true[nrtp]  = tptrue;
	    m_evtTree.tp_ppass[nrtp] = ppass;
	    m_evtTree.tp_l1[nrtp] = l1Trigger;
	    m_evtTree.tp_hlt[nrtp] = hltTrigger;
	    m_evtTree.tp_p[nrtp]  = p;
	    m_evtTree.tp_px[nrtp] = px;
	    m_evtTree.tp_py[nrtp] = py;
	    m_evtTree.tp_pz[nrtp] = pz;
	    m_evtTree.tp_pt[nrtp] = pt;
	    m_evtTree.tp_e[nrtp]  = e;
	    m_evtTree.tp_et[nrtp] = et;
	    m_evtTree.tp_mass[nrtp] = mass;

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

	    m_evtTree.tp_dp[nrtp][0]   = dp;
	    m_evtTree.tp_dpx[nrtp][0]  = dpx;
	    m_evtTree.tp_dpy[nrtp][0]  = dpy;
	    m_evtTree.tp_dpz[nrtp][0]  = dpz;
	    m_evtTree.tp_dpt[nrtp][0]  = dpt;
	    m_evtTree.tp_de[nrtp][0]   = de;
	    m_evtTree.tp_det[nrtp][0]  = det;
	    m_evtTree.tp_dq[nrtp][0]   = dq;
	    m_evtTree.tp_deta[nrtp][0] = deta;
	    m_evtTree.tp_dphi[nrtp][0] = dphi;

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

	    m_evtTree.tp_dp[nrtp][1]   = dp;
	    m_evtTree.tp_dpx[nrtp][1]  = dpx;
	    m_evtTree.tp_dpy[nrtp][1]  = dpy;
	    m_evtTree.tp_dpz[nrtp][1]  = dpz;
	    m_evtTree.tp_dpt[nrtp][1]  = dpt;
	    m_evtTree.tp_de[nrtp][1]   = de;
	    m_evtTree.tp_det[nrtp][1]  = det;
	    m_evtTree.tp_dq[nrtp][1]   = dq;
	    m_evtTree.tp_deta[nrtp][1] = deta;
	    m_evtTree.tp_dphi[nrtp][1] = dphi;

	    nrtp++;
	 }
      }
   }

   m_evtTree.nrtp = nrtp;
}
// ******************************************************************* //


// ********************* Fill Muon Efficiency Info *********************** //
void
MuonTagProbeAnalyzer::fillMuonEffInfo()
{

   // Should change this to get the eff info for all types of tag-probe!!
   Handle<CandidateCollection> genparticles;
   try{ m_event->getByLabel(genParticlesTag_,genparticles); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract gen particles with input tag " 
			<< genParticlesTag_;
   }

   Handle<CandMatchMap> tagmatch;
   try{ m_event->getByLabel(tagTruthMatchMapTags_[0],tagmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract muon match map "
		      << "with input tag " << tagTruthMatchMapTags_[0];
   }

   Handle<CandMatchMap> apmatch;
   try{ m_event->getByLabel(allProbeTruthMatchMapTags_[0],apmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract all probe match map "
		      << "with input tag " << allProbeTruthMatchMapTags_[0];
   }

   Handle<CandMatchMap> ppmatch;
   try{ m_event->getByLabel(passProbeTruthMatchMapTags_[0],ppmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract pass probe match map "
		      << "with input tag " << passProbeTruthMatchMapTags_[0];
   }
   
   // Fill some information about the muon efficiency
   int ngmc = 0;
   if( genparticles.isValid() )
   {
      for( unsigned int i=0; i<genparticles->size(); i++ )
      {
	 if( ngmc >= EvtTree::MAXNMU ) break;

	 int pdg_id = (*genparticles)[i].pdgId();

	 // If this is not a muon keep going!
	 if( abs( pdg_id ) != 13 ) continue;

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
		  const Candidate *gmu = &*(f->key);
		  rp   = gmu->p();
		  rpx  = gmu->px();
		  rpy  = gmu->py();
		  rpz  = gmu->pz(); 
		  rpt  = gmu->pt(); 
		  re   = gmu->energy();
		  ret  = gmu->et();
		  rq   = gmu->charge();
		  rphi = gmu->phi();
		  reta = gmu->eta();
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
		  fppb = 1;
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
		  fapb = 1;
		  if( ftag == 0 && fppb == 0 )
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

	 m_evtTree.gmu_tag[ngmc]    = ftag;
	 m_evtTree.gmu_aprobe[ngmc] = fapb;
	 m_evtTree.gmu_pprobe[ngmc] = fppb;
	 m_evtTree.gmu_moid[ngmc] = moid;
	 m_evtTree.gmu_gmid[ngmc] = gmoid;

	 m_evtTree.gmu_p[ngmc] = p;
	 m_evtTree.gmu_px[ngmc] = px;
	 m_evtTree.gmu_py[ngmc] = py;
	 m_evtTree.gmu_pz[ngmc] = pz;
	 m_evtTree.gmu_pt[ngmc] = pt;
	 m_evtTree.gmu_e[ngmc] = e;
	 m_evtTree.gmu_et[ngmc] = et;
	 m_evtTree.gmu_q[ngmc] = q;
	 m_evtTree.gmu_phi[ngmc] = phi;
	 m_evtTree.gmu_eta[ngmc] = eta;

	 m_evtTree.gmu_rp[ngmc] = rp;
	 m_evtTree.gmu_rpx[ngmc] = rpx;
	 m_evtTree.gmu_rpy[ngmc] = rpy;
	 m_evtTree.gmu_rpz[ngmc] = rpz;
	 m_evtTree.gmu_rpt[ngmc] = rpt;
	 m_evtTree.gmu_re[ngmc] = re;
	 m_evtTree.gmu_ret[ngmc] = ret;
	 m_evtTree.gmu_rq[ngmc] = rq;
	 m_evtTree.gmu_rphi[ngmc] = rphi;
	 m_evtTree.gmu_reta[ngmc] = reta;

	 ngmc++;
      }
   }
   m_evtTree.ngmu = ngmc;

}
// *********************************************************************** //


// ********************* Fill Reco Z Info *********************** //
void
MuonTagProbeAnalyzer::fillZInfo()
{
   Handle<CandMatchMap> tagmatch;
   try{ m_event->getByLabel(tagTruthMatchMapTags_[0],tagmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract tag muon match map "
		      << "with input tag " << tagTruthMatchMapTags_[0] << endl;
   }

   Handle<CandMatchMap> apmatch;
   try{ m_event->getByLabel(allProbeTruthMatchMapTags_[0],apmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract all-probe match map "
		      << "with input tag " << allProbeTruthMatchMapTags_[0] << endl;
   }

   Handle<CandMatchMap> ppmatch;
   try{ m_event->getByLabel(passProbeTruthMatchMapTags_[0],ppmatch); }
   catch(...)
   { 
      LogWarning("Z") << "Could not extract pass-probe match map "
		      << "with input tag " << passProbeTruthMatchMapTags_[0] << endl;
   }

   vector<const CandMatchMap *> matchmaps;

   int nrz = 0;
   for( int itype=0; itype<(int)zCandTags_.size(); ++itype )
   {
      Handle<CandidateCollection> zcands;
      try{ m_event->getByLabel(zCandTags_[itype],zcands); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract reco Z -> tag tag cands "
			 << "with input tag " << zCandTags_[itype] << endl;
      }

      
      //if( zcands.isValid() ) cout << "Number Z's " << zcands->size() << endl;

      matchmaps.clear();
      if( itype == 0 )
      {
	 if( tagmatch.isValid() ) matchmaps.push_back(tagmatch.product());
	 if( tagmatch.isValid() ) matchmaps.push_back(tagmatch.product());
      }
      else if( itype == 1 )
      {
	 if( tagmatch.isValid() ) matchmaps.push_back(tagmatch.product());
	 if( apmatch.isValid() ) matchmaps.push_back(apmatch.product());
      }
      else
      {
	 if( tagmatch.isValid() ) matchmaps.push_back(tagmatch.product());
	 if( ppmatch.isValid() ) matchmaps.push_back(ppmatch.product());
      }

      // Fill some information about the reconstructed Zmumu vertices
      if( zcands.isValid() && matchmaps.size() == 2 )
      {
	 MCCandMatcher<CandidateCollection> match( matchmaps );
	 
	 for( unsigned int i=0; i<zcands->size(); i++ )
	 {
	    if( nrz >= EvtTree::MAXNPAIR ) break;
	    
	    double rp    = (*zcands)[i].p();
	    double rpt   = (*zcands)[i].pt();
	    double rpx   = (*zcands)[i].px();
	    double rpy   = (*zcands)[i].py();
	    double rpz   = (*zcands)[i].pz();
	    double re    = (*zcands)[i].energy();
	    double rmass = (*zcands)[i].mass();
	    double rx    = (*zcands)[i].vx();
	    double ry    = (*zcands)[i].vy();
	    double rz    = (*zcands)[i].vz();

	    double drx = 0;
	    double dry = 0;
	    double drz = 0;

	    CandidateRef mcCand = match( (*zcands)[i] );
	    bool isZ  = mcCand.isNonnull();
	    cout << "Is Z? " << isZ << endl;
	    if( isZ )
	    {
	       drx = mcCand->vx()-rx;
	       dry = mcCand->vy()-ry;
	       drz = mcCand->vz()-rz;
	    }

	    double chi2 = -1.0;
	    double ndof = 0.0;
	    
	    m_evtTree.z_type[nrz] = itype;
	    m_evtTree.z_true[nrz] = (int)isZ;
	    m_evtTree.z_p[nrz] = rp;
	    m_evtTree.z_pt[nrz] = rpt;
	    m_evtTree.z_px[nrz] = rpx;
	    m_evtTree.z_py[nrz] = rpy;
	    m_evtTree.z_pz[nrz] = rpz;
	    m_evtTree.z_e[nrz] = re;
	    m_evtTree.z_mass[nrz] = rmass;
	    cout << "Z mass " << rmass << " true Z " << isZ << endl;

	    m_evtTree.z_vvalid[nrz] = 1;
	    m_evtTree.z_vchi2[nrz] = chi2;
	    m_evtTree.z_vndof[nrz] = ndof;
	    m_evtTree.z_vx[nrz] = rx;
	    m_evtTree.z_vy[nrz] = ry;
	    m_evtTree.z_vz[nrz] = rz;
	    m_evtTree.z_vdx[nrz] = drx;
	    m_evtTree.z_vdy[nrz] = dry;
	    m_evtTree.z_vdz[nrz] = drz;
	    m_evtTree.z_vd0[nrz] = sqrt(rx*rx + ry*ry + rz*rz);

	    // Get the daughter info!
	    for( int id=0; id<2; ++id )
	    {
	       const Candidate *d = (*zcands)[i].daughter(id);

	       m_evtTree.z_dp[nrz][id]   = d->p();
	       m_evtTree.z_dpx[nrz][id]  = d->px();
	       m_evtTree.z_dpy[nrz][id]  = d->py();
	       m_evtTree.z_dpz[nrz][id]  = d->pz();
	       m_evtTree.z_dpt[nrz][id]  = d->pt();
	       m_evtTree.z_de[nrz][id]   = d->energy();
	       m_evtTree.z_det[nrz][id]  = d->et();
	       m_evtTree.z_dq[nrz][id]   = d->charge();
	       m_evtTree.z_dphi[nrz][id] = d->phi();
	       m_evtTree.z_deta[nrz][id] = d->eta();
	    }
	    
	    ++nrz;
	 }
      }
   }
      
   m_evtTree.nrz = nrz;

}
// ************************************************************** //



// ********************* Fill Vertex Info *********************** //
void
MuonTagProbeAnalyzer::fillVertexInfo()
{
   // get primary vertex collection
   Handle<VertexCollection> vertices;
   try{ m_event->getByLabel(PvxtTag_, vertices);}
   catch(...)
   {
      LogWarning("Z") << "Could not extract Primery Vertex with input tag " 
			<< PvxtTag_ << endl;

   }

   int nVx(0);
   math::XYZPoint vtx(0.,0.,0.);

   if( vertices.isValid() && vertices->size() > 0 )
   {
      const VertexCollection * VertexData;
      VertexData = vertices.product();
     
      const Vertex* pVertex = &(*VertexData)[0];
      vtx = pVertex->position();
    
      if (DEBUGLVL == 1)
      {
	 cout << "Vertex collection: " << endl;
	 for (int i=0; i< (int) VertexData->size(); i++)
	 {
	    const Vertex* pVertex = &(*VertexData)[i];
	    cout << " Vertex index = " << i << ", x = " << pVertex->x()
		 << ", y = " << pVertex->y()<< ", z = " << pVertex->z() << endl;
	 }
      }
     

      nVx = VertexData->size();
      if (nVx <= 0)     m_evtTree.vtx_NormalizedChi2 = -1.;


      const Vertex* primVx;
      int PVnum(0);

      primVx = &(*VertexData)[PVnum];
      m_evtTree.vtx_IndexPrimaryVx = PVnum;
      m_evtTree.vtx_Chi2 = primVx->chi2() ;
      m_evtTree.vtx_Ndof = primVx->ndof() ;
      m_evtTree.vtx_NormalizedChi2 = primVx->normalizedChi2() ;
      m_evtTree.vtx_PVx = primVx->x() ;
      m_evtTree.vtx_PVy = primVx->y() ;
      m_evtTree.vtx_PVz = primVx->z();
      m_evtTree.vtx_PVdx = primVx->xError();
      m_evtTree.vtx_PVdy = primVx->yError();
      m_evtTree.vtx_PVdz = primVx->zError();
  
      // Store also the number and scalar pt sum of all associated tracks
      m_evtTree.vtx_PvnTracks = primVx->tracksSize();
      float ptsum = 0.;
      Vertex::trackRef_iterator itk = primVx->tracks_begin();
      for( ; itk != primVx->tracks_end(); ++itk ) 
      {
	 float eta = (*itk)->eta();
	 // make this a config var  
	 if (fabs(eta) < 3.5){
	    float pt = (*itk)->pt();
	    //     cout << "  track PT = " << pt << ", eta = " << eta << endl;
	    ptsum += pt;
	 }
      }
      m_evtTree.vtx_PvPtsum = ptsum;

   }
   m_evtTree.nvtx = nVx;

   m_vtx = &vtx;

}
// ************************************************************** //

// ***************** Trigger object matching ******************** //
bool MuonTagProbeAnalyzer::MatchObjects( CandidateBaseRef hltObj, const MuonRef& tagObj )
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
//DEFINE_FWK_MODULE(MuonTagProbeAnalyzer);
