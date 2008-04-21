// -*- C++ -*-
//
// Package:    TagAndProbe
// Class:      TagProbeAnalyzer
// 
/**\class TagAndProbe TagProbeAnalyzer.cc MuonAnalysis/TagAndProbe/src/TagProbeAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nadia Adam
//         Created:  Tue Mar  4 12:27:53 CST 2008
// $Id: TagProbeAnalyzer.cc,v 1.2 2008/04/14 19:21:42 neadam Exp $
//
//


// system include files
#include <cmath>

// user include files
#include "MuonAnalysis/TagAndProbe/interface/TagProbeAnalyzer.h"

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
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
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
TagProbeAnalyzer::TagProbeAnalyzer(const ParameterSet& iConfig)
{

   rootFile_ = iConfig.getUntrackedParameter<string>("rootFile");

   candType_ = iConfig.getUntrackedParameter<string>("tagProbeType","Muon");
   candPDGId_ = 13;
   if( candType_ != "Muon" ) candPDGId_ = 11;

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


TagProbeAnalyzer::~TagProbeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TagProbeAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
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

   // Fill Tag-Probe Info
   fillTagProbeInfo();

   // Fill Efficiency Info for true muons
   fillTrueEffInfo();

   // Fill Vertex Info
   fillVertexInfo();

   // Finally fill the tree for this event!
   m_fillEvt.fill();   

   return;
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
TagProbeAnalyzer::beginJob(const EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TagProbeAnalyzer::endJob() {

   m_fillEvt.finalize();
}

// ****************************************************************** //
// ********************* Tree Fill Functions! *********************** //
// ****************************************************************** //


// ********************* Fill Run & Event Info *********************** //
void
TagProbeAnalyzer::fillRunEventInfo()
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
TagProbeAnalyzer::fillTriggerInfo()
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
      cout << "Could not extract HLT object!!" << endl;
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
TagProbeAnalyzer::fillMCInfo()
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
TagProbeAnalyzer::fillTrackInfo()
{

   // Fill some information about the reconstructed tracks
   int nrtrk = 0;
   for( int tcoll=0; tcoll<(int)trackTags_.size(); ++tcoll )
   {
      // Fill info for all tracks
      Handle<TrackCollection> tracks;
      try{ m_event->getByLabel(trackTags_[tcoll],tracks); }
      catch(...)
      { 
	 LogWarning("Z") << "Could not extract reco tracks with input tag " 
			 << trackTags_[tcoll];
      }

      if( tracks.isValid() )
      {
	 for( unsigned int i=0; i<tracks->size(); i++ )
	 {
	    if( nrtrk >= EvtTree::MAXNTRK ) break;

	    int id = i+1;
	    int type = tcoll;
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
   }

   m_evtTree.nrtrk = nrtrk;
}
// ************************************************************* //




// ********************* Fill Tag-Probe Info *********************** //
void
TagProbeAnalyzer::fillTagProbeInfo()
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
	    if( nrtp >= EvtTree::MAXNPAIR ) break;

	    const CandidateRef &tag = tpItr->key;
	    const CandidateRefVector &probes = tpItr->val;
	    //cout << "Tag key " << tag.key() << endl;

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


// ********************* Fill True Efficiency Info *********************** //
void
TagProbeAnalyzer::fillTrueEffInfo()
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

	 m_evtTree.cnd_tag[ngmc]    = ftag;
	 m_evtTree.cnd_aprobe[ngmc] = fapb;
	 m_evtTree.cnd_pprobe[ngmc] = fppb;
	 m_evtTree.cnd_moid[ngmc] = moid;
	 m_evtTree.cnd_gmid[ngmc] = gmoid;

	 m_evtTree.cnd_p[ngmc] = p;
	 m_evtTree.cnd_px[ngmc] = px;
	 m_evtTree.cnd_py[ngmc] = py;
	 m_evtTree.cnd_pz[ngmc] = pz;
	 m_evtTree.cnd_pt[ngmc] = pt;
	 m_evtTree.cnd_e[ngmc] = e;
	 m_evtTree.cnd_et[ngmc] = et;
	 m_evtTree.cnd_q[ngmc] = q;
	 m_evtTree.cnd_phi[ngmc] = phi;
	 m_evtTree.cnd_eta[ngmc] = eta;

	 m_evtTree.cnd_rp[ngmc] = rp;
	 m_evtTree.cnd_rpx[ngmc] = rpx;
	 m_evtTree.cnd_rpy[ngmc] = rpy;
	 m_evtTree.cnd_rpz[ngmc] = rpz;
	 m_evtTree.cnd_rpt[ngmc] = rpt;
	 m_evtTree.cnd_re[ngmc] = re;
	 m_evtTree.cnd_ret[ngmc] = ret;
	 m_evtTree.cnd_rq[ngmc] = rq;
	 m_evtTree.cnd_rphi[ngmc] = rphi;
	 m_evtTree.cnd_reta[ngmc] = reta;

	 ngmc++;
      }
   }
   m_evtTree.ncnd = ngmc;

}
// *********************************************************************** //


// ********************* Fill Vertex Info *********************** //
void
TagProbeAnalyzer::fillVertexInfo()
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

// ***************** Are Candidates from a Z? ******************** //
bool TagProbeAnalyzer::CandFromZ( const reco::CandidateRef& cand, 
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
int TagProbeAnalyzer::ProbePassProbeOverlap( const CandidateRef& probe, 
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


// ***************** Trigger object matching ******************** //
bool TagProbeAnalyzer::MatchObjects( CandidateBaseRef hltObj, const CandidateRef& tagObj )
{
   double tEta = tagObj->eta();
   double tPhi = tagObj->phi();
   double hEta = hltObj->eta();
   double hPhi = hltObj->phi();

   double dRval = deltaR(tEta, tPhi, hEta, hPhi);
   return ( dRval < delRMatchingCut_ );
}
// ************************************************************** //

// ***************** Trigger object matching ******************** //
bool TagProbeAnalyzer::MatchObjects( const Candidate *hltObj, 
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
//DEFINE_FWK_MODULE(TagProbeAnalyzer);
