// -*- C++ -*-
//
// Package:    MuonTagProbeProducer
// Class:      MuonTagProbeProducer
// 
/**\class MuonTagProbeProducer MuonTagProbeProducer.cc MuonEffAnalysis/MuonTagProbeProducer/src/MuonTagProbeProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nadia Adam
//         Created:  Wed Dec 19 17:32:47 CST 2007
// $Id$
//
//


// system include files

// user include files
#include "MuonAnalysis/TagAndProbe/interface/MuonTagProbeProducer.h"
#include "AnalysisDataFormats/Muon/interface/MuonTagProbeAssociation.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Math/GenVector/VectorUtil.h"


//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonTagProbeProducer::MuonTagProbeProducer(const edm::ParameterSet& iConfig)
{
   tagCollection_   = iConfig.getParameter<edm::InputTag>("TagCollection");
   probeCollection_ = iConfig.getParameter<edm::InputTag>("ProbeCollection");

   massMinCut_      = iConfig.getUntrackedParameter<double>("MassMinCut",50.0);
   massMaxCut_      = iConfig.getUntrackedParameter<double>("MassMaxCut",120.0);
   delRMinCut_      = iConfig.getUntrackedParameter<double>("DelRMinCut",0.0);
   delRMaxCut_      = iConfig.getUntrackedParameter<double>("DelRMaxCut",10000.0);

   requireOS_       = iConfig.getUntrackedParameter<bool>("RequireOS",true);

   produces<reco::MuonTagProbeAssociationCollection>();
}


MuonTagProbeProducer::~MuonTagProbeProducer()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonTagProbeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   // We need the output Muon association collection to fill
   std::auto_ptr<MuonTagProbeAssociationCollection> 
      muonTPCollection( new MuonTagProbeAssociationCollection );
 
   // Read in the tag muons
   edm::Handle<MuonCollection> tags;
   try{ iEvent.getByLabel( tagCollection_, tags ); }
   catch(...)
   {
      LogWarning("MuonTagProbe") << "Could not extract tag muons with input tag "
				 << tagCollection_;
   }

   // Read in the probe muons
   edm::Handle<MuonCollection> probes;
   try{ iEvent.getByLabel( probeCollection_, probes ); }
   catch(...)
   {
      LogWarning("MuonTagProbe") << "Could not extract probe muons with input tag "
				 << probeCollection_;
   }


   // Loop over Tag and associate with Probes
   if( tags.isValid() && probes.isValid() )
   {
      int itag = 0;
      MuonCollection::const_iterator tag = (*tags).begin();
      for( ; tag != (*tags).end(); ++tag ) 
      {  
	 MuonRef tagRef(tags,itag);
	 ++itag;

	 int iprobe = 0;
	 MuonCollection::const_iterator probe = (*probes).begin();
	 for( ; probe != (*probes).end(); ++probe ) 
	 {
	    MuonRef probeRef(probes,iprobe);
	    ++iprobe;
	    
	    // Tag-Probe invariant mass cut
	    double invMass = ROOT::Math::VectorUtil::InvariantMass(tag->p4(), probe->p4());
	    if( invMass < massMinCut_ ) continue;
	    if( invMass > massMaxCut_ ) continue;

	    // Tag-Probe deltaR cut
	    double delR = deltaR<double>(tag->eta(),tag->phi(),probe->eta(),probe->phi());
	    if( delR < delRMinCut_ ) continue;
	    if( delR > delRMaxCut_ ) continue;

	    // Tag-Probe opposite sign
	    int sign = tag->charge() * probe->charge();
	    if( requireOS_ && sign > 0 ) continue;

            muonTPCollection->insert( tagRef, probeRef );

	 }	 
      }
   }

   // Finally put the tag probe collection in the event
   iEvent.put( muonTPCollection );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonTagProbeProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonTagProbeProducer::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(MuonTagProbeProducer);
