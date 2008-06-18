#include "MuonAnalysis/TagAndProbe/interface/eTriggerCandProducer.h"
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <string>

using namespace std; 
using namespace reco; 
using namespace edm;
using namespace l1extra;
using namespace trigger;


eTriggerCandProducer::eTriggerCandProducer(const edm::ParameterSet& iConfig )
{

  _inputProducer = iConfig.getParameter<std::string>("InputProducer");

   // **************** Trigger ******************* //
   const edm::InputTag dTriggerEventTag("triggerSummaryRAW");
   triggerEventTag_ = 
      iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);

   const edm::InputTag dHLTL1Tag("l1IsoElectronPixelSeeds");
   hltL1Tag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltL1Tag",dHLTL1Tag);

   const edm::InputTag dHLTTag("HLTElectronPixelMatchFilter");
   hltTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hltTag",dHLTTag);

   delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch",0.15);
   // ******************************************** //



   produces<reco::GsfElectronCollection>();
}




eTriggerCandProducer::~eTriggerCandProducer()
{

}


//
// member functions
//


// ------------ method called to produce the data  ------------

void eTriggerCandProducer::produce(edm::Event &event, const edm::EventSetup &eventSetup)
{
   // Create the output collection
   std::auto_ptr<reco::GsfElectronCollection> outCol(new reco::GsfElectronCollection);



   // Get the input GsfElectron collection
   edm::Handle<reco::GsfElectronCollection> eleCandidatesHandle;
   try
   {
      event.getByLabel(_inputProducer, eleCandidatesHandle);
   }
   catch(cms::Exception &ex)
   {
      edm::LogError("GsfElectron ") << "Error! Can't get collection " << _inputProducer;
      throw ex;
   }



   // Loop over electrons
   for(unsigned int i = 0; i < eleCandidatesHandle->size(); ++i) {
     // Get cut decision for each electron
     edm::Ref<reco::GsfElectronCollection> electronRef(eleCandidatesHandle, i);

     reco::GsfElectron electron = *electronRef;
     bool boolDecision = TriggerDecision(event, electron);

     if(boolDecision) outCol->push_back(*electronRef);
   } 

   event.put(outCol);
}







// ***************** Trigger object matching ******************** //
bool eTriggerCandProducer::TriggerDecision(edm::Event &event, 
					   reco::GsfElectron& electron) {


  // Trigger Info
  Handle<TriggerEventWithRefs> trgEvent;
  try{ event.getByLabel(triggerEventTag_,trgEvent); } 
  catch (...) { 
    LogWarning("TagAndProbe") << "Could not extract trigger event summary "
			      << "with tag " << triggerEventTag_;
  }
  

  // Did this tag cause a L1 and/or HLT trigger?
  bool boolDecision = false;
  
  if( trgEvent.isValid() ) {

    bool l1Trigger = false;
    bool hltTrigger = false;
       
    // L1 Trigger
    vector< L1EmParticleRef   > emL1CandRefs;
    size_type l1index = trgEvent->filterIndex(hltL1Tag_.label());
    if( l1index < trgEvent->size() ) {
      trgEvent->getObjects( l1index, trigger::TriggerL1IsoEG, emL1CandRefs );
      int npart = emL1CandRefs.size();
      for(int ipart = 0; ipart != npart; ++ipart) { 
	const L1EmParticle* l1e = (emL1CandRefs[ipart]).get();
	l1Trigger = MatchObjects( (const Candidate*)l1e, electron);
	if( l1Trigger ) break;
      }
    }
       

    // HLT Trigger
    vector< RecoChargedCandidateRef > theCandRefs;
    size_type index = trgEvent->filterIndex(hltTag_.label());
    if( index < trgEvent->size() ) {
      trgEvent->getObjects( index, trigger::TriggerElectron, theCandRefs );
      int npart = theCandRefs.size();
      for(int ipart = 0; ipart != npart; ++ipart) { 
	const RecoChargedCandidate* thecand = (theCandRefs[ipart]).get();
	hltTrigger = MatchObjects( (const Candidate*)thecand, electron);
	if( hltTrigger ) break;
      }
    }

    boolDecision = l1Trigger || hltTrigger;
  }
   
  return boolDecision;
}





// ***************** Trigger object matching ******************** //
bool eTriggerCandProducer::MatchObjects( const Candidate *hltObj, 
					 const reco::GsfElectron& tagObj)
{
   double tEta = tagObj.eta();
   double tPhi = tagObj.phi();
   double hEta = hltObj->eta();
   double hPhi = hltObj->phi();

   double dRval = deltaR(tEta, tPhi, hEta, hPhi);

   return ( dRval < delRMatchingCut_);
}
// ************************************************************** //









// ------------ method called once each job just before starting event loop  ---



void eTriggerCandProducer::beginJob(const edm::EventSetup &eventSetup) {
}




void eTriggerCandProducer::endJob() {
}



//define this as a plug-in
DEFINE_FWK_MODULE( eTriggerCandProducer );
