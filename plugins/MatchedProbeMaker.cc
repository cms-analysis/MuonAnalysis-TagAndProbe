#include "MuonAnalysis/TagAndProbe/src/MatchedProbeMaker.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

//
// constants, enums and typedefs
//
   using namespace edm;
   using namespace reco;
   using namespace std;

//
// constructors and destructor
//
MatchedProbeMaker::MatchedProbeMaker(const edm::ParameterSet& iConfig):
  m_candidateSource(iConfig.getUntrackedParameter<edm::InputTag>("CandidateSource")),
  m_referenceSource(iConfig.getUntrackedParameter<edm::InputTag>("ReferenceSource")),
  m_resMatchMapSource(iConfig.getUntrackedParameter<edm::InputTag>("ResMatchMapSource"))
{
  //register your products
  produces<CandidateCollection>(std::string("matched"));
  produces<CandidateCollection>(std::string("unmatched"));

  //now do what ever other initialization is needed

}

MatchedProbeMaker::~MatchedProbeMaker()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to produce the data  ------------
void
MatchedProbeMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<CandidateCollection> outputCollection_matched(new CandidateCollection());
  std::auto_ptr<CandidateCollection> outputCollection_unmatched(new CandidateCollection());
  
  // Get the candidates from the event
  Handle<CandidateCollection> Cands;
  iEvent.getByLabel(m_candidateSource,Cands);

  // Get the resolution matching map from the event
  Handle<CandMatchMap> ResMatchMap;
  iEvent.getByLabel(m_resMatchMapSource,ResMatchMap);
  
  LogDebug("ProbeProducer")<< "Cands size " << Cands->size();
  LogDebug("ProbeProducer")<< "Map size " << ResMatchMap->size();

  // Loop over the candidates looking for a match
  for (unsigned i=0; i<Cands->size(); i++){
    CandidateRef CandRef(Cands,i);
    LogDebug("ProbeProducer") << "MuonType sta: " <<  dynamic_cast<const Muon *>(&*CandRef)->isStandAloneMuon();

    // Loop over match map
    CandMatchMap::const_iterator f = ResMatchMap->find(CandRef);
    if (f!=ResMatchMap->end()){
      //const CandidateRef &CandMatch = f->val;
      LogDebug("ProbeProducer") << "Match found";
      outputCollection_matched->push_back(CandRef->clone());      
    }
    else {
      LogDebug("ProbeProducer")<< "Match not found";
      outputCollection_unmatched->push_back(CandRef->clone());      
    }
  }  

   iEvent.put(outputCollection_matched,"matched");
   iEvent.put(outputCollection_unmatched,"unmatched");

}

// ------------ method called once each job just before starting event loop  ------------
void
MatchedProbeMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MatchedProbeMaker::endJob() 
{
}

