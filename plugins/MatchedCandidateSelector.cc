//
// $Id: MatchedCandidateSelector.cc,v 1.1 2009/09/30 15:26:02 gpetrucc Exp $
//

/**
  \class    MatchedCandidateSelector MatchedCandidateSelector.h "MuonAnalysis/TagAndProbe/plugins/MatchedCandidateSelector.h"
  \brief    Separate a list of candidates in matched ones and unmatched ones
            Inputs are:
                - a list of Candidates (param. "src")
                - a ValueMap<CandidatePtr>  with the matches (param. "match")
            Outputs are:
                - a list of references to candidate containing the items in 'src' that have a non-null match
                  (type CandidateBaseRefVector)
                - a list of references to candidate containing the items in 'src' that have non-null match
                  (type CandidateBaseRefVector, instance label "unmatched")
            Note: 
                - the "src" can be a strict subset of the keys of the "match" map.
            
  \author   Giovanni Petrucciani
  \version  $Id: MatchedCandidateSelector.cc,v 1.1 2009/09/30 15:26:02 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <DataFormats/Candidate/interface/CandidateFwd.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/Common/interface/ValueMap.h>

class MatchedCandidateSelector : public edm::EDProducer {
    public:
        explicit MatchedCandidateSelector(const edm::ParameterSet & iConfig);
        virtual ~MatchedCandidateSelector() { }

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::InputTag src_, match_;
};

MatchedCandidateSelector::MatchedCandidateSelector(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")), 
    match_(iConfig.getParameter<edm::InputTag>("match")) 
{
    produces<reco::CandidateBaseRefVector>();
    produces<reco::CandidateBaseRefVector>("unmatched");
}

void 
MatchedCandidateSelector::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    Handle<View<reco::Candidate> > src;
    iEvent.getByLabel(src_, src);
    Handle<ValueMap<reco::CandidatePtr> > match;
    iEvent.getByLabel(match_, match);

    std::auto_ptr<reco::CandidateBaseRefVector >  out( new reco::CandidateBaseRefVector());
    std::auto_ptr<reco::CandidateBaseRefVector >  fail(new reco::CandidateBaseRefVector());
    for (size_t i = 0, n = src->size(); i < n; ++i) {
        reco::CandidateBaseRef cbr = src->refAt(i);
        if ((*match)[cbr].isNonnull()) out->push_back(cbr); else fail->push_back(cbr);
    }

    iEvent.put(out);
    iEvent.put(fail, "unmatched");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MatchedCandidateSelector);
