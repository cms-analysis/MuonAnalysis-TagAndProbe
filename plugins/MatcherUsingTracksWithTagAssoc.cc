//
// $Id: MatcherUsingTracksWithTagAssoc.cc,v 1.5 2010/07/12 20:56:11 gpetrucc Exp $
//

/**
  \class    pat::MatcherUsingTracksWithTagAssoc MatcherUsingTracksWithTagAssoc.h "MuonAnalysis/MuonAssociators/interface/MatcherUsingTracksWithTagAssoc.h"
  \brief    Matcher of reconstructed objects to other reconstructed objects using the tracks inside them 
            
  \author   Giovanni Petrucciani
  \version  $Id: MatcherUsingTracksWithTagAssoc.cc,v 1.5 2010/07/12 20:56:11 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/UserData.h"

#include "MuonAnalysis/MuonAssociators/interface/MatcherUsingTracksAlgorithm.h"

// template-related workaround for bug in OwnVector+Ptr
namespace edm { using std::advance; }

namespace pat {

  class MatcherUsingTracksWithTagAssoc : public edm::EDProducer {
    public:
      explicit MatcherUsingTracksWithTagAssoc(const edm::ParameterSet & iConfig);
      virtual ~MatcherUsingTracksWithTagAssoc() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      /// Labels for input collections
      edm::InputTag src_, matched_;

      /// Tags to preselect probes on 
      edm::InputTag tags_;
      double deltaZ_;

      /// The real workhorse
      MatcherUsingTracksAlgorithm algo_;

      /// Some extra configurables
      bool dontFailOnMissingInput_;
      bool writeExtraPATOutput_;

      /// Store extra information in a ValueMap
      template<typename T>
      void storeValueMap(edm::Event &iEvent, 
                     const edm::Handle<edm::View<reco::Candidate> > & handle,
                     const std::vector<T> & values,
                     const std::string    & label) const ;

  };
  
} // namespace

pat::MatcherUsingTracksWithTagAssoc::MatcherUsingTracksWithTagAssoc(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    matched_(iConfig.getParameter<edm::InputTag>("matched")),
    tags_(iConfig.getParameter<edm::InputTag>("tags")),
    deltaZ_(iConfig.getParameter<double>("tagDeltaZ")),
    algo_(iConfig),
    dontFailOnMissingInput_(iConfig.existsAs<bool>("dontFailOnMissingInput") ? iConfig.getParameter<bool>("dontFailOnMissingInput") : false)
{
    // this is the basic output (edm::Association is not generic)
    produces<edm::ValueMap<reco::CandidatePtr> >(); 
    if (writeExtraPATOutput_) {
        // this is the crazy stuff to get the same with UserData
        produces<edm::OwnVector<pat::UserData> >();
        produces<edm::ValueMap<edm::Ptr<pat::UserData> > >();
    }
    // this is the extra stuff
    if (algo_.hasMetrics()) {
        produces<edm::ValueMap<float> >("deltaR");
        produces<edm::ValueMap<float> >("deltaEta");
        produces<edm::ValueMap<float> >("deltaPhi");
        produces<edm::ValueMap<float> >("deltaLocalPos");
        produces<edm::ValueMap<float> >("deltaPtRel");
        if (algo_.hasChi2()) {
            produces<edm::ValueMap<float> >("chi2");
        }
    } else {
        produces<edm::ValueMap<int> >("matched");
    }
}

void 
pat::MatcherUsingTracksWithTagAssoc::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    algo_.init(iSetup);

    Handle<View<reco::Candidate> > src, matched, tags;

    iEvent.getByLabel(src_, src);
    iEvent.getByLabel(matched_, matched);
    iEvent.getByLabel(tags_, tags);

    // declare loop variables and some intermediate stuff
    View<reco::Candidate>::const_iterator itsrc, edsrc = src->end(); 
    View<reco::Candidate>::const_iterator itmatched, edmatched = matched->end();
    int isrc, imatched, nsrc = src->size();

    std::vector<double> tagZs;
    for (View<reco::Candidate>::const_iterator ittags = tags->begin(), edtags = tags->end(); ittags != edtags; ++ittags) {
        tagZs.push_back(ittags->vz());
    }

    // working and output variables
    vector<int>   match(nsrc, -1);
    vector<float> deltaRs(nsrc, 999);
    vector<float> deltaEtas(nsrc, 999);
    vector<float> deltaPhis(nsrc, 999);
    vector<float> deltaPtRel(nsrc, 999);
    vector<float> deltaLocalPos(nsrc, 999);
    vector<float> chi2(nsrc, 999999);

    // don't try matching if the input collection is missing and the module is configured to fail silently
    if (!(matched.failedToGet() && dontFailOnMissingInput_)) {
        // loop on the source collection, and request for the match
        for (itsrc = src->begin(), isrc = 0; itsrc != edsrc; ++itsrc, ++isrc) {
            match[isrc] = -1;
            for (itmatched = matched->begin(), imatched = 0; itmatched != edmatched; ++itmatched, ++imatched) {
                double matchedz = itmatched->vz(); 
                bool isAssoc = false;
                for (std::vector<double>::const_iterator itz = tagZs.begin(), edz = tagZs.end(); itz != edz; ++itz) {
                    if (fabs(*itz - matchedz) < deltaZ_) { isAssoc = true; break; }
                }
                if (!isAssoc) continue;
                if (algo_.match(*itsrc, *matched, deltaRs[isrc], deltaEtas[isrc], deltaPhis[isrc], deltaLocalPos[isrc], deltaPtRel[isrc], chi2[isrc])) {
                    match[isrc] =  imatched;
                }
            }
        }
    }

    std::vector<reco::CandidatePtr> ptrs(nsrc);
    for (isrc = 0; isrc < nsrc; ++isrc) {
        if (match[isrc] != -1) {
            ptrs[isrc] = matched->ptrAt(match[isrc]);
        }
    }
    storeValueMap<reco::CandidatePtr>(iEvent, src, ptrs, "");

    if (algo_.hasMetrics()) {
        storeValueMap<float>(iEvent, src, deltaRs,   "deltaR");
        storeValueMap<float>(iEvent, src, deltaEtas, "deltaEta");
        storeValueMap<float>(iEvent, src, deltaPhis, "deltaPhi");
        storeValueMap<float>(iEvent, src, deltaLocalPos, "deltaLocalPos");
        storeValueMap<float>(iEvent, src, deltaPtRel,    "deltaPtRel");
        if (algo_.hasChi2()) {
            storeValueMap<float>(iEvent, src, chi2, "chi2");
        }
    } else  {
        std::vector<int> ismatched(nsrc, 0);
        for (isrc = 0; isrc < nsrc; ++isrc) {
            ismatched[isrc] = (match[isrc] != -1);
        }
        storeValueMap<int>(iEvent, src, ismatched, "matched");
    }
}

template<typename T>
void
pat::MatcherUsingTracksWithTagAssoc::storeValueMap(edm::Event &iEvent,
                     const edm::Handle<edm::View<reco::Candidate> > & handle,
                     const std::vector<T> & values,
                     const std::string    & label) const {
    using namespace edm; using namespace std;
    auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}


#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(MatcherUsingTracksWithTagAssoc);
