//
// $Id: HighPtMuonsInfo.cc,v 1.6 2011/07/21 01:07:49 botta Exp $
//

/**
  \class    HighPtMuonsInfo HighPtMuonsInfo.h "PhysicsTools/PatAlgos/interface/HighPtMuonsInfo.h"
  \brief    Interface new TuneP to the T&P.
            https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#New_Version_recommended
            
  \author   Giovanni Petrucciani
  \version  $Id: HighPtMuonsInfo.cc,v 1.6 2011/07/21 01:07:49 botta Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include <DataFormats/Candidate/interface/CompositeCandidate.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

class HighPtMuonsInfo : public edm::EDProducer {
    public:
      explicit HighPtMuonsInfo(const edm::ParameterSet & iConfig);
      virtual ~HighPtMuonsInfo() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
    private:
        edm::InputTag src_;

        /// Write a ValueMap<float> in the event
        void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Candidate> > & handle,
                const std::vector<float> & values,
                const std::string    & label) const ;

};

HighPtMuonsInfo::HighPtMuonsInfo(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<edm::ValueMap<float> >("pt");
    produces<edm::ValueMap<float> >("ptRelError");
    produces<edm::ValueMap<float> >("mass");
    produces<edm::ValueMap<float> >("trackType");
}

void 
HighPtMuonsInfo::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<reco::Candidate> > src;
    iEvent.getByLabel(src_, src);

    size_t n = src->size();
    std::vector<float> pt(n,0), ptRelError(n,0), mass(n, -999), trackType(n,-999);
    for (size_t i = 0; i < n; ++i) {
        const reco::Candidate & ci = (*src)[i];
        if (ci.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << 
            "HighPtMuonsInfo should be used on composite candidates with two daughters, this one has " << ci.numberOfDaughters() << "\n";
        const reco::Candidate &d1 = *ci.daughter(0), &d2 = *ci.daughter(1);

        const reco::Muon *mu1 = dynamic_cast<const reco::Muon *>(&*d1.masterClone());
        if (mu1 == 0) throw cms::Exception("CorruptData") << "First daughter of candidate is not a ShallowClone of a reco::Muon\n";
        reco::Muon::MuonTrackTypePair tuneP1 = muon::tevOptimized(*mu1, 200, 40., 17., 0.25);

        const reco::Muon *mu2 = dynamic_cast<const reco::Muon *>(&*d2.masterClone());
        if (mu2 == 0) throw cms::Exception("CorruptData") << "Second daughter of candidate is not a ShallowClone of a reco::Muon\n";
        reco::Muon::MuonTrackTypePair tuneP2 = muon::tevOptimized(*mu2, 200, 40., 17., 0.25);

        // Momentum and relative uncertainty
        pt[i] = tuneP2.first->pt();
        ptRelError[i] = tuneP2.first->ptError()/pt[i];

        // Invariant mass
        mass[i] = ( reco::Particle::PolarLorentzVector(tuneP2.first->pt(), tuneP2.first->eta(), tuneP2.first->phi(), d2.mass()) +
                    reco::Particle::PolarLorentzVector(tuneP2.first->pt(), tuneP2.first->eta(), tuneP2.first->phi(), d2.mass())  ).M();
 
        // convert the enumerator into int by hand: numbers are not given explicitly in the header file, and might change */
        switch (tuneP2.second) {
            case reco::Muon::None: trackType[i] = 0; break;
            case reco::Muon::InnerTrack: trackType[i] = 1; break;
            case reco::Muon::OuterTrack: trackType[i] = 2; break;
            case reco::Muon::CombinedTrack: trackType[i] = 3; break;
            case reco::Muon::TPFMS: trackType[i] = 4; break;
            case reco::Muon::Picky: trackType[i] = 5; break;
            case reco::Muon::DYT: trackType[i] = 6; break;
        }
    }

    writeValueMap(iEvent, src, pt, "pt");
    writeValueMap(iEvent, src, ptRelError, "ptRelError");
    writeValueMap(iEvent, src, mass, "mass");
    writeValueMap(iEvent, src, trackType, "trackType");
}

void
HighPtMuonsInfo::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Candidate> > & handle,
        const std::vector<float> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;
    auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HighPtMuonsInfo);
