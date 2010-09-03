//
// $Id: RedefineMuonP4FromTrackT.cc,v 1.1 2009/05/28 17:33:26 gpetrucc Exp $
//

/**
  \class    RedefineMuonP4FromTrackT RedefineMuonP4FromTrackT.h "PhysicsTools/PatAlgos/interface/RedefineMuonP4FromTrackT.h"
  \brief    Use another track to define the 4-momentum of a Muon object
            
  \author   Giovanni Petrucciani
  \version  $Id: RedefineMuonP4FromTrackT.cc,v 1.1 2009/05/28 17:33:26 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

template<typename T>
class RedefineMuonP4FromTrackT : public edm::EDProducer {
    public:
      explicit RedefineMuonP4FromTrackT(const edm::ParameterSet & iConfig);
      virtual ~RedefineMuonP4FromTrackT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      enum TkType { TrackerTk, MuonTk, GlobalTk } ;
      edm::InputTag src_;
      TkType        track_;
      
};

template<typename T>
RedefineMuonP4FromTrackT<T>::RedefineMuonP4FromTrackT(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    std::string track = iConfig.getParameter<std::string>("track");
    if (track == "tracker" || track == "inner" || track == "silicon") track_ = TrackerTk;
    else if (track == "muon" || track == "outer" || track == "standAlone") track_ = MuonTk;
    else if (track == "global" || track == "combined") track_ = GlobalTk;
    else throw cms::Exception("Configuration") << "Track type '" << track << "' not suported. Allowed values are 'tracker', 'muon', 'global'.\n";
    produces<std::vector<T> >();
}

template<typename T>
void 
RedefineMuonP4FromTrackT<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<T> >      src;
    iEvent.getByLabel(src_, src);

    auto_ptr<vector<T> >  out(new vector<T>());
    out->reserve(src->size());

    for (typename View<T>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
        const T & original = *it;
        reco::TrackRef tkref;
        switch (track_) {
            case TrackerTk: tkref = original.track();                   break;
            case MuonTk:    tkref = original.standAloneMuon(); break;
            case GlobalTk:  tkref = original.combinedMuon();   break;
        }
        if (tkref.isNull()) continue;
        const reco::Track &trk = *tkref;
        out->push_back(original);
        T & mu = out->back();
        mu.setCharge(trk.charge());
        mu.setVertex(trk.vertex());
        mu.setP4(reco::Candidate::PolarLorentzVector(trk.pt(), trk.eta(), trk.phi(), original.mass()));
    }

    iEvent.put(out);
}

typedef RedefineMuonP4FromTrackT<reco::Muon> RedefineMuonP4FromTrack;
typedef RedefineMuonP4FromTrackT< pat::Muon> RedefineMuonP4FromTrackPAT;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RedefineMuonP4FromTrack);
DEFINE_FWK_MODULE(RedefineMuonP4FromTrackPAT);
