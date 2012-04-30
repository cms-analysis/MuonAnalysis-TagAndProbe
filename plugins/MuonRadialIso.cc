// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//
// class declaration
//

class MuonRadialIso : public edm::EDProducer {
public:
  explicit MuonRadialIso(const edm::ParameterSet&);
  ~MuonRadialIso();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
  const edm::InputTag pfCandidates_;
  const double photonPtMin_, neutralHadPtMin_;


  /// Store extra information in a ValueMap
  template<typename Hand, typename T>
  void storeMap(edm::Event &iEvent, 
                 const Hand & handle,
                 const std::vector<T> & values,
                 const std::string    & label) const ;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonRadialIso::MuonRadialIso(const edm::ParameterSet& iConfig):
    probes_(iConfig.getParameter<edm::InputTag>("probes")),
    pfCandidates_(iConfig.getParameter<edm::InputTag>("pfCandidates")),
    photonPtMin_(iConfig.getParameter<double>("photonPtMin")),
    neutralHadPtMin_(iConfig.getParameter<double>("neutralHadPtMin"))
{
  produces<edm::ValueMap<float> >();
}


MuonRadialIso::~MuonRadialIso()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonRadialIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByLabel(probes_, probes);

  Handle<View<reco::PFCandidate> > pfCandidates;
  iEvent.getByLabel(pfCandidates_, pfCandidates);

  View<reco::Muon>::const_iterator probe, endprobes=probes->end();
  View<reco::PFCandidate>::const_iterator iP, beginpf = pfCandidates->begin(), endpf=pfCandidates->end();
  unsigned int n = probes->size();

  std::vector<float> iso(n,0);

  // loop on PROBES
  unsigned int imu = 0;
  for (probe = probes->begin(); probe != endprobes; ++probe, ++imu) {
    const reco::Muon &mu = *probe;

    for (iP = beginpf; iP != endpf; ++iP) {
        double dr = deltaR(*iP, mu);

        bool hastk = iP->trackRef().isNonnull();
        if (hastk && mu.track() == iP->trackRef()) continue;

        if (dr < 0.01 || dr > 0.3) continue;
        
        if (hastk) {
            if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
        } else if (iP->particleId() == reco::PFCandidate::gamma) {
            if (iP->pt() <= photonPtMin_) continue;
        } else {
            if (iP->pt() <= neutralHadPtMin_) continue;
        }
        iso[imu] += iP->pt() * ( 1 - 3*dr) / mu.pt();
    }
  }// end loop on probes

  storeMap(iEvent, probes, iso, "");
}

template<typename Hand, typename T>
void
MuonRadialIso::storeMap(edm::Event &iEvent,
                     const Hand & handle,
                     const std::vector<T> & values,
                     const std::string    & label) const {
    using namespace edm; using namespace std;
    auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonRadialIso);
