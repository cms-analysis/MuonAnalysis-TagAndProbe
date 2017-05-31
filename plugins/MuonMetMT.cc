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
#include "DataFormats/METReco/interface/MET.h"

//
// class declaration
//

class MuonMetMT : public edm::EDProducer {
public:

  typedef edm::View<reco::MET> METColl;

  explicit MuonMetMT(const edm::ParameterSet&);
  ~MuonMetMT();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  const edm::EDGetTokenT<edm::View<reco::Muon>> probes_;    
  const edm::EDGetTokenT<METColl> pfMet_;

  template<typename Hand, typename T>
  void storeMap(edm::Event &iEvent, const Hand & handle, const std::vector<T> & values, const std::string    & label) const ; 
};

MuonMetMT::MuonMetMT(const edm::ParameterSet& iConfig):
probes_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("probes"))),
pfMet_(consumes<METColl>(iConfig.getParameter<edm::InputTag>("met")))
{
  produces<edm::ValueMap<float> >("met");
  produces<edm::ValueMap<float> >("mt");
}


MuonMetMT::~MuonMetMT()
{
}

void
MuonMetMT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<View<reco::Muon> > probes;
  iEvent.getByToken(probes_, probes);

  Handle<METColl> mets;
  iEvent.getByToken(pfMet_, mets);
  const auto & met = mets->front();

  unsigned int n = probes->size();
  std::vector<float> vmet(n,met.pt());
  std::vector<float> vmt(n,0);

  unsigned int imu = 0;
  for (auto probe = probes->begin(); imu < n; ++probe, ++imu) {
    const reco::Muon &mu = *probe;
    vmt[imu] = std::sqrt(2 * mu.pt() * met.pt() * ( 1 - std::cos(mu.phi()-met.phi()) ) );
  }

  storeMap(iEvent, probes, vmet, "met");
  storeMap(iEvent, probes, vmt,  "mt");
}

template<typename Hand, typename T>
void
MuonMetMT::storeMap(edm::Event &iEvent,
const Hand & handle,
const std::vector<T> & values,
const std::string    & label) const {
  using namespace edm; using namespace std;
  unique_ptr<ValueMap<T> > valMap(new ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valMap), label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMetMT);
