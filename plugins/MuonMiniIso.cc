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

class MuonMiniIso : public edm::EDProducer {
public:

  typedef std::vector< edm::FwdPtr<reco::PFCandidate> > PFCollection;

  explicit MuonMiniIso(const edm::ParameterSet&);
  ~MuonMiniIso();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const edm::EDGetTokenT<edm::View<reco::Muon>> probes_;    
  const edm::EDGetTokenT<PFCollection> pfCandidates_;
  double dRCandProbeVeto_;
  double dRCandSoftActivityCone_;
  double CandPtThreshold_;
  double ChargedPVdZ_;
  bool usePUcands_;

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
MuonMiniIso::MuonMiniIso(const edm::ParameterSet& iConfig):
probes_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("probes"))),
pfCandidates_(consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
dRCandProbeVeto_(iConfig.getParameter<double>("dRCandProbeVeto")),
dRCandSoftActivityCone_(iConfig.getParameter<double>("dRCandSoftActivityCone")),
CandPtThreshold_(iConfig.getParameter<double>("CandPtThreshold"))
{
  produces<edm::ValueMap<float> >("miniIso");
  produces<edm::ValueMap<float> >("activity");
}


MuonMiniIso::~MuonMiniIso()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonMiniIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByToken(probes_, probes);

  //Handle<View<reco::PFCandidate> > pfCandidates;
  Handle<PFCollection> pfCandidates;
  iEvent.getByToken(pfCandidates_, pfCandidates);

  View<reco::Muon>::const_iterator probe, endprobes=probes->end();
  PFCollection::const_iterator iP, beginpf = pfCandidates->begin(), endpf=pfCandidates->end();
  unsigned int n = probes->size();
  
  std::vector<float> iso(n,0);
  std::vector<float> activity(n,0);

  // loop on PROBES
  unsigned int imu = 0;
  for (probe = probes->begin(); probe != endprobes; ++probe, ++imu) {
    const reco::Muon &mu = *probe;

    // loop on PF candidates
    int ipf=0;
    for (iP = beginpf; iP != endpf; ++iP, ++ipf) {
      
      // check pf candidate threshold
      if(iP->get()->pt() < CandPtThreshold_) continue;
      
      double dr = deltaR( *(iP->get() ) , mu );
      // check dr min
      if (dr < dRCandProbeVeto_) continue;

      // miniiso and soft activity definitions
      if (dr > dRCandSoftActivityCone_) continue;
      
      if ( (mu.pt()<=50 && dr <= 0.2) 
        || (50 < mu.pt() && mu.pt() < 200 && dr <= 10/mu.pt())
        || (mu.pt()>200 && dr <= 0.05) ){
        
        iso[imu] += iP->get()->pt();
        
      }else{
        
        activity[imu] += iP->get()->pt();
        
      }
      
    }
  }// end loop on probes

  storeMap(iEvent, probes, iso, "miniIso");
  storeMap(iEvent, probes, activity, "activity");
}

template<typename Hand, typename T>
void
MuonMiniIso::storeMap(edm::Event &iEvent,
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
DEFINE_FWK_MODULE(MuonMiniIso);
