#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/Common/interface/ValueMap.h"

class MuonPFIsoSingleTypeMapProd : public edm::EDProducer {
public:
  explicit MuonPFIsoSingleTypeMapProd(const edm::ParameterSet&);
  ~MuonPFIsoSingleTypeMapProd() ; 

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag muonLabel_;
  edm::InputTag pfLabel_;
  StringCutObjectSelector<reco::PFCandidate> pfSelection_;
  double deltaR_;

};



MuonPFIsoSingleTypeMapProd::MuonPFIsoSingleTypeMapProd(const edm::ParameterSet& iConfig) :
  muonLabel_(iConfig.getParameter<edm::InputTag>("muonLabel")),
  pfLabel_(iConfig.getParameter<edm::InputTag>("pfLabel")),
  pfSelection_(iConfig.getParameter<std::string>("pfSelection")),
  deltaR_(iConfig.getParameter<double>("deltaR"))
{
  produces<edm::ValueMap<float> >();
}

void MuonPFIsoSingleTypeMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon> > muH;
  iEvent.getByLabel(muonLabel_, muH);

  edm::Handle<edm::View<reco::PFCandidate> > pfH;
  iEvent.getByLabel(pfLabel_, pfH);

  std::vector<float> isoV;
  std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler isoF(*isoM);

  for(size_t i=0; i<muH->size();++i) {
      const reco::Muon &mu = muH->at(i);
      double ptSum = 0;

      for (size_t j=0; j<pfH->size();j++) {   
          const reco::PFCandidate &pf = pfH->at(j);

          double dr = deltaR(mu, pf);

          // In cone
          if (dr >= deltaR_) continue;            

          // Passes selection
          if (!pfSelection_(pf)) continue;

          // exclude self
          if(pf.trackRef().isNonnull() && mu.track().isNonnull() &&
                  pf.trackRef() == mu.track()) continue;

          // scalar sum
          ptSum += pf.pt();

      }
      isoV.push_back(ptSum);
  }

  isoF.insert(muH,isoV.begin(),isoV.end());

  isoF.fill();
  iEvent.put(isoM);

}

MuonPFIsoSingleTypeMapProd::~MuonPFIsoSingleTypeMapProd() { }
DEFINE_FWK_MODULE(MuonPFIsoSingleTypeMapProd);
