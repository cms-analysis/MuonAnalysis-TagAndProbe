// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
//
// class declaration
//

class ComputeTkClusterInfos : public edm::EDProducer {
public:
  explicit ComputeTkClusterInfos(const edm::ParameterSet&);
  ~ComputeTkClusterInfos();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const edm::InputTag probesLabel_;
  const edm::InputTag stripLabel_, pixelLabel_;

  void writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) ;
};

ComputeTkClusterInfos::ComputeTkClusterInfos(const edm::ParameterSet& iConfig):
  probesLabel_(iConfig.getParameter<edm::InputTag>("probes")),
  stripLabel_(iConfig.getParameter<edm::InputTag>("stripClusters")),
  pixelLabel_(iConfig.getParameter<edm::InputTag>("pixelClusters"))
{

  produces<edm::ValueMap<float> >("siStripClusterCount");
  produces<edm::ValueMap<float> >("siPixelClusterCount");
}


ComputeTkClusterInfos::~ComputeTkClusterInfos() {}


void ComputeTkClusterInfos::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Read input
  Handle<View<reco::Candidate> > probes;
  iEvent.getByLabel(probesLabel_, probes);

  edm::Handle<edmNew::DetSetVector<SiStripCluster> > stripC; 
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelC; 
  iEvent.getByLabel(stripLabel_, stripC); 
  iEvent.getByLabel(pixelLabel_, pixelC); 

  writeGlobalFloat(iEvent, probes, stripC->dataSize(), "siStripClusterCount");
  writeGlobalFloat(iEvent, probes, pixelC->dataSize(), "siPixelClusterCount");
}

void ComputeTkClusterInfos::writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) { 
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  std::vector<float> values(probes->size(), value);
  filler.insert(probes, values.begin(), values.end());
  filler.fill();
  iEvent.put(out, label);
}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeTkClusterInfos);
