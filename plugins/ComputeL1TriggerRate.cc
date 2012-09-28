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

#include "DataFormats/Scalers/interface/Level1TriggerScalers.h"

//
// class declaration
//

class ComputeL1TriggerRate : public edm::EDProducer {
public:
  explicit ComputeL1TriggerRate(const edm::ParameterSet&);
  ~ComputeL1TriggerRate();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const edm::InputTag probesLabel_;
  const edm::InputTag scalersLabel_;

  void writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) ;
};

ComputeL1TriggerRate::ComputeL1TriggerRate(const edm::ParameterSet& iConfig):
  probesLabel_(iConfig.getParameter<edm::InputTag>("probes")),
  scalersLabel_(iConfig.getParameter<edm::InputTag>("scalers"))
{
  produces<edm::ValueMap<float> >();
  produces<edm::ValueMap<float> >("bx");
}


ComputeL1TriggerRate::~ComputeL1TriggerRate() {}


void ComputeL1TriggerRate::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Read input
  Handle<View<reco::Candidate> > probes;
  iEvent.getByLabel(probesLabel_, probes);

  edm::Handle<std::vector<Level1TriggerScalers> > scalers; 
  iEvent.getByLabel(scalersLabel_, scalers); 

  writeGlobalFloat(iEvent, probes, scalers->empty() ? -1 : scalers->front().gtTriggersRate(), "");
  writeGlobalFloat(iEvent, probes, iEvent.bunchCrossing(), "bx");
}

void ComputeL1TriggerRate::writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) { 
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  std::vector<float> values(probes->size(), value);
  filler.insert(probes, values.begin(), values.end());
  filler.fill();
  iEvent.put(out, label);
}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeL1TriggerRate);
