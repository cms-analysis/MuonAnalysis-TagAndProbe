//
// ComputeIsoCorrections.cc,  2011/03/3 
//

/**
  \brief    

  \author   Clara Jorda (adapted and ruined by D. Trocino), Luca Perrozzi
*/


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

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <Math/VectorUtil.h>
#include "TMath.h"

//
// class declaration
//

class ComputeIsoCorrections : public edm::EDProducer {
public:
  explicit ComputeIsoCorrections(const edm::ParameterSet&);
  ~ComputeIsoCorrections();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::EDGetTokenT<edm::View<pat::Muon>> probesLabel;
  edm::EDGetTokenT<double> fixedGridRhoFastjetAllCalo_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetAll_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralChargedPileUp_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralNeutral_;
  edm::EDGetTokenT<double> fixedGridRhoFastjetCentralCalo_;
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
ComputeIsoCorrections::ComputeIsoCorrections(const edm::ParameterSet& iConfig):
  probesLabel(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("probes"))),
  fixedGridRhoFastjetAllCalo_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAllCalo"))),
  fixedGridRhoFastjetAll_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  fixedGridRhoFastjetCentralChargedPileUp_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"))),
  fixedGridRhoFastjetCentralNeutral_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"))),
  fixedGridRhoFastjetCentralCalo_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralCalo"))) 
{

  produces<edm::ValueMap<float> >("fixedGridRhoFastjetAllCalo");
  produces<edm::ValueMap<float> >("fixedGridRhoFastjetAll");
  produces<edm::ValueMap<float> >("fixedGridRhoFastjetCentralChargedPileUp");
  produces<edm::ValueMap<float> >("fixedGridRhoFastjetCentralNeutral");
  produces<edm::ValueMap<float> >("fixedGridRhoFastjetCentralCalo");
  
}


ComputeIsoCorrections::~ComputeIsoCorrections() {}


//
// member functions
//

void ComputeIsoCorrections::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Read input
  Handle<View<pat::Muon> > probes;
  iEvent.getByToken(probesLabel, probes);
  
  // Rho values
  edm::Handle<double> fixedGridRhoFastjetAllH;
  iEvent.getByToken(fixedGridRhoFastjetAll_, fixedGridRhoFastjetAllH);

  edm::Handle<double> fixedGridRhoFastjetAllCaloH;
  iEvent.getByToken(fixedGridRhoFastjetAllCalo_, fixedGridRhoFastjetAllCaloH);

  edm::Handle<double> fixedGridRhoFastjetCentralCaloH;
  iEvent.getByToken(fixedGridRhoFastjetCentralCalo_, fixedGridRhoFastjetCentralCaloH);

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUpH;
  iEvent.getByToken(fixedGridRhoFastjetCentralChargedPileUp_, fixedGridRhoFastjetCentralChargedPileUpH);
  
  edm::Handle<double> fixedGridRhoFastjetCentralNeutralH;
  iEvent.getByToken(fixedGridRhoFastjetCentralNeutral_, fixedGridRhoFastjetCentralNeutralH);


  // Prepare vector for output 
  std::vector<double> ecalIsoRhoCorr;
  std::vector<double> hcalIsoRhoCorr;
  std::vector<double> combRelIsoRhoCorr;

  View<pat::Muon>::const_iterator probe;

  // Store rho
  std::unique_ptr<ValueMap<float> > fixedGridRhoFastjetAll(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*fixedGridRhoFastjetAll);
  std::vector<float> myfixedGridRhoFastjetAll(probes->size(), *fixedGridRhoFastjetAllH);
  filler1.insert(probes, myfixedGridRhoFastjetAll.begin(), myfixedGridRhoFastjetAll.end());
  filler1.fill();
  iEvent.put(std::move(fixedGridRhoFastjetAll), "fixedGridRhoFastjetAll");

  // Store rho
  std::unique_ptr<ValueMap<float> > fixedGridRhoFastjetAllCalo(new ValueMap<float>());
  ValueMap<float>::Filler filler11(*fixedGridRhoFastjetAllCalo);
  std::vector<float> myfixedGridRhoFastjetAllCalo(probes->size(), *fixedGridRhoFastjetAllCaloH);
  filler11.insert(probes, myfixedGridRhoFastjetAllCalo.begin(), myfixedGridRhoFastjetAllCalo.end());
  filler11.fill();
  iEvent.put(std::move(fixedGridRhoFastjetAllCalo), "fixedGridRhoFastjetAllCalo");

  // Store rho
  std::unique_ptr<ValueMap<float> > fixedGridRhoFastjetCentralCalo(new ValueMap<float>());
  ValueMap<float>::Filler filler110(*fixedGridRhoFastjetCentralCalo);
  std::vector<float> myfixedGridRhoFastjetCentralCalo(probes->size(), *fixedGridRhoFastjetCentralCaloH);
  filler110.insert(probes, myfixedGridRhoFastjetCentralCalo.begin(), myfixedGridRhoFastjetCentralCalo.end());
  filler110.fill();
  iEvent.put(std::move(fixedGridRhoFastjetCentralCalo), "fixedGridRhoFastjetCentralCalo");

  // Store rho
  std::unique_ptr<ValueMap<float> > fixedGridRhoFastjetCentralChargedPileUp(new ValueMap<float>());
  ValueMap<float>::Filler filler12(*fixedGridRhoFastjetCentralChargedPileUp);
  std::vector<float> myfixedGridRhoFastjetCentralChargedPileUp(probes->size(), *fixedGridRhoFastjetCentralChargedPileUpH);
  filler12.insert(probes, myfixedGridRhoFastjetCentralChargedPileUp.begin(), myfixedGridRhoFastjetCentralChargedPileUp.end());
  filler12.fill();
  iEvent.put(std::move(fixedGridRhoFastjetCentralChargedPileUp), "fixedGridRhoFastjetCentralChargedPileUp");

  // Store rho
  std::unique_ptr<ValueMap<float> > fixedGridRhoFastjetCentralNeutral(new ValueMap<float>());
  ValueMap<float>::Filler filler13(*fixedGridRhoFastjetCentralNeutral);
  std::vector<float> myfixedGridRhoFastjetCentralNeutral(probes->size(), *fixedGridRhoFastjetCentralNeutralH);
  filler13.insert(probes, myfixedGridRhoFastjetCentralNeutral.begin(), myfixedGridRhoFastjetCentralNeutral.end());
  filler13.fill();
  iEvent.put(std::move(fixedGridRhoFastjetCentralNeutral), "fixedGridRhoFastjetCentralNeutral");

}

void 
ComputeIsoCorrections::beginJob() {}

void 
ComputeIsoCorrections::endJob() {}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeIsoCorrections);
