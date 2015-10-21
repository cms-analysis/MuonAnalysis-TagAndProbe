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
  const edm::InputTag probesLabel;

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
  probesLabel(iConfig.getParameter<edm::InputTag>("probes"))
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
  iEvent.getByLabel(probesLabel, probes);
  
  // Rho values
  edm::Handle<double> fixedGridRhoFastjetAllH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetAll", ""), fixedGridRhoFastjetAllH);

  edm::Handle<double> fixedGridRhoFastjetAllCaloH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetAllCalo", ""), fixedGridRhoFastjetAllCaloH);

  edm::Handle<double> fixedGridRhoFastjetCentralCaloH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetCentralCalo", ""), fixedGridRhoFastjetCentralCaloH);

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUpH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetCentralChargedPileUp", ""), fixedGridRhoFastjetCentralChargedPileUpH);
  
  edm::Handle<double> fixedGridRhoFastjetCentralNeutralH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetCentralNeutral", ""), fixedGridRhoFastjetCentralNeutralH);


  // Prepare vector for output 
  std::vector<double> ecalIsoRhoCorr;
  std::vector<double> hcalIsoRhoCorr;
  std::vector<double> combRelIsoRhoCorr;

  View<pat::Muon>::const_iterator probe;

  // Store rho
  std::auto_ptr<ValueMap<float> > fixedGridRhoFastjetAll(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*fixedGridRhoFastjetAll);
  std::vector<float> myfixedGridRhoFastjetAll(probes->size(), *fixedGridRhoFastjetAllH);
  filler1.insert(probes, myfixedGridRhoFastjetAll.begin(), myfixedGridRhoFastjetAll.end());
  filler1.fill();
  iEvent.put(fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll");

  // Store rho
  std::auto_ptr<ValueMap<float> > fixedGridRhoFastjetAllCalo(new ValueMap<float>());
  ValueMap<float>::Filler filler11(*fixedGridRhoFastjetAllCalo);
  std::vector<float> myfixedGridRhoFastjetAllCalo(probes->size(), *fixedGridRhoFastjetAllCaloH);
  filler11.insert(probes, myfixedGridRhoFastjetAllCalo.begin(), myfixedGridRhoFastjetAllCalo.end());
  filler11.fill();
  iEvent.put(fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo");

  // Store rho
  std::auto_ptr<ValueMap<float> > fixedGridRhoFastjetCentralCalo(new ValueMap<float>());
  ValueMap<float>::Filler filler110(*fixedGridRhoFastjetCentralCalo);
  std::vector<float> myfixedGridRhoFastjetCentralCalo(probes->size(), *fixedGridRhoFastjetCentralCaloH);
  filler110.insert(probes, myfixedGridRhoFastjetCentralCalo.begin(), myfixedGridRhoFastjetCentralCalo.end());
  filler110.fill();
  iEvent.put(fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo");

  // Store rho
  std::auto_ptr<ValueMap<float> > fixedGridRhoFastjetCentralChargedPileUp(new ValueMap<float>());
  ValueMap<float>::Filler filler12(*fixedGridRhoFastjetCentralChargedPileUp);
  std::vector<float> myfixedGridRhoFastjetCentralChargedPileUp(probes->size(), *fixedGridRhoFastjetCentralChargedPileUpH);
  filler12.insert(probes, myfixedGridRhoFastjetCentralChargedPileUp.begin(), myfixedGridRhoFastjetCentralChargedPileUp.end());
  filler12.fill();
  iEvent.put(fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp");

  // Store rho
  std::auto_ptr<ValueMap<float> > fixedGridRhoFastjetCentralNeutral(new ValueMap<float>());
  ValueMap<float>::Filler filler13(*fixedGridRhoFastjetCentralNeutral);
  std::vector<float> myfixedGridRhoFastjetCentralNeutral(probes->size(), *fixedGridRhoFastjetCentralNeutralH);
  filler13.insert(probes, myfixedGridRhoFastjetCentralNeutral.begin(), myfixedGridRhoFastjetCentralNeutral.end());
  filler13.fill();
  iEvent.put(fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral");

}

void 
ComputeIsoCorrections::beginJob() {}

void 
ComputeIsoCorrections::endJob() {}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeIsoCorrections);
