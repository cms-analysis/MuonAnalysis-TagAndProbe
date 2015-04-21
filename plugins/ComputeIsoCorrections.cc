//
// ComputeIsoCorrections.cc,  2011/03/3 
//

/**
  \brief    

  \author   Clara Jorda (adapted and ruined by D. Trocino)
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

  produces<edm::ValueMap<float> >("RhoAllCalo");
  produces<edm::ValueMap<float> >("RhoAll");
  produces<edm::ValueMap<float> >("RhoPU");
  produces<edm::ValueMap<float> >("RhoNeu05");
  produces<edm::ValueMap<float> >("RhoNeu1");
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
  edm::Handle<double> rhoAllH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetAll", ""), rhoAllH);

  edm::Handle<double> rhoAllCaloH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetAllCalo", ""), rhoAllCaloH);

  edm::Handle<double> rhoPUH;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetCentralChargedPileUp", ""), rhoPUH);
  
  edm::Handle<double> rhoNeu05H;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetCentralNeutral", ""), rhoNeu05H);
  
  edm::Handle<double> rhoNeu1H;
  iEvent.getByLabel(InputTag("fixedGridRhoFastjetAll", ""), rhoNeu1H); //MM FIXME, was tight in 53X


  // Prepare vector for output 
  std::vector<double> ecalIsoRhoCorr;
  std::vector<double> hcalIsoRhoCorr;
  std::vector<double> combRelIsoRhoCorr;

  View<pat::Muon>::const_iterator probe;

  // Store rho
  std::auto_ptr<ValueMap<float> > RhoAll(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*RhoAll);
  std::vector<float> myRhosAll(probes->size(), *rhoAllH);
  filler1.insert(probes, myRhosAll.begin(), myRhosAll.end());
  filler1.fill();
  iEvent.put(RhoAll, "RhoAll");

  // Store rho
  std::auto_ptr<ValueMap<float> > RhoAllCalo(new ValueMap<float>());
  ValueMap<float>::Filler filler11(*RhoAllCalo);
  std::vector<float> myRhosAllCalo(probes->size(), *rhoAllCaloH);
  filler11.insert(probes, myRhosAllCalo.begin(), myRhosAllCalo.end());
  filler11.fill();
  iEvent.put(RhoAllCalo, "RhoAllCalo");

  // Store rho
  std::auto_ptr<ValueMap<float> > RhoPU(new ValueMap<float>());
  ValueMap<float>::Filler filler12(*RhoPU);
  std::vector<float> myRhosPU(probes->size(), *rhoPUH);
  filler12.insert(probes, myRhosPU.begin(), myRhosPU.end());
  filler12.fill();
  iEvent.put(RhoPU, "RhoPU");

  // Store rho
  std::auto_ptr<ValueMap<float> > RhoNeu05(new ValueMap<float>());
  ValueMap<float>::Filler filler13(*RhoNeu05);
  std::vector<float> myRhosNeu05(probes->size(), *rhoNeu05H);
  filler13.insert(probes, myRhosNeu05.begin(), myRhosNeu05.end());
  filler13.fill();
  iEvent.put(RhoNeu05, "RhoNeu05");

  // Store rho
  std::auto_ptr<ValueMap<float> > RhoNeu1(new ValueMap<float>());
  ValueMap<float>::Filler filler14(*RhoNeu1);
  std::vector<float> myRhosNeu1(probes->size(), *rhoNeu1H);
  filler14.insert(probes, myRhosNeu1.begin(), myRhosNeu1.end());
  filler14.fill();
  iEvent.put(RhoNeu1, "RhoNeu1");
}

void 
ComputeIsoCorrections::beginJob() {}

void 
ComputeIsoCorrections::endJob() {}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeIsoCorrections);
