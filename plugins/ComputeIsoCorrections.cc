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
  const double effAreaEcalBar;
  const double effAreaEcalEnd;
  const double effAreaHcalBar;
  const double effAreaHcalEnd;

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
  probesLabel(iConfig.getParameter<edm::InputTag>("probes")),
  effAreaEcalBar(iConfig.getParameter<double>("EffAreaEcalBar")),
  effAreaEcalEnd(iConfig.getParameter<double>("EffAreaEcalEnd")),
  effAreaHcalBar(iConfig.getParameter<double>("EffAreaHcalBar")),
  effAreaHcalEnd(iConfig.getParameter<double>("EffAreaHcalEnd"))
{
  produces<edm::ValueMap<float> >("Rho");
  produces<edm::ValueMap<float> >("ecalIsoRhoCorrected");
  produces<edm::ValueMap<float> >("hcalIsoRhoCorrected");
  produces<edm::ValueMap<float> >("combRelIsoRhoCorrected");
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
  
  // Rho value
  edm::Handle<double> rhoH;
  iEvent.getByLabel(InputTag("kt6PFJets", "rho"), rhoH);

  // Prepare vector for output 
  std::vector<double> ecalIsoRhoCorr;
  std::vector<double> hcalIsoRhoCorr;
  std::vector<double> combRelIsoRhoCorr;

  View<pat::Muon>::const_iterator probe, endprobes=probes->end();

  // Loop on probes
  for(probe=probes->begin(); probe!=endprobes; ++probe) {

    // Isolation variables
    double trackIso = probe->trackIso();
    double ecalIso  = probe->ecalIso();
    double hcalIso  = probe->hcalIso();

    if( fabs(probe->eta())<1.479 ) {
      ecalIso=std::max(0., ecalIso-effAreaEcalBar*(*rhoH));
      hcalIso=std::max(0., hcalIso-effAreaHcalBar*(*rhoH));
    }
    else {
      ecalIso=std::max(0., ecalIso-effAreaEcalEnd*(*rhoH));
      hcalIso=std::max(0., hcalIso-effAreaHcalEnd*(*rhoH));
    }
    double probeRhoRelativeIsolation=(trackIso+ecalIso+hcalIso)/std::max(0.5, probe->pt());

    ecalIsoRhoCorr.push_back(ecalIso);
    hcalIsoRhoCorr.push_back(hcalIso);
    combRelIsoRhoCorr.push_back(probeRhoRelativeIsolation);

  } // end loop on probes

  
  // Store rho
  std::auto_ptr<ValueMap<float> > Rho(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*Rho);
  std::vector<float> myRhos(probes->size(), *rhoH);
  filler1.insert(probes, myRhos.begin(), myRhos.end());
  filler1.fill();
  iEvent.put(Rho, "Rho");

  // Store corrected EcalIso 
  std::auto_ptr<ValueMap<float> > EcalIsoRhoCorrected(new ValueMap<float>());
  ValueMap<float>::Filler filler2(*EcalIsoRhoCorrected);
  filler2.insert(probes, ecalIsoRhoCorr.begin(), ecalIsoRhoCorr.end());
  filler2.fill();
  iEvent.put(EcalIsoRhoCorrected, "ecalIsoRhoCorrected");
  
  // Store corrected HcalIso 
  std::auto_ptr<ValueMap<float> > HcalIsoRhoCorrected(new ValueMap<float>());
  ValueMap<float>::Filler filler3(*HcalIsoRhoCorrected);
  filler3.insert(probes, hcalIsoRhoCorr.begin(), hcalIsoRhoCorr.end());
  filler3.fill();
  iEvent.put(HcalIsoRhoCorrected, "hcalIsoRhoCorrected");
  
  // Store corrected CombRelIso 
  std::auto_ptr<ValueMap<float> > RelIsoRhoCorrected(new ValueMap<float>());
  ValueMap<float>::Filler filler4(*RelIsoRhoCorrected);
  filler4.insert(probes, combRelIsoRhoCorr.begin(), combRelIsoRhoCorr.end());
  filler4.fill();
  iEvent.put(RelIsoRhoCorrected, "combRelIsoRhoCorrected");
  
}

void 
ComputeIsoCorrections::beginJob() {}

void 
ComputeIsoCorrections::endJob() {}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeIsoCorrections);
