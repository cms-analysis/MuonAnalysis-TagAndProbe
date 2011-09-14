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

//
// class declaration
//

class MuonSIP : public edm::EDProducer {
public:
  explicit MuonSIP(const edm::ParameterSet&);
  ~MuonSIP();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
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
MuonSIP::MuonSIP(const edm::ParameterSet& iConfig):
probes_(iConfig.getParameter<edm::InputTag>("probes"))

{
  produces<edm::ValueMap<float> >("IP");
  produces<edm::ValueMap<float> >("IPError");
  produces<edm::ValueMap<float> >("SIP");
}


MuonSIP::~MuonSIP()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonSIP::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<pat::Muon> > probes;
  iEvent.getByLabel(probes_, probes);

  // prepare vector for output    
  std::vector<double> muon_IP;
  std::vector<double> muon_IPError;
  std::vector<double> muon_SIP;

  View<pat::Muon>::const_iterator probe, endprobes=probes->end();


  // loop on PROBES
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    Double_t IP      = fabs(probe->dB(pat::Muon::PV3D));
    Double_t IPError = probe->edB(pat::Muon::PV3D);	  
    Double_t SIP     = IP/IPError;
    
    muon_IP.push_back(IP);
    muon_IPError.push_back(IPError);
    muon_SIP.push_back(SIP);
  
  }// end loop on probes

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > IP(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > IPError(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > SIP(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*IP);
  ValueMap<float>::Filler filler2(*IPError);
  ValueMap<float>::Filler filler3(*SIP);
  filler1.insert(probes, muon_IP.begin(), muon_IP.end());
  filler2.insert(probes, muon_IPError.begin(), muon_IPError.end());
  filler3.insert(probes, muon_SIP.begin(), muon_SIP.end());
  filler1.fill();
  filler2.fill();
  filler3.fill();

  iEvent.put(IP, "IP");
  iEvent.put(IPError, "IPError");
  iEvent.put(SIP, "SIP");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonSIP::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonSIP::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonSIP);
