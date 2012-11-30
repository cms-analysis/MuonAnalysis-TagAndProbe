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

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class declaration
//

class ComputeL1HLTPrescales : public edm::EDProducer {
public:
  explicit ComputeL1HLTPrescales(const edm::ParameterSet&);
  ~ComputeL1HLTPrescales();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginRun(edm::Run&, edm::EventSetup const&);

  // ----------member data ---------------------------
  const edm::InputTag probesLabel_;
  const std::string hltConfigLabel_;
  std::vector<std::string> hltPaths_;
  HLTConfigProvider hltConfig_;

  void writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const float value, const std::string &label) ;
};

ComputeL1HLTPrescales::ComputeL1HLTPrescales(const edm::ParameterSet& iConfig):
  probesLabel_(iConfig.getParameter<edm::InputTag>("probes")),
  hltConfigLabel_(iConfig.getParameter<std::string>("hltConfig")),
  hltPaths_(iConfig.getParameter<std::vector<std::string> >("hltPaths"))
{
  for(unsigned int j=0; j<hltPaths_.size(); ++j) {
    // Trigger name should include "_v" at the end
    std::string tmpStr = hltPaths_[j].substr(0, hltPaths_[j].find("_v")); 
    while(tmpStr.find("_")!=std::string::npos) {tmpStr.replace(tmpStr.find("_"),1,"");}
    produces<edm::ValueMap<float> >( (tmpStr+"TotalPrescale").c_str() );
  }
}


ComputeL1HLTPrescales::~ComputeL1HLTPrescales() {}


void ComputeL1HLTPrescales::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Read input
  Handle<View<reco::Candidate> > probes;
  iEvent.getByLabel(probesLabel_, probes);

  if(!hltConfig_.inited()) {
    for(unsigned int j=0; j<hltPaths_.size(); ++j) {
      std::string tmpStr = hltPaths_[j].substr(0, hltPaths_[j].find("_v")); 
      while(tmpStr.find("_")!=std::string::npos) {tmpStr.replace(tmpStr.find("_"),1,"");}
      writeGlobalFloat(iEvent, probes, -2., tmpStr+"TotalPrescale");
    }
    return;
  }

  for(unsigned int j=0; j<hltPaths_.size(); ++j) {
    std::string trigName = "";

    for(unsigned int iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
      if( hltConfig_.triggerName(iHltPath).find(hltPaths_[j]) != std::string::npos ) {
	trigName = hltConfig_.triggerName(iHltPath);
	break;
      }
    } // end for(iHltPath)

    float totPs = -1.;
    if(trigName.length()>0) {
      std::pair<int,int> pss = hltConfig_.prescaleValues(iEvent, iSetup, trigName);
      totPs = pss.first * pss.second;
    }

    // Works if trigger name does not contain "_v" in the middle... 
    // Otherwise use HLTConfigProvider::removeVersion() or regexp "_v[0-9]+$"
    trigName = trigName.substr(0, trigName.find("_v")); 
    while(trigName.find("_")!=std::string::npos) {trigName.replace(trigName.find("_"),1,"");}
    writeGlobalFloat(iEvent, probes, totPs, trigName+"TotalPrescale");
  } // end for(j)
}

void ComputeL1HLTPrescales::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if(!hltConfig_.init(iRun, iSetup, hltConfigLabel_, changed)) 
    std::cout << "Warning, didn't find HLTConfigProvider with label " 
	      << hltConfigLabel_.c_str() << " in run " << iRun.run() << std::endl;
}

void ComputeL1HLTPrescales::writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const float value, const std::string &label) { 
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  std::vector<float> values(probes->size(), value);
  filler.insert(probes, values.begin(), values.end());
  filler.fill();
  iEvent.put(out, label);
}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeL1HLTPrescales);
