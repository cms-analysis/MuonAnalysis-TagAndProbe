// system include files
#include <memory>
#include <cmath>
#include <set>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/MessageLogger/interface/ELseverityLevel.h"
#include "FWCore/MessageLogger/interface/ELstring.h"
#include "FWCore/MessageLogger/interface/ErrorSummaryEntry.h"

//
// class declaration
//

class ComputeLogErrorTotals : public edm::EDProducer {
public:
  explicit ComputeLogErrorTotals(const edm::ParameterSet&);
  ~ComputeLogErrorTotals();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  edm::InputTag probesLabel_;
  edm::InputTag logErrorLabel_;

  std::map<std::string, std::set<std::string> > counters_;

  void writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) ;
};

ComputeLogErrorTotals::ComputeLogErrorTotals(const edm::ParameterSet& iConfig):
  probesLabel_(iConfig.getParameter<edm::InputTag>("probes")),
  logErrorLabel_(iConfig.getParameter<edm::InputTag>("logErrors"))
{
    edm::ParameterSet idps = iConfig.getParameter<edm::ParameterSet>("counters");
    std::vector<std::string> names = idps.getParameterNamesForType<std::vector<std::string> >();
    for (std::vector<std::string>::const_iterator it = names.begin(), ed = names.end(); it != ed; ++it) {
        std::vector<std::string> modules = idps.getParameter<std::vector<std::string> >(*it);
        counters_[*it] = std::set<std::string>(modules.begin(), modules.end());
        produces<edm::ValueMap<float> >(*it);
    }
}


ComputeLogErrorTotals::~ComputeLogErrorTotals() {}


void ComputeLogErrorTotals::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Read input
  Handle<View<reco::Candidate> > probes;
  iEvent.getByLabel(probesLabel_, probes);

  edm::Handle<std::vector<edm::ErrorSummaryEntry> >  errors;
  iEvent.getByLabel(logErrorLabel_,errors);

  // Compare severity level of error with ELseveritylevel instance el : "-e" should be the lowest error
  edm::ELseverityLevel el("-e");

  for (std::map<std::string,std::set<std::string> >::const_iterator it = counters_.begin(), ed = counters_.end(); it != ed; ++it) {
      int nerr = 0;
      for( size_t i = 0, n = errors->size();  i < n ; i++){    
          const edm::ErrorSummaryEntry &err = (*errors)[i];
          if(err.severity.getLevel() < el.getLevel()) continue;
          std::string s = err.module;
          size_t pos = s.find(':'); if (pos == std::string::npos) continue;
          s = s.substr(pos+1,s.size());
          if (it->second.count(s)) nerr += err.count;
      } 
      if (nerr) std::cout << " found " << nerr << " for " << it->first << " (" << it->second.size() << " modules)" << std::endl;
      writeGlobalFloat(iEvent, probes, float(nerr), it->first);
  }
}

void ComputeLogErrorTotals::writeGlobalFloat(edm::Event &iEvent, const edm::Handle<edm::View<reco::Candidate> > &probes, const double value, const std::string &label) { 
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  std::vector<float> values(probes->size(), value);
  filler.insert(probes, values.begin(), values.end());
  filler.fill();
  iEvent.put(out, label);
}

// Define this module as plugin
DEFINE_FWK_MODULE(ComputeLogErrorTotals);
