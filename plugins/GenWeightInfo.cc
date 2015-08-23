//
// $Id: GenWeightInfo.cc,v 1.1 2015/08/23 13:49:16 battilan Exp $
//

/**
  \class    GenWeightInfo 
  \brief    Adds generator level information to the T&P according to 74X interface.
            https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf
            
  \Author   Carlo Battilana
  \version  $Id: GenWeightInfo.cc,v 1.1 2015/08/23 13:49:16 battilan Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

class GenWeightInfo : public edm::EDProducer
{
  
public:

  explicit GenWeightInfo(const edm::ParameterSet & iConfig);
  virtual ~GenWeightInfo() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) final;

private:

  edm::InputTag m_genInfoTag;
  edm::InputTag m_pairTag;

  /// Write a ValueMap<float> in the event
  void writeValueMap(edm::Event &ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
		     const std::vector<float> & values, const std::string    & label) const ;

};


GenWeightInfo::GenWeightInfo(const edm::ParameterSet & iConfig) :
  m_genInfoTag(iConfig.getParameter<edm::InputTag>("genInfoTag")),
  m_pairTag(iConfig.getParameter<edm::InputTag>("pairTag"))
{
  
  produces<edm::ValueMap<float> >("genWeight");

}

void GenWeightInfo::produce(edm::Event & ev, const edm::EventSetup & iSetup)
{

  float weight = 1.;

  if (!ev.isRealData())
    {
      edm::Handle<GenEventInfoProduct> genInfo;
      ev.getByLabel(m_genInfoTag, genInfo);

      weight = genInfo->weight();
    } 
  
  edm::Handle<edm::View<reco::Candidate> > pairs;
  ev.getByLabel(m_pairTag, pairs);

  size_t n = pairs->size();
  std::vector<float> genWeight(n,0);
  
  for (size_t iPair = 0; iPair < n; ++iPair)
    {
      const reco::Candidate & pair = (*pairs)[iPair];
      if (pair.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << 
	"[GenWeightInfo::produce] GenWeightInfo should be used on composite candidates with two daughters, this one has " << pair.numberOfDaughters() << "\n";

      genWeight[iPair] = weight;
    }

  writeValueMap(ev, pairs, genWeight, "genWeight");

}

void GenWeightInfo::writeValueMap(edm::Event & ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
				  const std::vector<float> & values, const std::string & label) const 
{
  
  std::auto_ptr<edm::ValueMap<float> > valMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  ev.put(valMap, label);
  
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenWeightInfo);
