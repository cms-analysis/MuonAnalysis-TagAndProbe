//
// $Id: GenAdditionalInfo.cc,v 1.1 2015/08/23 13:49:16 battilan Exp $
//

/**
  \class    GenAdditionalInfo 
  \brief    Adds generator level information to the T&P according to 74X interface.
            https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf
            
  \Author   Carlo Battilana
  \version  $Id: GenAdditionalInfo.cc,v 1.1 2015/08/23 13:49:16 battilan Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

class GenAdditionalInfo : public edm::EDProducer
{
  
public:

  explicit GenAdditionalInfo(const edm::ParameterSet & iConfig);
  virtual ~GenAdditionalInfo() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) final;

private:

  edm::EDGetTokenT<GenEventInfoProduct> m_genInfoTag;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> m_pileUpInfoTag;
  edm::EDGetTokenT<edm::View<reco::Candidate>> m_pairTag;

  /// Write a ValueMap<float> in the event
  template<typename T> void writeValueMap(edm::Event &ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
					  const std::vector<T> & values, const std::string    & label) const ;

};


GenAdditionalInfo::GenAdditionalInfo(const edm::ParameterSet & iConfig) :
  m_genInfoTag(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoTag"))),
  m_pileUpInfoTag(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileUpInfoTag"))),
  m_pairTag(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairTag")))
{
  
  produces<edm::ValueMap<float> >("genWeight");
  produces<edm::ValueMap<float> >("truePileUp");
  produces<edm::ValueMap<float> >("actualPileUp");

}

void GenAdditionalInfo::produce(edm::Event & ev, const edm::EventSetup & iSetup)
{

  float weight = 1.;

  float truePU   = -1.;
  float actualPU = -1.;

  if (!ev.isRealData())
    {
      edm::Handle<GenEventInfoProduct> genInfo;
      ev.getByToken(m_genInfoTag, genInfo);

      edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
      ev.getByToken(m_pileUpInfoTag, puInfo);

      weight = genInfo->weight();

      for( auto & puInfoEntry : *puInfo.product() ) 
	{
	  int bx = puInfoEntry.getBunchCrossing();
	  
	  if(bx == 0) 
	    { 
	      truePU   = puInfoEntry.getTrueNumInteractions();
	      actualPU = puInfoEntry.getPU_NumInteractions();
	      continue;
	    }
	}

    }

  edm::Handle<edm::View<reco::Candidate> > pairs;
  ev.getByToken(m_pairTag, pairs);

  size_t n = pairs->size();
  std::vector<float> genWeight(n,0);
  std::vector<float>   truePileUp(n,0);
  std::vector<float>   actualPileUp(n,0);
  
  for (size_t iPair = 0; iPair < n; ++iPair)
    {
      const reco::Candidate & pair = (*pairs)[iPair];
      if (pair.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << 
	"[GenAdditionalInfo::produce] GenAdditionalInfo should be used on composite candidates with two daughters, this one has " << pair.numberOfDaughters() << "\n";

      genWeight[iPair]    = weight;
      truePileUp[iPair]   = truePU;
      actualPileUp[iPair] = actualPU;

    }

  writeValueMap<float>(ev, pairs, genWeight,    "genWeight");
  writeValueMap<float>(ev, pairs, truePileUp,   "truePileUp");
  writeValueMap<float>(ev, pairs, actualPileUp, "actualPileUp");

}

template<typename T> void GenAdditionalInfo::writeValueMap(edm::Event & ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
						       const std::vector<T> & values, const std::string & label) const 
{
  
  std::unique_ptr<edm::ValueMap<T> > valMap(new edm::ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  ev.put(std::move(valMap), label);
  
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenAdditionalInfo);
