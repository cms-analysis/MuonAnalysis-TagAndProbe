//
// $Id: AdditionalEventInfo.cc,v 1.1 2015/08/23 13:49:16 battilan Exp $
//

/**
  \class    AdditionalEventInfo 
  \brief    Adds bx, orbit and inst lumi info (the last from scalers)
            
  \Author   Carlo Battilana
  \version  $Id: AdditionalEventInfo.cc,v 1.1 2015/01/06 10:41:10 battilan Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"




class AdditionalEventInfo : public edm::EDProducer
{
  
public:

  explicit AdditionalEventInfo(const edm::ParameterSet & iConfig);
  virtual ~AdditionalEventInfo() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) final;

private:

  edm::EDGetTokenT<LumiScalersCollection> m_lumiScalerTag;
  edm::EDGetTokenT<edm::View<reco::Candidate> > m_muonTag;
  
  /// Write a ValueMap<T> in the event
  template<typename T> void writeValueMap(edm::Event &ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
					  const std::vector<T> & values, const std::string    & label) const ;
  
};

AdditionalEventInfo::AdditionalEventInfo(const edm::ParameterSet & iConfig) :
  m_lumiScalerTag(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalerTag"))),
  m_muonTag(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("muonTag")))
{
  
  // produces<edm::ValueMap<ULong64_t> >("orbit");
  produces<edm::ValueMap<float> >("instLumi");

}

void AdditionalEventInfo::produce(edm::Event & ev, const edm::EventSetup & iSetup)
{

  // ULong64_t orbit = 0;
  float instLumi = -999.;

  if (ev.isRealData())
    {
      // orbit = ev.orbitNumber();
      
      if (!m_lumiScalerTag.isUninitialized())
	{
	  edm::Handle<LumiScalersCollection> lumiScaler;
	  ev.getByToken(m_lumiScalerTag, lumiScaler);
	  
	  if (lumiScaler->begin() != lumiScaler->end())
	    instLumi= lumiScaler->begin()->instantLumi();
	} 
      else
	{
	  throw cms::Exception("CorruptData") << 
	    "[AdditionalEventInfo::produce] AdditionalEventInfo requires a valid LumiScalerCollecion InpuTag" << std::endl;
	}
    } 
  
  edm::Handle<edm::View<reco::Candidate> > muons;
  ev.getByToken(m_muonTag, muons);

  size_t n = muons->size();
  
  // std::vector<ULong64_t> orbits(n, orbit);
  std::vector<float> instLumis(n, instLumi);
  

  // writeValueMap<ULong64_t>(ev, muons, orbits, "orbit");
  writeValueMap<float>(ev, muons, instLumis, "instLumi");

}

template<typename T> void AdditionalEventInfo::writeValueMap(edm::Event & ev, const edm::Handle<edm::View<reco::Candidate> > & handle,
							     const std::vector<T> & values, const std::string & label) const 
{
  
  std::unique_ptr<edm::ValueMap<T> > valMap(new edm::ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  ev.put(std::move(valMap), label);
  
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AdditionalEventInfo);
