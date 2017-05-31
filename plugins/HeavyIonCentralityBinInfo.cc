// -*- C++ -*-
//
// Package:    MuonAnalysis/TagAndProbe
// Class:      HeavyIonCentralityBinInfo
// 
/**\class HeavyIonCentralityBinInfo HeavyIonCentralityBinInfo.cc MuonAnalysis/TagAndProbe/plugins/HeavyIonCentralityBinInfo.cc

 Description: Add centrality bin information on the tag muons
*/
//
// Original Author:  Anna Zsigmond
//         Created:  Mon, 01 Feb 2016 13:03:29 GMT
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

class HeavyIonCentralityBinInfo : public edm::stream::EDProducer<> {
   public:
      explicit HeavyIonCentralityBinInfo(const edm::ParameterSet&);
      ~HeavyIonCentralityBinInfo();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<edm::View<reco::Candidate> > src_;
      edm::EDGetTokenT<int> CentralityBinTag_;
};

HeavyIonCentralityBinInfo::HeavyIonCentralityBinInfo(const edm::ParameterSet& iConfig) :
   src_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"))),
   CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc")))
{
   produces<edm::ValueMap<float> >();
}


HeavyIonCentralityBinInfo::~HeavyIonCentralityBinInfo()
{
}

void
HeavyIonCentralityBinInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<View<reco::Candidate> > src;
   iEvent.getByToken(src_, src);

   Handle<int> cbin_;
   iEvent.getByToken(CentralityBinTag_,cbin_);
   int hiBin = *cbin_;

   vector<float> values(src->size(), (float)hiBin);

   unique_ptr<ValueMap<float> > valMap(new ValueMap<float>());
   ValueMap<float>::Filler filler(*valMap);
   filler.insert(src, values.begin(), values.end());
   filler.fill();
   iEvent.put(std::move(valMap));
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyIonCentralityBinInfo);
