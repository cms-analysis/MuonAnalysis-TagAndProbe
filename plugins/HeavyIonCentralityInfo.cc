// -*- C++ -*-
//
// Package:    MuonAnalysis/TagAndProbe
// Class:      HeavyIonCentralityInfo
// 
/**\class HeavyIonCentralityInfo HeavyIonCentralityInfo.cc MuonAnalysis/TagAndProbe/plugins/HeavyIonCentralityInfo.cc

 Description: Add centrality variables information on the tag muons
*/
//
// Original Author:  Anna Zsigmond
//         Created:  Mon, 01 Feb 2016 14:45:12 GMT
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

class HeavyIonCentralityInfo : public edm::stream::EDProducer<> {
   public:
      explicit HeavyIonCentralityInfo(const edm::ParameterSet&);
      ~HeavyIonCentralityInfo();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<edm::View<reco::Candidate> > src_;
      edm::EDGetTokenT<reco::Centrality> CentralityTag_;

      void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Candidate> > & handle,
                const std::vector<float> & values,
                const std::string    & label) const ;
};

HeavyIonCentralityInfo::HeavyIonCentralityInfo(const edm::ParameterSet& iConfig) :
   src_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"))),
   CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc")))
{
   produces<edm::ValueMap<float> >("HFtowers");
   produces<edm::ValueMap<float> >("HFtowersPlus");
   produces<edm::ValueMap<float> >("HFtowersMinus");
   produces<edm::ValueMap<float> >("HFtowersTrunc");
   produces<edm::ValueMap<float> >("HFtowersPlusTrunc");
   produces<edm::ValueMap<float> >("HFtowersMinusTrunc");
   produces<edm::ValueMap<float> >("HFhits");
   produces<edm::ValueMap<float> >("PixelHits");
   produces<edm::ValueMap<float> >("PixelTracks");
   produces<edm::ValueMap<float> >("Tracks");
   produces<edm::ValueMap<float> >("EB");
   produces<edm::ValueMap<float> >("EE");
   produces<edm::ValueMap<float> >("ET");
}


HeavyIonCentralityInfo::~HeavyIonCentralityInfo()
{
}

void
HeavyIonCentralityInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<View<reco::Candidate> > src;
   iEvent.getByToken(src_, src);

   Handle<reco::Centrality> centrality;
   iEvent.getByToken(CentralityTag_, centrality);

   vector<float> hiHF(src->size(), (float)centrality->EtHFtowerSum());
   vector<float> hiHFplus(src->size(), (float)centrality->EtHFtowerSumPlus());
   vector<float> hiHFminus(src->size(), (float)centrality->EtHFtowerSumMinus());
   vector<float> hiHFeta4(src->size(), (float)centrality->EtHFtruncatedPlus()+centrality->EtHFtruncatedMinus());
   vector<float> hiHFplusEta4(src->size(), (float)centrality->EtHFtruncatedPlus());
   vector<float> hiHFminusEta4(src->size(), (float)centrality->EtHFtruncatedMinus());
   vector<float> hiHFhit(src->size(), (float)centrality->EtHFhitSum());
   vector<float> hiNpix(src->size(), (float)centrality->multiplicityPixel());
   vector<float> hiNpixelTracks(src->size(), (float)centrality->NpixelTracks());
   vector<float> hiNtracks(src->size(), (float)centrality->Ntracks());
   vector<float> hiEB(src->size(), (float)centrality->EtEBSum());
   vector<float> hiEE(src->size(), (float)centrality->EtEESum());
   vector<float> hiET(src->size(), (float)centrality->EtMidRapiditySum());

   writeValueMap(iEvent, src, hiHF,	"HFtowers");
   writeValueMap(iEvent, src, hiHFplus,	"HFtowersPlus");
   writeValueMap(iEvent, src, hiHFminus,	"HFtowersMinus");
   writeValueMap(iEvent, src, hiHFeta4,	"HFtowersTrunc");
   writeValueMap(iEvent, src, hiHFplusEta4,	"HFtowersPlusTrunc");
   writeValueMap(iEvent, src, hiHFminusEta4,	"HFtowersMinusTrunc");
   writeValueMap(iEvent, src, hiHFhit,	"HFhits");
   writeValueMap(iEvent, src, hiNpix,	"PixelHits");
   writeValueMap(iEvent, src, hiNpixelTracks,	"PixelTracks");
   writeValueMap(iEvent, src, hiNtracks,	"Tracks");
   writeValueMap(iEvent, src, hiEB,	"EB");
   writeValueMap(iEvent, src, hiEE,	"EE");
   writeValueMap(iEvent, src, hiET,	"ET");

}

void
HeavyIonCentralityInfo::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Candidate> > & handle,
        const std::vector<float> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;
    unique_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(std::move(valMap), label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyIonCentralityInfo);
