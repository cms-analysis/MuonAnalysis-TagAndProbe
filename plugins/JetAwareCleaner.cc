#ifndef JetAwareCleaner_h
#define JetAwareCleaner_h

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

#include <map>
#include <vector>

class JetAwareCleaner : public edm::global::EDProducer<> {

public:

  explicit JetAwareCleaner(const edm::ParameterSet& iConfig);
  ~JetAwareCleaner();

  virtual void produce(edm::StreamID streamID, edm::Event & iEvent,
  const edm::EventSetup & iSetup) const;

private:

  edm::EDGetTokenT<std::vector<reco::PFJet> > jetCollectionTag_;
  edm::EDGetTokenT<reco::CandidateView> leptonCollectionTag_;
  
  edm::EDGetTokenT<reco::JetCorrector> tagL1Corrector_;
  edm::EDGetTokenT<reco::JetCorrector> tagL1L2L3ResCorrector_;

  double dRmax_, dRCandProbeVeto_;
};

#endif



JetAwareCleaner::JetAwareCleaner(const edm::ParameterSet& iConfig) {

  edm::InputTag jetcollection = iConfig.getParameter<edm::InputTag>("RawJetCollection");
  jetCollectionTag_ = consumes<reco::PFJetCollection>(jetcollection);

  edm::InputTag leptoncollection = iConfig.getParameter<edm::InputTag>("LeptonCollection");
  leptonCollectionTag_ = consumes<reco::CandidateView>(leptoncollection);

  edm::InputTag l1Cortag = iConfig.getParameter<edm::InputTag>("L1Corrector");
  tagL1Corrector_ = consumes<reco::JetCorrector>(l1Cortag);
  edm::InputTag l1l2l3ResCortag = iConfig.getParameter<edm::InputTag>("L1L2L3ResCorrector");
  tagL1L2L3ResCorrector_ = consumes<reco::JetCorrector>(l1l2l3ResCortag);

  dRmax_  = iConfig.getParameter<double>("dRmax");
  dRCandProbeVeto_  = iConfig.getParameter<double>("dRCandProbeVeto");

  // produces vector of pfJet
  produces<std::vector<reco::PFJet> >();
}


JetAwareCleaner::~JetAwareCleaner() {

}


void JetAwareCleaner::produce(edm::StreamID streamID, edm::Event & iEvent,
const edm::EventSetup & iSetup) const {

  edm::Handle<reco::PFJetCollection> jets;       
  iEvent.getByToken (jetCollectionTag_, jets);    

  edm::Handle<reco::CandidateView> leptons;               
  iEvent.getByToken (leptonCollectionTag_, leptons);       

  edm::Handle<reco::JetCorrector> correctorL1L2L3Res;
  iEvent.getByToken(tagL1L2L3ResCorrector_, correctorL1L2L3Res);

  edm::Handle<reco::JetCorrector> correctorL1;
  iEvent.getByToken(tagL1Corrector_, correctorL1);

  std::map<unsigned int, std::pair<int,unsigned int> > match;
  std::map<unsigned int, std::pair<int,unsigned int> >::const_iterator it;

  //super ugly....
  for(unsigned int il=0;il<leptons->size();il++) {
    match[il]=std::make_pair(1000,-1);

    for(unsigned int ij=0;ij<jets->size();ij++) {
      reco::PFJetRef jet(jets,ij);

      float dR2=deltaR2( (*leptons)[il].p4(), jet->p4());
      if(dR2<(dRmax_*dRmax_)) {
        if(dR2<match[il].first) {
          match[il]=std::make_pair(dR2,ij);
        }
      }
    }
  }


  std::vector<reco::PFJet> *jetCol = new std::vector<reco::PFJet>();
  
  for(unsigned int ij=0;ij<jets->size();ij++) {
    reco::PFJetRef jet(jets,ij);
    
    double jecL1L2L3Res = correctorL1L2L3Res->correction(*jet);
    double jecL1 = correctorL1->correction(*jet);
    
    reco::Candidate::LorentzVector p4J=jet->p4();
    bool mlep=false;
    for(it=match.begin();it!=match.end();++it) {
      if(it->second.second == ij ) {
        mlep=true;
        
        if( (p4J-(*leptons)[it->first].p4()).Rho()<dRCandProbeVeto_ ) { p4J=(*leptons)[it->first].p4(); break;} //protection against jets built only with a lepton
        p4J -= (*leptons)[it->first].p4()/jecL1;
        p4J *= jecL1L2L3Res;
        p4J += (*leptons)[it->first].p4();
        break;
      }
    }//matched lepton loop

    if(!mlep) p4J *= jecL1L2L3Res;

    //build the cleaned jet
    reco::PFJet cleanJet(p4J, jet->vertex(), jet->getSpecific());
    jetCol->push_back(cleanJet);
  }//jet loop
  
  std::auto_ptr<std::vector<reco::PFJet> > recoCleanedJets(jetCol);
  iEvent.put(recoCleanedJets);

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(JetAwareCleaner);
