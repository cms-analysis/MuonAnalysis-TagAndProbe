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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//
// class declaration
//

class MuonMVAIsoVariables : public edm::EDProducer {
public:
  explicit MuonMVAIsoVariables(const edm::ParameterSet&);
  ~MuonMVAIsoVariables();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
  const bool doOverlapRemoval_;
  const edm::InputTag goodElectrons_;    
  const edm::InputTag goodMuons_;    
  const edm::InputTag pfCandidates_;
  const edm::InputTag vertices_;
  const double dzCut_, photonPtMin_, neutralHadPtMin_;

  /// Store extra information in a ValueMap
  template<typename Hand, typename T>
  void storeMap(edm::Event &iEvent, 
                 const Hand & handle,
                 const std::vector<T> & values,
                 const std::string    & label) const ;
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
MuonMVAIsoVariables::MuonMVAIsoVariables(const edm::ParameterSet& iConfig):
    probes_(iConfig.getParameter<edm::InputTag>("probes")),
    doOverlapRemoval_(iConfig.getParameter<bool>("doOverlapRemoval")),
    goodElectrons_(doOverlapRemoval_ ? iConfig.getParameter<edm::InputTag>("goodElectrons") : edm::InputTag("NONE")),
    goodMuons_(doOverlapRemoval_ ? iConfig.getParameter<edm::InputTag>("goodMuons") : edm::InputTag("NONE")),
    pfCandidates_(iConfig.getParameter<edm::InputTag>("pfCandidates")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    dzCut_(iConfig.getParameter<double>("dzCut")),
    photonPtMin_(iConfig.getParameter<double>("photonPtMin")),
    neutralHadPtMin_(iConfig.getParameter<double>("neutralHadPtMin"))
{
  produces<edm::ValueMap<float> >("ChargedIsoDR0p0To0p1");
  produces<edm::ValueMap<float> >("ChargedIsoDR0p1To0p2");
  produces<edm::ValueMap<float> >("ChargedIsoDR0p2To0p3");
  produces<edm::ValueMap<float> >("ChargedIsoDR0p3To0p4");
  produces<edm::ValueMap<float> >("ChargedIsoDR0p4To0p5");
  produces<edm::ValueMap<float> >("GammaIsoDR0p0To0p1");
  produces<edm::ValueMap<float> >("GammaIsoDR0p1To0p2");
  produces<edm::ValueMap<float> >("GammaIsoDR0p2To0p3");
  produces<edm::ValueMap<float> >("GammaIsoDR0p3To0p4");
  produces<edm::ValueMap<float> >("GammaIsoDR0p4To0p5");
  produces<edm::ValueMap<float> >("NeutralHadronIsoDR0p0To0p1");
  produces<edm::ValueMap<float> >("NeutralHadronIsoDR0p1To0p2");
  produces<edm::ValueMap<float> >("NeutralHadronIsoDR0p2To0p3");
  produces<edm::ValueMap<float> >("NeutralHadronIsoDR0p3To0p4");
  produces<edm::ValueMap<float> >("NeutralHadronIsoDR0p4To0p5");
  produces<edm::ValueMap<float> >("PFCharged");
  produces<edm::ValueMap<float> >("PFNeutral");
  produces<edm::ValueMap<float> >("PFPhotons");
  produces<edm::ValueMap<float> >("SumDeltaR");
  produces<edm::ValueMap<float> >("DeltaRMean");
  produces<edm::ValueMap<float> >("Density");
}


MuonMVAIsoVariables::~MuonMVAIsoVariables()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonMVAIsoVariables::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByLabel(probes_, probes);

  Handle<View<reco::GsfElectron> > goodElectrons;
  Handle<View<reco::Muon> >        goodMuons;
  if (doOverlapRemoval_) {
      iEvent.getByLabel(goodMuons_,     goodMuons);
      iEvent.getByLabel(goodElectrons_, goodElectrons);
  }

  Handle<View<reco::PFCandidate> > pfCandidates;
  Handle<reco::VertexCollection  > vertices;
  iEvent.getByLabel(pfCandidates_, pfCandidates);
  iEvent.getByLabel(vertices_, vertices);
  const reco::Vertex &vtx = vertices->at(0);

  View<reco::Muon>::const_iterator probe, endprobes=probes->end();
  View<reco::PFCandidate>::const_iterator iP, beginpf = pfCandidates->begin(), endpf=pfCandidates->end();
  unsigned int n = probes->size();

  std::vector<float> chargedIso[5], gammaIso[5], neutralHadronIso[5];
  for (int i = 0; i < 5; ++i) {
    chargedIso[i] = std::vector<float>(n);
    gammaIso[i] = std::vector<float>(n);
    neutralHadronIso[i] = std::vector<float>(n);
  }
  std::vector<float> pfCharged(n), pfNeutral(n), pfPhotons(n), sumDeltaR(n), deltaRMean(n), density(n);  

  // loop on PROBES
  unsigned int imu = 0;
  for (probe = probes->begin(); probe != endprobes; ++probe, ++imu) {
    const reco::Muon &mu = *probe;

    double z0 = 0; 
    if (mu.track().isNonnull()) z0 = mu.track()->dz(vtx.position()); 

    int nPF = 0;
    for (iP = beginpf; iP != endpf; ++iP) {
        double dr = deltaR(*iP, mu);
        int idr = floor(dr/0.1);
        if (idr >= 5) continue; 

        bool hastk = iP->trackRef().isNonnull();
        if (hastk && mu.track() == iP->trackRef()) continue;
        
        if (doOverlapRemoval_) {
            bool veto = false;
            for (View<reco::GsfElectron>::const_iterator iE = goodElectrons->begin(), elend = goodElectrons->end(); iE != elend; ++iE) {
                if (iP->gsfTrackRef().isNonnull() && iP->gsfTrackRef() == iE->gsfTrack()) { veto = true; break; }
                if (hastk && iP->trackRef() == iE->closestCtfTrackRef()) { veto = true; break; }
                if (fabs(iE->superCluster()->eta()) >= 1.479) {
                    double drE = deltaR(*iE, *iP);
                    if (hastk && drE < 0.015) { veto = true; break; }
                    if (iP->particleId() == reco::PFCandidate::gamma && drE < 0.08) { veto = true; break; }
                }
            }
            if (veto) continue;
            for (View<reco::Muon>::const_iterator iM = goodMuons->begin(), muend = goodMuons->end(); iM != muend; ++iM) {
                if (hastk && iP->trackRef() == iM->innerTrack()) { veto = true; break; }
                double drM = deltaR(*iM, *iP);
                if (drM < 0.01)  { veto = true; break; }
            }
            if (veto) continue;
        }

        if (hastk) {
            if (fabs(iP->trackRef()->dz(vtx.position()) - z0) >= dzCut_) continue;
            if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
            if (fabs(mu.eta()) && dr < 0.01) continue;
            chargedIso[idr][imu] += iP->pt();
        } else if (iP->particleId() == reco::PFCandidate::gamma) {
            if (iP->pt() <= photonPtMin_) continue;
            gammaIso[idr][imu] += iP->pt();
        } else {
            if (iP->pt() <= neutralHadPtMin_) continue;
            neutralHadronIso[idr][imu] += iP->pt();
        }
        nPF++;
        sumDeltaR[imu] += dr;
        density[imu] += iP->pt() / dr; 
    }
    for (int idr = 0; idr < 5; ++idr) {
        pfCharged[imu] += chargedIso[idr][imu];
        pfNeutral[imu] += gammaIso[idr][imu];
        pfPhotons[imu] += neutralHadronIso[idr][imu];
    }
    deltaRMean[imu] = sumDeltaR[imu] / std::max<int>(1, nPF);
  }// end loop on probes

  storeMap(iEvent, probes, pfCharged, "PFCharged");
  storeMap(iEvent, probes, pfPhotons, "PFPhotons");
  storeMap(iEvent, probes, pfNeutral, "PFNeutral");
  storeMap(iEvent, probes, sumDeltaR, "SumDeltaR");
  storeMap(iEvent, probes, deltaRMean, "DeltaRMean");
  storeMap(iEvent, probes, density, "Density");
  for (int idr = 0; idr < 5; ++idr) {
      char buff[50]; sprintf(buff, "DR0p%dTo0p%d", idr, idr+1);
      storeMap(iEvent, probes, chargedIso[idr], std::string("ChargedIso")+buff);
      storeMap(iEvent, probes, gammaIso[idr], std::string("GammaIso")+buff);
      storeMap(iEvent, probes, neutralHadronIso[idr], std::string("NeutralHadronIso")+buff);
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMVAIsoVariables::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMVAIsoVariables::endJob() {
}

template<typename Hand, typename T>
void
MuonMVAIsoVariables::storeMap(edm::Event &iEvent,
                     const Hand & handle,
                     const std::vector<T> & values,
                     const std::string    & label) const {
    using namespace edm; using namespace std;
    auto_ptr<ValueMap<T> > valMap(new ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMVAIsoVariables);
