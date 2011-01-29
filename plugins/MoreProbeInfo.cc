// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include <DataFormats/MuonReco/interface/Muon.h>

#include <FWCore/Common/interface/TriggerNames.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>
#include "DataFormats/Math/interface/deltaPhi.h"

#include <iostream>
#include <map>

#include "TH1D.h"


// MARIO NTUPLIZER
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
// #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"




//
// class declaration
//

class MoreProbeInfo : public edm::EDProducer {
public:
  explicit MoreProbeInfo(const edm::ParameterSet&);
  ~MoreProbeInfo();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  float Phi( float&, float&);

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
  const edm::InputTag DT4Segments_;            

  bool makeControlHisto_;
  int MBx_;

  const edm::InputTag trigRes_;

  std::map<std::string,TH1D*> h1_;

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
MoreProbeInfo::MoreProbeInfo(const edm::ParameterSet& iConfig):
probes_(iConfig.getParameter<edm::InputTag>("probes"))

{

  produces<edm::ValueMap<float> >("dxyPVdzmin");
  //register your products
  /* Examples
    produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
  */
  //now do what ever other initialization is needed

}


MoreProbeInfo::~MoreProbeInfo()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MoreProbeInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<pat::Muon> > probes;
  iEvent.getByLabel(probes_,  probes);
  
  edm::Handle<std::vector<reco::Vertex> > primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primaryVerticesHandle);
      
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot" ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;
  math::XYZPoint BSPosition;
  BSPosition = bs.position();

  // prepare vector for output    
  std::vector<double> muon_dxy;

  // fill
  View<pat::Muon>::const_iterator probe, endprobes = probes->end();

  // loop on PROBES
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    Double_t dxy_ivtx_dzmin=65535;
    
    if(probe->innerTrack().isNonnull()){
      
      edm::ESHandle<TransientTrackBuilder> ttBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
      
      reco::TrackRef muonTrack = probe->innerTrack();
      reco::TransientTrack tTrack = ttBuilder->build(*muonTrack);
      
      Double_t dzmin = 1e6;
      Int_t ivtx_dzmin = -1;
      
      for(unsigned int iVtx=0; iVtx<primaryVerticesHandle->size(); iVtx++){
        Double_t pvx,pvy,pvz,bsx,bsy,bsz;
        pvx=primaryVerticesHandle->at(iVtx).x();
        pvy=primaryVerticesHandle->at(iVtx).y();
        pvz=primaryVerticesHandle->at(iVtx).z();
        bsx = BSPosition.x();
        bsy = BSPosition.y();
        bsz = BSPosition.z();
        Double_t dz_pv_bs = TMath::Abs(bsz-pvz);
        Double_t rxy_pv_bs = TMath::Sqrt( (bsx-pvx)*(bsx-pvx) + (bsy-pvy)*(bsy-pvy) );
        std::pair<bool,Measurement1D> dxy_fromIP = IPTools::absoluteTransverseImpactParameter(tTrack,primaryVerticesHandle->at(iVtx));
        std::pair<bool,Measurement1D> dxyz_fromIP = IPTools::absoluteImpactParameter3D(tTrack,primaryVerticesHandle->at(iVtx));
        Double_t dzmin_loop = sqrt(dxyz_fromIP.second.value()*dxyz_fromIP.second.value() - dxy_fromIP.second.value()*dxy_fromIP.second.value());
        
        if(!primaryVerticesHandle->at(iVtx).isFake()
           && primaryVerticesHandle->at(iVtx).ndof()>4
           && dz_pv_bs <=24
           && rxy_pv_bs <=2
           && dzmin_loop < dzmin
          ){
          dzmin = dzmin_loop;
          ivtx_dzmin = iVtx;
          dxy_ivtx_dzmin = dxy_fromIP.second.value();
        }
      }
      if(dzmin<1 && dxy_ivtx_dzmin<100){
        muon_dxy.push_back(dxy_ivtx_dzmin);
      }else{
        muon_dxy.push_back(65535);
      }
    }else{
      muon_dxy.push_back(dxy_ivtx_dzmin);
    }

  }// end loop on probes

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > dxyPVdzmin(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*dxyPVdzmin);
  filler1.insert(probes, muon_dxy.begin(), muon_dxy.end());
  filler1.fill();
  iEvent.put(dxyPVdzmin, "dxyPVdzmin");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MoreProbeInfo::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MoreProbeInfo::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MoreProbeInfo);
