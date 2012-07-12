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

#include "DataFormats/Common/interface/ValueMap.h"
#include <DataFormats/MuonReco/interface/Muon.h>

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


//
// class declaration
//

class MuonDxyPVdzmin : public edm::EDProducer {
public:
  explicit MuonDxyPVdzmin(const edm::ParameterSet&);
  ~MuonDxyPVdzmin();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
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
MuonDxyPVdzmin::MuonDxyPVdzmin(const edm::ParameterSet& iConfig):
probes_(iConfig.getParameter<edm::InputTag>("probes"))

{
  produces<edm::ValueMap<float> >("dxyBS");
  produces<edm::ValueMap<float> >("dxyPVdzmin");
  produces<edm::ValueMap<float> >("dzPV");
}


MuonDxyPVdzmin::~MuonDxyPVdzmin()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonDxyPVdzmin::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByLabel(probes_,  probes);
  
  edm::Handle<std::vector<reco::Vertex> > primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primaryVerticesHandle);
      
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot" ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;
  math::XYZPoint BSPosition;
  BSPosition = bs.position();

  // prepare vector for output    
  std::vector<double> muon_dxyBS;
  std::vector<double> muon_dxy;
  std::vector<double> muon_dz;

  // fill
  View<reco::Muon>::const_iterator probe, endprobes = probes->end();

  // loop on PROBES
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    Double_t dxyBS = -99999.;
    Double_t dxy_ivtx_dzmin=65535;
    Double_t dzPV = -99999.;

    if(probe->innerTrack().isNonnull()){
      
      edm::ESHandle<TransientTrackBuilder> ttBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
      
      reco::TrackRef muonTrack = probe->innerTrack();
      reco::TransientTrack tTrack = ttBuilder->build(*muonTrack);
      
      Double_t dzmin = 1e6;

      dxyBS = muonTrack->dxy(BSPosition);
      muon_dxyBS.push_back(dxyBS);
      dzPV = muonTrack->dz(primaryVerticesHandle->at(0).position());
      muon_dz.push_back(dzPV);

      for(unsigned int iVtx=0; iVtx<primaryVerticesHandle->size(); iVtx++){

        Double_t pvx,pvy,pvz,bsx,bsy,bsz;
        pvx=primaryVerticesHandle->at(iVtx).x();
        pvy=primaryVerticesHandle->at(iVtx).y();
        pvz=primaryVerticesHandle->at(iVtx).z();
        bsx = BSPosition.x();
        bsy = BSPosition.y();
        bsz = BSPosition.z();
        Double_t dz_pv_bs = std::abs(bsz-pvz);
        Double_t rxy_pv_bs = std::sqrt( (bsx-pvx)*(bsx-pvx) + (bsy-pvy)*(bsy-pvy) );
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
  std::auto_ptr<ValueMap<float> > dxyBS(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > dxyPVdzmin(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > dzPV(new ValueMap<float>());
  ValueMap<float>::Filler filler0(*dxyBS);
  ValueMap<float>::Filler filler1(*dxyPVdzmin);
  ValueMap<float>::Filler filler2(*dzPV);
  filler0.insert(probes, muon_dxyBS.begin(), muon_dxyBS.end());
  filler1.insert(probes, muon_dxy.begin(), muon_dxy.end());
  filler2.insert(probes, muon_dz.begin(), muon_dz.end());
  filler0.fill();
  filler1.fill();
  filler2.fill();

  iEvent.put(dxyBS, "dxyBS");
  iEvent.put(dxyPVdzmin, "dxyPVdzmin");
  iEvent.put(dzPV, "dzPV");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonDxyPVdzmin::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonDxyPVdzmin::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonDxyPVdzmin);
