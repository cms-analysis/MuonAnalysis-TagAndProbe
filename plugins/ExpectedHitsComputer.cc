#include "FWCore/Framework/interface/EDProducer.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include <DataFormats/RecoCandidate/interface/RecoCandidate.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h>
#include <RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h>
#include <Geometry/Records/interface/GlobalTrackingGeometryRecord.h>


#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include <TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h>
#include <TrackingTools/DetLayers/interface/GeometricSearchDet.h> 
#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>
#include <TrackingTools/MeasurementDet/interface/MeasurementDet.h>
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include <TrackingTools/TransientTrack/interface/TrackTransientTrack.h>
#include <TrackingTools/TransientTrack/interface/GsfTransientTrack.h>

#include "DataFormats/DetId/interface/DetId.h"
#include <DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h>

#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
#include <DataFormats/SiPixelDetId/interface/PXBDetId.h>

class ExpectedHitsComputer : public edm::EDProducer {
public:
  explicit ExpectedHitsComputer(const edm::ParameterSet & iConfig);
  virtual ~ExpectedHitsComputer() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

private:
  edm::InputTag input_;
  bool useGsfTrack_;
  StringCutObjectSelector<reco::Candidate,true> objCut_; 
  std::string thePropName;      
  std::string theNavSchoolName; 
  std::string theMeasTkName;    
};

ExpectedHitsComputer::ExpectedHitsComputer(const edm::ParameterSet & iConfig) :
  input_(iConfig.getParameter<edm::InputTag>("inputColl")),
  useGsfTrack_(iConfig.getParameter<bool>("useGsfTrack")),
  objCut_(iConfig.existsAs<std::string>("objectSelection") ? iConfig.getParameter<std::string>("objectSelection") : "", true)
{
  produces<edm::ValueMap<int> >("in").setBranchAlias("in");
  produces<edm::ValueMap<int> >("out").setBranchAlias("out");

  thePropName      = iConfig.getParameter<std::string>("propagator");  
  theNavSchoolName = iConfig.getParameter<std::string>("navigationSchool");
  theMeasTkName    = iConfig.getParameter<std::string>("measurementTracker");
}


ExpectedHitsComputer::~ExpectedHitsComputer()
{
}

void 
ExpectedHitsComputer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;  
  using namespace std;

  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<Propagator> theProp;
  edm::ESHandle<NavigationSchool> theNavSchool;
  edm::ESHandle<MeasurementTracker> theMeasTk;
  edm::ESHandle<GlobalTrackingGeometry> theGeo;

  iSetup.get<IdealMagneticFieldRecord>().get(theMF);
  iSetup.get<TrackingComponentsRecord>().get(thePropName,theProp);
  iSetup.get<NavigationSchoolRecord>().get(theNavSchoolName, theNavSchool); 
  NavigationSetter setter( *theNavSchool );

  iSetup.get<CkfComponentsRecord>().get(theMeasTkName,theMeasTk); 
  iSetup.get<GlobalTrackingGeometryRecord>().get(theGeo); 

  Chi2MeasurementEstimator estimator(30.,-3.0);
  //Chi2MeasurementEstimator estimator(30.,3.0);
    


  // read input
  Handle<View<reco::RecoCandidate> > inputCands;
  iEvent.getByLabel(input_,  inputCands);

  // prepare vector for output    
  std::vector<int> innerValues;
  std::vector<int> outerValues;
    
  // fill
  View<reco::RecoCandidate>::const_iterator cand, endcands = inputCands->end();
  for (cand = inputCands->begin(); cand != endcands; ++cand) {

    TrajectoryStateOnSurface tsosInner;
    TrajectoryStateOnSurface tsosOuter;
    DetId idInner;
    DetId idOuter;

    if(useGsfTrack_){   
      if(cand->gsfTrack().isNull()) {
	//cout << "ERROR: null track for ELE" << endl;
	innerValues.push_back(0);
	outerValues.push_back(0);
	continue;
      }
      //cout << "is ele" << endl;
      reco::GsfTransientTrack tt(*cand->gsfTrack(),theMF.product(),theGeo);
      tsosInner = tt.innermostMeasurementState();
      tsosOuter = tt.outermostMeasurementState();
      idInner = DetId(tt.innerDetId());
      idOuter = DetId(tt.outerDetId());
    }else{
      if(cand->track().isNull()) {
	///cout << "ERROR: null track for MUON" << endl;
	innerValues.push_back(0);
	outerValues.push_back(0);
	continue;
      }
      //cout << "is mu" << endl;
      reco::TrackTransientTrack tt(*cand->track(),theMF.product(),theGeo);
      tsosInner = tt.innermostMeasurementState();
      tsosOuter = tt.outermostMeasurementState();
      idInner = DetId(tt.innerDetId());
      idOuter = DetId(tt.outerDetId());
    }
    
    //cout << "tsos pos.perp: " << tsos.globalPosition().perp() << endl;
    //cout << "DetId subdetId: " << id.subdetId() << endl;

    /*
    if(id.subdetId() == 1){
      PXBDetId tmpId(id);
      cout << "is on pxb layer: " << tmpId.layer() << endl;
    }
    */


    const DetLayer* innerLayer = theMeasTk->geometricSearchTracker()->idToLayer(idInner);
    const DetLayer* outerLayer = theMeasTk->geometricSearchTracker()->idToLayer(idOuter);
    //cout << "innerLayer radius: " << innerLayer->position().perp() << endl;
    

    PropagationDirection dirForInnerLayers = oppositeToMomentum;
    PropagationDirection dirForOuterLayers = alongMomentum;



    std::vector< const DetLayer * > innerCompLayers = 
      innerLayer->compatibleLayers(*tsosInner.freeState(),dirForInnerLayers);

    std::vector< const DetLayer * > outerCompLayers = 
      outerLayer->compatibleLayers(*tsosOuter.freeState(),dirForOuterLayers);

    //cout << "innerCompLayers size: " << innerCompLayers.size() << endl; 

    int counter(0);
    for(vector<const DetLayer *>::const_iterator it=innerCompLayers.begin(); it!=innerCompLayers.end();
	++it){
      vector< GeometricSearchDet::DetWithState > detWithState = (*it)->compatibleDets(tsosInner,
										      *theProp.product(),
										      estimator);
      if(!detWithState.size()) continue;
      DetId id = detWithState.front().first->geographicalId();
      const MeasurementDet* measDet = theMeasTk->idToDet(id);	
      if(measDet->isActive()){	
      //if(1){
	counter++;
	//InvalidTrackingRecHit  tmpHit(id,TrackingRecHit::missing);
	////track.setTrackerExpectedHitsInner(tmpHit,counter); 	 
	//cout << "WARNING: this hit is marked as lost because the detector was marked as active" << endl;
      }else{
	//cout << "WARNING: this hit is NOT marked as lost because the detector was marked as inactive" << endl;
      }
    }//end loop over layers
    innerValues.push_back(counter);

    counter=0;
    for(vector<const DetLayer *>::const_iterator it=outerCompLayers.begin(); it!=outerCompLayers.end();
	++it){
      vector< GeometricSearchDet::DetWithState > detWithState = (*it)->compatibleDets(tsosOuter,
										      *theProp.product(),
										      estimator);
      if(!detWithState.size()) continue;
      DetId id = detWithState.front().first->geographicalId();
      const MeasurementDet* measDet = theMeasTk->idToDet(id);	
      if(measDet->isActive()){	
      //if(1){
	counter++;
	//InvalidTrackingRecHit  tmpHit(id,TrackingRecHit::missing);
	////track.setTrackerExpectedHitsInner(tmpHit,counter); 	 
	//cout << "WARNING: this hit is marked as lost because the detector was marked as active" << endl;
      }else{
	//cout << "WARNING: this hit is NOT marked as lost because the detector was marked as inactive" << endl;
      }
    }//end loop over layers
    outerValues.push_back(counter);

    //cout << "counter: " << counter << endl;
  }
    
  // convert into ValueMap and store
  std::auto_ptr<ValueMap<int> > valMapInner(new ValueMap<int>());
  std::auto_ptr<ValueMap<int> > valMapOuter(new ValueMap<int>());

  ValueMap<int>::Filler fillerInner(*valMapInner);
  ValueMap<int>::Filler fillerOuter(*valMapOuter);

  fillerInner.insert(inputCands, innerValues.begin(), innerValues.end());
  fillerOuter.insert(inputCands, outerValues.begin(), outerValues.end());

  fillerInner.fill();
  fillerOuter.fill();

  iEvent.put(valMapInner,"in");
  iEvent.put(valMapOuter,"out");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ExpectedHitsComputer);
