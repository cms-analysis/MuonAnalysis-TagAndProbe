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
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "../interface/PtRel.h"


using namespace std;


//
// class declaration
//

class AddPtRatioPtRel : public edm::EDProducer{
public:
  explicit AddPtRatioPtRel(const edm::ParameterSet&);
  ~AddPtRatioPtRel();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
  const edm::InputTag jets_;    
  const double dRmax_;
  const bool addLepToJet_;
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
AddPtRatioPtRel::AddPtRatioPtRel(const edm::ParameterSet& iConfig):
probes_(iConfig.getParameter<edm::InputTag>("probes")),
jets_(iConfig.getParameter<edm::InputTag>("jets")),
dRmax_(iConfig.getParameter<double>("dRmax")),
addLepToJet_(iConfig.getParameter<double>("addLepToJet"))
{

  //now do what ever initialization is needed
  produces<edm::ValueMap<float> >("PtRatio");
  produces<edm::ValueMap<float> >("PtRel");

}

AddPtRatioPtRel::~AddPtRatioPtRel()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//




// ------------ method called for each event  ------------
void
AddPtRatioPtRel::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<View<reco::Candidate> > probes;
  iEvent.getByLabel(probes_, probes);

  //Aussi  mettre jet dans le fichier de configuration
  Handle<std::vector<reco::PFJet> > Jets;
  iEvent.getByLabel(jets_, Jets);
  //iEvent.getByLabel("ak4PFJetsCHS", Jets);

  //Output
  std::vector<double > ptratio,ptrel;

  for( reco::CandidateView::const_iterator icand = probes->begin(); icand != probes->end(); ++ icand){

    //Initialise loop variables
    double mupt_loop = 9999;
    double jetpt_loop = -9999;
    float dR = 9999;
    double lepPx = 9999, lepPy = 9999, lepPz = 9999;
    double jetPt = 9999, jetEta = 9999, jetPhi = 9999, jetE = 9999;

    for( std::vector<reco::PFJet>::const_iterator ijet = Jets->begin(); ijet != Jets->end(); ++ ijet){

      
      if(dR > deltaR(*ijet, *icand)){

        dR = deltaR(*ijet, *icand);
        mupt_loop = icand->pt();
        jetpt_loop = ijet->pt();
        
        jetPt =  ijet->pt();
        jetEta = ijet->eta();
        jetPhi = ijet->phi();
        jetE =   ijet->energy();
        lepPx = icand->px();
        lepPy = icand->py();
        lepPz = icand->pz();
        
      }

    }//end jet loop

    //
    //Fill the pt ratio and pt rel
    //
    
    //No jets found
    if(dR > dRmax_){

      ptratio.push_back(1);
      ptrel.push_back(0);

    }else{
      
      ptratio.push_back(mupt_loop/jetpt_loop);
      ptrel.push_back(jetmet::getPtRel(jetPt, jetEta, jetPhi, jetE, lepPx, lepPy, lepPz, addLepToJet_));

    }


  }//end muon loop


  std::auto_ptr<ValueMap<float> > PtRatio(new ValueMap<float>());
  ValueMap<float>::Filler filler(*PtRatio);
  filler.insert(probes, ptratio.begin(), ptratio.end()); 
  filler.fill();
  iEvent.put(PtRatio,"PtRatio");

  std::auto_ptr<ValueMap<float> > PtRel(new ValueMap<float>());
  ValueMap<float>::Filler filler1(*PtRel);
  filler1.insert(probes, ptrel.begin(), ptrel.end()); 
  filler1.fill();
  iEvent.put(PtRel,"PtRel");

}


// ------------ method called once each job just before starting event loop  ------------
void 
AddPtRatioPtRel::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddPtRatioPtRel::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddPtRatioPtRel);
