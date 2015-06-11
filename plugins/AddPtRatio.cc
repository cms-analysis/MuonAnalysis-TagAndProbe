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

#include "TVector3.h"


using namespace std;


//
// class declaration
//

class AddPtRatio : public edm::EDProducer{
	public:
		explicit AddPtRatio(const edm::ParameterSet&);
		~AddPtRatio();

	private:
		virtual void beginJob();
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob();

		// ----------member data ---------------------------
		const edm::InputTag probes_;    
		const edm::InputTag jets_;    
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
AddPtRatio::AddPtRatio(const edm::ParameterSet& iConfig):
	probes_(iConfig.getParameter<edm::InputTag>("probes")),
	jets_(iConfig.getParameter<edm::InputTag>("jets"))
{

	//now do what ever initialization is needed
	produces<edm::ValueMap<float> >("PtRatio");

}

AddPtRatio::~AddPtRatio()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//




// ------------ method called for each event  ------------
	void
AddPtRatio::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	Handle<View<reco::Candidate> > probes;
	iEvent.getByLabel(probes_, probes);

	//Aussi  mettre jet dans le fichier de configuration
	Handle<std::vector<reco::PFJet> > Jets;
	iEvent.getByLabel(jets_, Jets);
	//iEvent.getByLabel("ak4PFJetsCHS", Jets);

	//Output
	std::vector<double > ptratio;

	for( reco::CandidateView::const_iterator k = probes->begin(); k != probes->end(); ++ k){

		//Initialise loop variables
		double mupt_loop = 9999;
		double jetpt_loop = -9999;
		TVector3 pmu(0,0,0);
		TVector3 pjet(0,0,0);
		float dR = 9999;

		pmu.SetPtEtaPhi( k->pt(), k->eta(), k->phi());

		for( std::vector<reco::PFJet>::const_iterator l = Jets->begin(); l != Jets->end(); ++ l){

			pjet.SetPtEtaPhi( l->pt(), l->eta(), l->phi());

			if(dR > pmu.DrEtaPhi(pjet)){

				dR = pmu.DrEtaPhi(pjet); 
				mupt_loop = k->pt();
				jetpt_loop = l->pt();

			}

		}//end jet loop

		//
		//Fill the pt ratio
		//
		
		//No jets found
		if(dR == 9999){

			ptratio.push_back(1);

		}else{
		
		ptratio.push_back(mupt_loop/jetpt_loop);

		}


	}//end muon loop


	std::auto_ptr<ValueMap<float> > PtRatio(new ValueMap<float>());
	ValueMap<float>::Filler filler(*PtRatio);
	filler.insert(probes, ptratio.begin(), ptratio.end()); 
	filler.fill();
	iEvent.put(PtRatio,"PtRatio");

}


// ------------ method called once each job just before starting event loop  ------------
	void 
AddPtRatio::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
AddPtRatio::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddPtRatio);
