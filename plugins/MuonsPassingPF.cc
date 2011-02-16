/** \class MuonsPassingPF
 *  Selector of muons identified by PF
 *
 *  $Date: 2011/01/30 12:09:24 $
 *  $Revision: 1.1 $
 *  \author G. Petrucciani (UCSD)
 */

// Base Class Headers
#include <cmath>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class MuonsPassingPF : public edm::EDProducer {

    public:
        /// constructor
        MuonsPassingPF(const edm::ParameterSet &iConfig) ;
        /// destructor
        ~MuonsPassingPF() ;
    
        /// method to be called at each event
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;

    private:
        /// Input collection of muons and of partice flow
        edm::InputTag muons_, pf_;

        StringCutObjectSelector<reco::PFCandidate> pfCut_;

        /// Perform matching by reference (works only if muons_ = "muons")
        bool matchByReference_;

}; // C++ note: you need a ';' at the end of the class declaration.


MuonsPassingPF::MuonsPassingPF(const edm::ParameterSet &iConfig) :
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    pf_(iConfig.getParameter<edm::InputTag>("pf")),
    pfCut_(iConfig.getParameter<std::string>("pfCut")),
    matchByReference_(iConfig.getParameter<bool>("matchByReference"))
{
    produces<edm::RefToBaseVector<reco::Muon> >();    
}

MuonsPassingPF::~MuonsPassingPF() 
{
}

void 
MuonsPassingPF::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    using namespace edm; 

    Handle<View<reco::Muon> > muons;
    iEvent.getByLabel(muons_, muons);

    Handle<View<reco::PFCandidate> > pf;
    iEvent.getByLabel(pf_, pf);
    
    /// prepare the vector for the output
    std::auto_ptr<RefToBaseVector<reco::Muon> > out(new RefToBaseVector<reco::Muon>());

    View<reco::PFCandidate>::const_iterator pfit, pfbegin = pf->begin(), pfend = pf->end();

    /// Now loop
    for (size_t i = 0, n = muons->size(); i < n; ++i) {
        // read the edm reference to the muon
        RefToBase<reco::Muon> muonRef = muons->refAt(i);

        /// Loop on PF
        for (pfit = pfbegin; pfit != pfend; ++pfit) {
            // Get muon ref, skip those that don't have one
            reco::MuonRef pfRef = pfit->muonRef(); 
            if (pfRef.isNull()) continue;
            if (!pfCut_(*pfit)) continue;

            // Perform the matching
            if (matchByReference_) {
                if (pfRef.id() != muonRef.id()) throw cms::Exception("Configuration") 
                    << "Cannot match by reference the muons from " << muons_.encode() << " (id " << muonRef.id() << ")" 
                    << " with the ones referenced by " << pf_.encode() << " (id " << pfRef.id() << ")\n";
                if (pfRef.key() == muonRef.key()) {
                    out->push_back(muonRef);
                    break;
                }
            } else if (std::abs(muonRef->eta() - pfRef->eta()) < 1e-5 && 
                       std::abs(muonRef->phi() - pfRef->phi()) < 1e-5 &&
                       std::abs(muonRef->pt()  - pfRef->pt() ) < 1e-5) {
                out->push_back(muonRef);
                break;
            }
        }
    }

    // Write the output to the event
    iEvent.put(out);
}

/// Register this as a CMSSW Plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonsPassingPF);
