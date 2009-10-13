/** \class MuonSelectorJPsi
 *  Example selector of muons.
 *
 *  $Date: 2009/10/12 06:55:33 $
 *  $Revision: 1.1 $
 *  \author G. Petrucciani (SNS Pisa)
 */

// Base Class Headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

class MuonSelectorJPsi : public edm::EDProducer {

    public:
        /// constructor
        MuonSelectorJPsi(const edm::ParameterSet &iConfig) ;
        /// destructor
        ~MuonSelectorJPsi() ;
    
        /// method to be called at each event
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;

    private:
        //--- put here your data members ---

        /// Input collection of muons
        edm::InputTag src_;
        
        bool selectGlobalMuons_;

}; // C++ note: you need a ';' at the end of the class declaration.


/// Constructor: read the configuration, initialize data members, declare what to produce
MuonSelectorJPsi::MuonSelectorJPsi(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    selectGlobalMuons_(iConfig.getParameter<bool>("selectGlobalMuons"))
{
    // declare what we produce: a vector of references to Muon-like objects (can be reco::Muon or pat::Muons)
    // subsequent modules should read this with View<reco::Muon>
    produces<edm::RefToBaseVector<reco::Muon> >();    
}

/// Destructor: here there is usually nothing to do
MuonSelectorJPsi::~MuonSelectorJPsi() 
{
}

/// Produce: the method where we do something
void 
MuonSelectorJPsi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    using namespace edm; // avoid the prefix 'edm::' in front of Handle, View, Ref..., 

    // read the input from the event.
    // View<reco::Muon> can read any collection of reco::Muon or pat::Muon, or any collection of references to them
    Handle<View<reco::Muon> > src;
    iEvent.getByLabel(src_, src);
    
    // if you need to read other stuff, e.g. other collections, associations or condition data, do it here.

    /// prepare the vector for the output
    std::auto_ptr<RefToBaseVector<reco::Muon> > out(new RefToBaseVector<reco::Muon>());

    /// Now loop
    for (size_t i = 0, n = src->size(); i < n; ++i) {
        // read the edm reference to the muon
        RefToBase<reco::Muon> muonRef = src->refAt(i);

        // get also the muon, to make the code more readable
        const reco::Muon & mu = *muonRef;

        /// now we perform some example selection...
        if (selectGlobalMuons_) {
            if (mu.isGlobalMuon() && mu.globalTrack()->chi2()/mu.globalTrack()->ndof()< 20.0) {
                out->push_back(muonRef);
            }
        } else {
            if (!mu.isGlobalMuon() && mu.isTrackerMuon() && mu.track()->numberOfValidHits() > 12 && 
                    (muon::isGoodMuon(mu, muon::TM2DCompatibilityTight) || muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtLoose)) &&  
                    mu.track()->chi2()/mu.track()->ndof()< 5.0) {
                out->push_back(muonRef);
            }
        }
        
    }

    // Write the output to the event
    iEvent.put(out);
}

/// Register this as a CMSSW Plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonSelectorJPsi);
