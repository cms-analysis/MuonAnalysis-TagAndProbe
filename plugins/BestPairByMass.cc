// -*- C++ -*-
//
// Package:    BestPairByMass
// Class:      BestPairByMass
// 
/**\class BestPairByMass BestPairByMass.cc MuonAnalysis/TagAndProbe/src/BestPairByMass.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: BestPairByMass.cc,v 1.1 2010/05/27 08:30:39 gpetrucc Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

// class decleration
class BestPairByMass : public edm::EDProducer {
    public:
        explicit BestPairByMass(const edm::ParameterSet&);
        ~BestPairByMass();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        /// The Tags
        edm::InputTag pairs_;

        /// Mass value and window
        double mass_;

        

};

BestPairByMass::BestPairByMass(const edm::ParameterSet &iConfig) :
    pairs_(iConfig.getParameter<edm::InputTag>("pairs")),
    mass_(iConfig.getParameter<double>("mass"))
{
    produces<edm::RefToBaseVector<reco::Candidate> >();
}

BestPairByMass::~BestPairByMass() 
{
}

void
BestPairByMass::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::Candidate> > pairs; 
    iEvent.getByLabel(pairs_, pairs);

    std::auto_ptr<edm::RefToBaseVector<reco::Candidate> > out(new edm::RefToBaseVector<reco::Candidate>());

    // Filter and fill
    size_t i, n = pairs->size(); 
    std::vector<int8_t> marked(n, 0); // 0 = unprocessed; 1 = processed, good; -1 = processed, bad
    edm::View<reco::Candidate>::const_iterator it, ed = pairs->end();
    for (it = pairs->begin(), i = 0; it != ed; ++it, ++i) {
        if (marked[i] != 0) continue;
        const reco::Candidate & pair = *it;
        if (pair.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "Must have 2 daughters, while pair " << i << " has " << pair.numberOfDaughters() << ".\n";
        reco::CandidateBaseRef tag = pair.daughter(0)->masterClone();
        if (tag.isNull()) throw cms::Exception("CorruptData") << "Null tag reference for pair " << i << ".\n";

        marked[i] = 1;
        size_t found = i; 
        double dm = fabs(pair.mass() - mass_);
        for  (size_t j = i+1; j < n; ++j) {
            if (marked[j] != 0) continue;
            const reco::Candidate & pair2 = (*pairs)[j];
            if (pair2.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "Must have 2 daughters, while pair " << j << " has " << pair2.numberOfDaughters() << ".\n";
            reco::CandidateBaseRef tag2 = pair2.daughter(0)->masterClone();
            if (tag2.isNull()) throw cms::Exception("CorruptData") << "Null tag reference for pair " << j << ".\n";
            if (tag2 == tag) {
                double dm2 = fabs(pair2.mass() - mass_);
                if (dm2 < dm) {
                    marked[found] = -1;
                    marked[j]     = 1;
                    found = j;
                    dm = dm2;
                } else {
                    marked[j] = -1;
                }
            }
        }

        out->push_back(pairs->refAt(i));
    }

    iEvent.put(out);
}

DEFINE_FWK_MODULE(BestPairByMass);
