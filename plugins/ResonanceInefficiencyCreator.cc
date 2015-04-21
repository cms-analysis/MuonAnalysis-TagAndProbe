// -*- C++ -*-
//
// Package:    ResonanceInefficiencyCreator
// Class:      ResonanceInefficiencyCreator
// 
/**\class ResonanceInefficiencyCreator ResonanceInefficiencyCreator.cc MuonAnalysis/TagAndProbe/src/ResonanceInefficiencyCreator.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nov 16 16:12 (lxplus231.cern.ch)
//         Created:  Sun Nov 16 16:14:09 CET 2008
// $Id: ResonanceInefficiencyCreator.cc,v 1.2 2009/09/30 14:03:35 gpetrucc Exp $
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

template<typename T> 
struct ValueVectorTrait { 
    typedef std::vector<T> type; 
};
template<> struct ValueVectorTrait<reco::Candidate> { 
    typedef edm::OwnVector<reco::Candidate> type; 
};
//
// class decleration
template<typename T>
class ResonanceInefficiencyCreator : public edm::EDProducer {
    public:
        explicit ResonanceInefficiencyCreator(const edm::ParameterSet&);
        ~ResonanceInefficiencyCreator();

        typedef typename ValueVectorTrait<T>::type PlainVecT;
        typedef edm::Ref<PlainVecT>       RefT;
        typedef edm::RefVector<PlainVecT> RefVecT;
        typedef edm::RefToBaseVector<T>   RefBaseVecT;

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        /// The RECO objects
        edm::InputTag src_;

        /// The Tags
        edm::InputTag tags_;

        /// Mass value and window
        double mass_, massMin_, massMax_;

        /// Drop only best match or all those in window
        bool onlyBestMatch_;

        typedef StringObjectFunction<T> Function;
        /// The Probability
        Function probability_;

        /// The Random number generator
  //std::auto_ptr<CLHEP::RandFlat> flatDistribution_;

        enum OutputMode { Values, Refs, RefBases };
        OutputMode outputMode_;

  edm::Service<edm::RandomNumberGenerator> rng_;
  bool useProb_;
  
};

template<typename T>
ResonanceInefficiencyCreator<T>::ResonanceInefficiencyCreator(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    tags_(iConfig.getParameter<edm::InputTag>("tags")),
    mass_(iConfig.getParameter<double>("mass")),
    massMin_(iConfig.getParameter<double>("massMin")),
    massMax_(iConfig.getParameter<double>("massMax")),
    onlyBestMatch_(iConfig.getParameter<bool>("onlyBestMatch")),
    probability_(iConfig.existsAs<std::string>("probability") ? iConfig.getParameter<std::string>("probability") : "1")
{
  useProb_ =false;
    if (iConfig.existsAs<std::string>("probability")) {
      useProb_ =true;
      if ( ! rng_.isAvailable()) {
	throw cms::Exception("Configuration")
	  << "XXXXXXX requires the RandomNumberGeneratorService\n"
	  "which is not present in the configuration file.  You must add the service\n"
	  "in the configuration file or remove the modules that require it.";
      }

        //CLHEP::HepRandomEngine& engine = rng->getEngine();

        // engine MUST be a reference here, if a pointer is used the
        // distribution will destroy the engine in its destructor, a major
        // problem because the service owns the engine and will destroy it 
	//        flatDistribution_.reset(new CLHEP::RandFlat(engine, 0, 1));
    }

    std::string outputMode = iConfig.getParameter<std::string>("outputMode");
    if (outputMode == "vector") {
        outputMode_ = Values;
        produces<PlainVecT>(); 
    } else if (outputMode == "RefVector") {
        outputMode_ = Refs;
        produces<RefVecT>(); 
    } else if (outputMode == "RefToBaseVector") {
        outputMode_ = RefBases;
        produces<RefBaseVecT>();
    }
}

template<typename T>
ResonanceInefficiencyCreator<T>::~ResonanceInefficiencyCreator() 
{
}

template<typename T>
void
ResonanceInefficiencyCreator<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<T> > src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<edm::View<reco::Candidate> > tags; 
    iEvent.getByLabel(tags_, tags);

    CLHEP::HepRandomEngine& engine = rng_->getEngine( iEvent.streamID() );

    // Prepare output
    std::auto_ptr<PlainVecT>   vec;
    std::auto_ptr<RefVecT>     refvec;
    std::auto_ptr<RefBaseVecT> rbvec;
    switch (outputMode_) {
        case Values:   vec.reset(new PlainVecT());      break;
        case Refs:     refvec.reset(new RefVecT());     break;
        case RefBases: rbvec.reset(new RefBaseVecT());  break;
    }

    // Filter and fill
    size_t i, n = src->size(); 
    std::vector<uint8_t> marked(n, 0);
    edm::View<reco::Candidate>::const_iterator it, ed;
    for (it = tags->begin(), ed = tags->end(); it != ed; ++it) {
        const reco::Candidate & tag = *it;
            
        int match = -1, j;
        double bestMassDiff = massMax_ - massMin_;
        typename edm::View<T>::const_iterator itp, edp;
        for (itp = src->begin(), edp = src->end(), j = 0; itp != edp; ++itp, ++j) {
            const reco::Candidate & probe = *itp;
            double mass = (tag.p4()+probe.p4()).M();
            if (massMin_ < mass && mass < massMax_) {
                if (onlyBestMatch_) {
                    if (fabs(mass - mass_) < bestMassDiff) {
                        match = j;
                        bestMassDiff = fabs(mass - mass_);
                    }
                } else {
                    marked[j] = 1;
                }
            }
        }
        if (onlyBestMatch_ && (match != -1)) marked[match] = 1;
    }

    // then remove marked elements
    for (i = 0; i < n; ++i) {
        const T &t = (*src)[i];
        bool  ok = true;
        if (marked[i]) { // it's targeted
	  if (useProb_) { //flatDistribution_.get()
	   // random chance
	    ok = engine.flat() <= probability_(t);
	  } else {
                // kill all
                ok = false;
            }
        }        
        if (ok) {
            switch (outputMode_) {
                case Values:   vec->push_back(t); break;
                case Refs:     refvec->push_back(src->refAt(i).template castTo<RefT>()); break;
                case RefBases: rbvec->push_back(src->refAt(i));  break;
            }
        }
    }

    /// Write out output
    switch (outputMode_) {
        case Values:   iEvent.put(vec);    break;
        case Refs:     iEvent.put(refvec); break;
        case RefBases: iEvent.put(rbvec);  break;
    }
}

typedef ResonanceInefficiencyCreator<reco::Candidate>   CandidateResonanceInefficiencyCreator;
DEFINE_FWK_MODULE(CandidateResonanceInefficiencyCreator);
