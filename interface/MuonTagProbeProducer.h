#ifndef MuonAnalysis_TagAndProbe_MuonTagProbeProducer_h
#define MuonAnalysis_TagAndProbe_MuonTagProbeProducer_h
// -*- C++ -*-
//
// Package:     TagAndProbe
// Class  :     MuonTagProbeProducer
// 
/**\class MuonTagProbeProducer MuonTagProbeProducer.h MuonAnalysis/TagAndProbe/interface/MuonTagProbeProducer.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Wed Mar  5 11:07:48 CST 2008
// $Id$
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"



//
// class decleration
//

class MuonTagProbeProducer : public edm::EDProducer 
{
   public:
      explicit MuonTagProbeProducer(const edm::ParameterSet&);
      ~MuonTagProbeProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      
      edm::InputTag tagCollection_;
      edm::InputTag probeCollection_;

      double massMinCut_;
      double massMaxCut_;
      double delRMinCut_;
      double delRMaxCut_;

      bool requireOS_;
};


#endif
