#ifndef MuonAnalysis_TagAndProbe_TagProbeEDMFilter_h
#define MuonAnalysis_TagAndProbe_TagProbeEDMFilter_h
// -*- C++ -*-
//
// Package:     TagAndProbe
// Class  :     TagProbeEDMFilter
// 
/**\class TagProbeEDMFilter TagProbeEDMFilter.h MuonAnalysis/TagAndProbe/interface/TagProbeEDMFilter.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author: Nadia Adam (Princeton) 
//         Created:  Fri Jun  6 09:13:10 CDT 2008
// $Id$
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class TagProbeEDMFilter : public edm::EDFilter 
{
   public:
      explicit TagProbeEDMFilter(const edm::ParameterSet&);
      ~TagProbeEDMFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
};
#endif
