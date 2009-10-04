#ifndef MuonPerformanceReadback_h
#define MuonPerformanceReadback_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoMuon/Records/interface/MuonPerformanceRecord.h" 
#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h"  
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"  
#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"  
#include <TRandom3.h>

class MuonPerformanceReadback {
public:
  MuonPerformanceReadback();
  ~MuonPerformanceReadback();

  double getEff(double, double, double, int, const MuonPerformance &);
  double getEffError(double, double, double, int, const MuonPerformance &); 
  bool passesPIDKilling(double, double, double, int, const MuonPerformance &);

  MuonPerformance getPerformanceRecord(std::string, const edm::EventSetup& eventSetup);

 private:
  TRandom3 *rnd; 

};

#endif
