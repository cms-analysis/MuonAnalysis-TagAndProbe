#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"

float MuonPerformance::getResult(PerformanceResult::ResultType r, BinningPointByMap p) const {
  return pl.getResult(r,p);
}

bool MuonPerformance::isResultOk(PerformanceResult::ResultType r, BinningPointByMap p) const {
  return pl.isInPayload(r,p);
}


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/Utilities/interface/typelookup.h"

TYPELOOKUP_DATA_REG(MuonPerformance);
