#ifndef MuonPerformance_h
#define MuonPerformance_h

#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayload.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformanceWorkingPoint.h"


#include <string>
#include <vector>

class MuonPerformance {
public:
  MuonPerformance(const PerformancePayload& p, const PerformanceWorkingPoint& w) : pl(p), wp(w) {}

  virtual float getResult(PerformanceResult::ResultType, BinningPointByMap) const ;

  virtual bool isResultOk(PerformanceResult::ResultType, BinningPointByMap) const ;
  
  virtual const PerformanceWorkingPoint& workingPoint() const {return wp;}

  private:
  const PerformancePayload& pl;
  const PerformanceWorkingPoint& wp;

};


#endif

