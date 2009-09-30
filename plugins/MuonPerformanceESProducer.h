#ifndef MuonAnalysis_TagAndProbe_MuonPerformanceESProducer_H
#define MuonAnalysis_TagAndProbe_MuonPerformanceESProducer_H

#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"
#include "RecoMuon/Records/interface/MuonPerformanceRecord.h"


#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <boost/shared_ptr.hpp>

class  MuonPerformanceESProducer : public edm::ESProducer{
 public:
  MuonPerformanceESProducer(const edm::ParameterSet & p);
  virtual ~MuonPerformanceESProducer(); 
  boost::shared_ptr<MuonPerformance> produce(const  MuonPerformanceRecord &);
 private:
  boost::shared_ptr<MuonPerformance> _perf;
  edm::ParameterSet pset_;
  std::string mypl;
  std::string mywp;
};


#endif




