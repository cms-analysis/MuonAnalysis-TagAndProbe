#include "MuonAnalysis/TagAndProbe/plugins/MuonPerformanceESProducer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"


#include <memory>
#include <string>

using namespace edm;

MuonPerformanceESProducer::MuonPerformanceESProducer(const edm::ParameterSet & p) 
{
  std::string myname = p.getParameter<std::string>("ComponentName");
  mypl = p.getParameter<std::string>("PayloadName"); 
  mywp = p.getParameter<std::string>("WorkingPointName");
  
  pset_ = p;
  setWhatProduced(this,myname);
}

MuonPerformanceESProducer::~MuonPerformanceESProducer() {}

boost::shared_ptr<MuonPerformance> 
MuonPerformanceESProducer::produce(const MuonPerformanceRecord & iRecord){ 
   ESHandle<PerformancePayload> pl;
   //ESHandle<PhysicsPerformancePayload> pl;
   ESHandle<PerformanceWorkingPoint> wp;
   iRecord.getRecord<PerformancePayloadRecord>().get(mypl,pl);
   
   std::cout <<"HERE "<<std::endl;
   iRecord.getRecord<PerformanceWPRecord>().get(mywp,wp);
   std::cout <<"HERE "<<std::endl;
   
   std::cout <<" Got the payload, which is a  "<<typeid(*(pl.product())).name()<<std::endl;
   
   _perf  = boost::shared_ptr<MuonPerformance>(new MuonPerformance(*((pl.product())), *((wp.product()))));

   return _perf;
}


DEFINE_FWK_EVENTSETUP_MODULE(MuonPerformanceESProducer);

