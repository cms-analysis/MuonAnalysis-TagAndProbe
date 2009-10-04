
#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h"

// constructors and destructor
//
MuonPerformanceReadback::MuonPerformanceReadback()
{
  rnd = new TRandom3();
}

MuonPerformanceReadback::~MuonPerformanceReadback()
{
}

//
// member functions
//
double
MuonPerformanceReadback::getEff(double pt, double eta, double phi, int charge, const MuonPerformance & theperf)
{
  BinningPointByMap p; 
  p.insert(BinningVariables::MuonPt,pt);  
  p.insert(BinningVariables::MuonEta,eta); 
  //  p.insert(BinningVariables::MuonPhi,phi); 
  //  p.insert(BinningVariables::MuonCharge,charge);  

  double theeff = theperf.getResult(PerformanceResult::MUEFF,p);  

  /* Debugging printout. */ 
  std::cout <<" Discriminant is "<<theperf.workingPoint().discriminantName()<<std::endl;   
  std::cout <<"\teta=" << eta << ", pt=" << pt << ", phi = " << phi << ", " << "charge = " << charge << std::endl;   
  std::cout <<"\tmueff/muerr ="<<theperf.getResult(PerformanceResult::MUEFF,p)<<"/"<<theperf.getResult(PerformanceResult::MUERR,p)<<std::endl;  
 
  return theeff;
}

double 
MuonPerformanceReadback::getEffError(double pt, double eta, double phi, int charge, const MuonPerformance & theperf) 
{ 
  BinningPointByMap p;  
  p.insert(BinningVariables::MuonPt,pt);   
  p.insert(BinningVariables::MuonEta,eta);  
  //  p.insert(BinningVariables::MuonPhi,phi);  
  //  p.insert(BinningVariables::MuonCharge,charge);   
 
  double theefferr = theperf.getResult(PerformanceResult::MUERR,p);  
  
  return theefferr; 
} 

bool
MuonPerformanceReadback::passesPIDKilling(double pt, double eta, double phi, int charge, const MuonPerformance & theperf)
{
  bool passes = false;

  BinningPointByMap p;   
  p.insert(BinningVariables::MuonPt,pt);    
  p.insert(BinningVariables::MuonEta,eta);   
  //  p.insert(BinningVariables::MuonPhi,phi);   
  //  p.insert(BinningVariables::MuonCharge,charge);    
  
  double theeff = theperf.getResult(PerformanceResult::MUEFF,p);   

  double randthrow = rnd->Rndm(); 
 
  if(randthrow < theeff) 
    { 
      passes = true;
    }

  return passes;
}

MuonPerformance 
MuonPerformanceReadback::getPerformanceRecord(std::string name, const edm::EventSetup& eventSetup)
{
  edm::ESHandle<MuonPerformance> performanceRecord; 
  eventSetup.get<MuonPerformanceRecord>().get(name,performanceRecord); 
  const MuonPerformance & thePerformance = *(performanceRecord.product()); 

  return thePerformance;
}
