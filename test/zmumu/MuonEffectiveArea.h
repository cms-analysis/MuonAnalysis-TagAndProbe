//--------------------------------------------------------------------------------------------------
// $Id $
//
// MuonEffectiveArea
//
// Helper Class for storing effective areas
//
// Authors: S.Xie, E. DiMarco, L. Perrozzi
//--------------------------------------------------------------------------------------------------


/// --> NOTE if you want to use this class as standalone without the CMSSW part 
///  you need to uncomment the below line and compile normally with scramv1 b 
///  Then you need just to load it in your root macro the lib with the correct path
//#define STANDALONE   // <---- this line

#ifndef MuonEffectiveArea_H
#define MuonEffectiveArea_H

#ifndef STANDALONE
#endif

using namespace std;

class MuonEffectiveArea{
public:
  MuonEffectiveArea();
  ~MuonEffectiveArea(); 
  
  enum MuonEffectiveAreaType {
    kMuMiniIso03
  };
  
  enum MuonEffectiveAreaTarget {
    kMuEANoCorr,
    kMuEASpring15_25ns
  };

  static Double_t GetMuonEffectiveArea(MuonEffectiveAreaType type, Double_t MuEta, 
  MuonEffectiveAreaTarget EffectiveAreaTarget = kMuEASpring15_25ns) {
    
    Double_t EffectiveArea = 0;
    
    
    if (EffectiveAreaTarget == kMuEANoCorr) {
      return 0.0;
    }
    
    //2015 Spring15 Effective Areas
    else if (EffectiveAreaTarget == kMuEASpring15_25ns) {
      if (type == kMuMiniIso03){
        if (fabs(MuEta) >= 0.0 && fabs(MuEta) < 0.8 )   EffectiveArea = 0.0735;
        if (fabs(MuEta) >= 0.8 && fabs(MuEta) < 1.3 ) EffectiveArea = 0.0619;
        if (fabs(MuEta) >= 1.3 && fabs(MuEta) < 2.0 ) EffectiveArea = 0.0465;
        if (fabs(MuEta) >= 2.0 && fabs(MuEta) < 2.2 )   EffectiveArea = 0.0433;
        if (fabs(MuEta) >= 2.2 && fabs(MuEta) <= 2.5 )   EffectiveArea = 0.0577;
      }

    }
    
    return EffectiveArea;  
  }
};

#endif
