#ifndef ScaledPts_h
#define ScaledPts_h

struct ScaledPts {
  double ScaledPt, ScaeldPt_Up, ScaeldPt_Down;
};

#endif


#ifndef GeneralizedEndpoint_h
#define GeneralizedEndpoint_h

#include <iostream>
#include <map>

//==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceSelectionAndCalibrationsRun2
//==== https://github.com/cms-analysis/SUSYBSMAnalysis-Zprime2muAnalysis/blob/mini-AOD/src/GeneralizedEndpoint.cc

using namespace std;

class GeneralizedEndpoint{

public:
  GeneralizedEndpoint();
  virtual ~GeneralizedEndpoint();

  ScaledPts GeneralizedEndpointPt(float MuonPt, int MuonCharge, float MuonEta, float MuonPhi, int seed);

private:
  std::map<int,std::map<int,float> > _Correction;
  std::map<int,std::map<int,float> > _CorrectionError;   

};

#endif
