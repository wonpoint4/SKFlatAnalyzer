#ifndef SKFlatValidation_h
#define SKFlatValidation_h

#include "AnalyzerCore.h"

class SKFlatValidation : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);

  void executeEvent();

  double weight_Prefire;

  vector<TString> Triggers_POG_Electron, Triggers_POG_Muon;
  TString TriggerNameForSF_POG_Electron, TriggerNameForSF_POG_Muon;
  double TriggerSafePt_POG_Electron, TriggerSafePt_POG_Muon;

  vector<TString> Triggers_POGHighPt_Electron, Triggers_POGHighPt_Muon;
  TString TriggerNameForSF_POGHighPt_Electron, TriggerNameForSF_POGHighPt_Muon;
  double TriggerSafePt_POGHighPt_Electron, TriggerSafePt_POGHighPt_Muon;
  double DoublePhotonSafePtCut;

};



#endif

