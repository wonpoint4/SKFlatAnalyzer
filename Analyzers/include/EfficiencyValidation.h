#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "SMPAnalyzerCore.h"

class EfficiencyValidation : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventWithParameter(Parameter p);
  void executeEvent();
  void FillHistsEfficiency(TString pre,TString suf,const vector<Lepton*>& leps,const map<TString,double>& weights);
  double LeptonTrigger_SF_OR(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys);

  EfficiencyValidation();
  ~EfficiencyValidation();

};



#endif

