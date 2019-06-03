#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "AFBAnalyzer.h"

class EfficiencyValidation : public AFBAnalyzer {

public:

  void executeEventFromParameter(TString channelname,Event* ev);
  void executeEvent();
  void FillHistsEfficiency(TString pre,TString suf,const vector<Lepton*>& leps,const map<TString,double>& weights);

  EfficiencyValidation();
  ~EfficiencyValidation();

};



#endif

