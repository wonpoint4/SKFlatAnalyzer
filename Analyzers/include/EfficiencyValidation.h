#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "SMPAnalyzerCore.h"

class EfficiencyValidation : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventWithChannelName(TString channelname);
  void executeEvent();
  void FillHistsEfficiency(TString pre,TString suf,const vector<Lepton*>& leps,const map<TString,double>& weights);

  EfficiencyValidation();
  ~EfficiencyValidation();

};



#endif

