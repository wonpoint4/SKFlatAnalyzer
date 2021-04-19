#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "SMPAnalyzerCore.h"

class EfficiencyValidation : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void FillHists(Parameter& p);
  void FillHistsEfficiency(Parameter& p,TString region);

  EfficiencyValidation();
  ~EfficiencyValidation();

};



#endif

