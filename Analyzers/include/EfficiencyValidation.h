#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "SMPAnalyzerCore.h"

class EfficiencyValidation : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual void EvalWeights(Parameter& p);
  virtual void FillHists(Parameter& p);
  virtual void FillHistsEfficiency(Parameter& p,TString region);
  virtual bool PassSelection(Parameter& p);

  EfficiencyValidation();
  ~EfficiencyValidation();

  static const int mbinnum=42;
  const double mbin[mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,1500,2000,3000};
};



#endif

