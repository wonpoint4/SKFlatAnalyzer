#ifndef AFBAnalyzerSyst_h
#define AFBAnalyzerSyst_h

#include "AFBAnalyzer.h"

class AFBAnalyzerSyst : public AFBAnalyzer {

public:
  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual void executeEventGen();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual Variations MakeVariations(const Parameter& p);
  virtual void FillHistsSyst(Parameter p,Variations& vs);
  virtual void FillHistsUnfold(Parameter& preco,Parameter& pgen);
  virtual bool PassSelection(Parameter& p,bool cutflow=false);
  using SMPAnalyzerCore::executeEventWithParameter;
  virtual void executeEventWithParameter(Parameter& p){
    SMPAnalyzerCore::executeEventWithParameter(p);
  }
  //virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;executeEventWithParameter(pp);}

  AFBAnalyzerSyst();
  ~AFBAnalyzerSyst();

};



#endif

