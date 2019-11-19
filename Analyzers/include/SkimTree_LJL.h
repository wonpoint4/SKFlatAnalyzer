#ifndef SkimTree_LJL_h
#define SkimTree_LJL_h

#include "SMPAnalyzerCore.h"

class SkimTree_LJL : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_LJL();
  ~SkimTree_LJL();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif

