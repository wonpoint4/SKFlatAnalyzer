#ifndef SkimTree_SS2lOR3l_h
#define SkimTree_SS2lOR3l_h

#include "AnalyzerCore.h"

class SkimTree_SS2lOR3l : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_SS2lOR3l();
  ~SkimTree_SS2lOR3l();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif

