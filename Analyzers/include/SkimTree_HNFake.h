#ifndef SkimTree_HNFake_h
#define SkimTree_HNFake_h

#include "AnalyzerCore.h"

class SkimTree_HNFake : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNFake();
  ~SkimTree_HNFake();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif
