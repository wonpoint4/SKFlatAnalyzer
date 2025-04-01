#ifndef SkimTree_1B2L_h
#define SkimTree_1B2L_h

#include "AnalyzerCore.h"

class SkimTree_1B2L : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_1B2L();
  ~SkimTree_1B2L();

  TTree *newtree;

  vector<TString> double_triggers;
  vector<TString> single_muon_triggers;
  vector<TString> single_electron_triggers;
  void WriteHist();

};



#endif
