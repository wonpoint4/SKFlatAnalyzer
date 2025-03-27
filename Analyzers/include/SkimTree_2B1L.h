#ifndef SkimTree_2B1L_h
#define SkimTree_2B1L_h

#include "AnalyzerCore.h"

class SkimTree_2B1L : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_2B1L();
  ~SkimTree_2B1L();

  TTree *newtree;

  vector<TString> single_muon_triggers;
  vector<TString> single_electron_triggers;
  void WriteHist();

};



#endif
