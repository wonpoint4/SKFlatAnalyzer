#ifndef testGen_h
#define testGen_h

#include "AnalyzerCore.h"

class testGen : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<Muon> AllMuons;
  vector<Gen> AllGens;

  testGen();
  ~testGen();

};



#endif

