#ifndef ControlplotISR_h
#define ControlplotISR_h

#include "AnalyzerCore.h"

class ControlplotISR : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  TString IsoMuTriggerName;
  vector<Muon> AllMuons;
  double weight_Prefire;

  ControlplotISR();
  ~ControlplotISR();

};



#endif

