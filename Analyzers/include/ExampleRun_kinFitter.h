#ifndef ExampleRun_kinFitter_h
#define ExampleRun_kinFitter_h

#include "AnalyzerCore.h"
#include "TKinFitterDriver.h"

class ExampleRun_kinFitter : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  double weight_Prefire;

  ExampleRun_kinFitter();
  ~ExampleRun_kinFitter();

private:

  TKinFitterDriver* fitter;

};



#endif

