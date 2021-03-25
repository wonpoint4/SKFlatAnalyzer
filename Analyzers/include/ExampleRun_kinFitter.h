#ifndef ExampleRun_kinFitter_h
#define ExampleRun_kinFitter_h

#include "SMPAnalyzerCore.h"
#include "TKinFitterDriver.h"

class ExampleRun_kinFitter : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void GetTTLJGenParticles(const vector<Gen>& gens, Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& b0,Gen& b1,Gen& j0,Gen& j1, int mode);
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys, MuonTrigSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  double weight_Prefire;

  ExampleRun_kinFitter();
  ~ExampleRun_kinFitter();

private:

  TKinFitterDriver* fitter;

};



#endif

