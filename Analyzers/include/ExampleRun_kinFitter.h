#ifndef ExampleRun_kinFitter_h
#define ExampleRun_kinFitter_h

#include "SMPAnalyzerCore.h"
#include "dybAnalyzer.h"
#include "TKinFitterDriver.h"

class ExampleRun_kinFitter : public dybAnalyzer {

 public:

  virtual void initializeAnalyzer();
  virtual void GetTTLJGenParticles(const vector<Gen>& gens, Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& b0,Gen& b1,Gen& j0,Gen& j1, int mode);
  virtual void executeEventWithParameter(TString channel);
  virtual void executeEvent();
  virtual bool Hasleptons(TString channel);

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  vector<TString> MuonIDs, MuonIDSFKeys, MuonTrigSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  TLorentzVector met;
  double weight_Prefire;
  double btagweight;

  ExampleRun_kinFitter();
  ~ExampleRun_kinFitter();

 private:

  TKinFitterDriver* fitter;

};



#endif
