#ifndef ExampleRun_kinFitter_h
#define ExampleRun_kinFitter_h

#include "SMPAnalyzerCore.h"
#include "dybAnalyzer.h"
#include "TKinFitterDriver.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

class ExampleRun_kinFitter : public dybAnalyzer {

 public:

  virtual void initializeAnalyzer();
  virtual void GetTTLJGenParticles(const vector<Gen>& gens, Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& b0,Gen& b1,Gen& j0,Gen& j1, int mode);
  virtual void executeEventWithParameter(TString channel);
  virtual void executeEvent();
  virtual bool Hasleptons(TString channel);
  void SetupJetResolution();
  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;

  Lepton* lepton0 = NULL;
  Jet* bjet0 = NULL;
  Jet* bjet1 = NULL;
  Jet* ajet0 = NULL;
  Jet* ajet1 = NULL;
  double bcharge0 = 0;
  double bcharge1 = 0;
  double acharge0 = 0;
  double acharge1 = 0;
  TLorentzVector met;

  ExampleRun_kinFitter();
  ~ExampleRun_kinFitter();

 private:

  TKinFitterDriver* fitter;

};



#endif
