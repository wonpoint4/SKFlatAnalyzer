#ifndef ExampleRun_Higgs_h
#define ExampleRun_Higgs_h

#include "AnalyzerCore.h"

class ExampleRun_Higgs : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  ExampleRun_Higgs();
  ~ExampleRun_Higgs();

  vector<Muon> muons;
  vector<Electron> electrons;
  vector<Photon> photons;

  void FillHists_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons, double weight, TString prefix="");
  bool Charge_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons);
  bool Mass_4leptons(const vector<Muon>& muons, const vector<Electron>& electrons, double mass_cut);
};



#endif

