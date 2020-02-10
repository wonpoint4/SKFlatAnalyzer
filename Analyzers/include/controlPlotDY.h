#ifndef controlPlotDY_h
#define controlPlotDY_h

#include "AnalyzerCore.h"

class controlPlotDY : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  double GetCosineCS(Muon mu1, Muon mu2);

  vector<TString> MuonIDISOSFs, MuonTrigSFs, Samples;
  bool IsSingleMuonData, IsSingleMuon;
  bool IsDoubleMuonData, IsDoubleMuon;
  bool IsDYsamples;

  TString IsFromTau;
  vector<Muon> AllMuons;
  vector<Gen> gens;

  double weight_Prefire;

  controlPlotDY();
  ~controlPlotDY();

};



#endif

