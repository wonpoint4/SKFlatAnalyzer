#ifndef tWAnalyzer_h
#define tWAnalyzer_h

#include "SMPAnalyzerCore.h"

class tWAnalyzer : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  double Lepton_SF(TString histkey,const Lepton* lep,int sys);
  double LeptonTrigger_SF(TString triggerSF_key,const vector<Lepton*>& leps,int sys);

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys, MuonTrigSFKeys;
  double weight_Prefire;

  tWAnalyzer();
  ~tWAnalyzer();

};



#endif

