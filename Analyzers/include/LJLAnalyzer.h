#ifndef LJLAnalyzer_h
#define LJLAnalyzer_h

#include "AFBAnalyzer.h"

class LJLAnalyzer : public AFBAnalyzer {

public:

  void executeEventFromParameter(TString channelname,Event* ev);
  void executeEvent();
  void FillLJLHists(TString pre,TString suffix,const vector<Lepton*>& isoleps,const vector<Lepton*>& nonisoleps,const vector<Jet>& jets,const map<TString,double>& weights);

  LJLAnalyzer();
  ~LJLAnalyzer();

};



#endif

