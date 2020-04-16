#ifndef ZptWeight_h
#define ZptWeight_h

#include "SMPAnalyzerCore.h"

class ZptWeight : public SMPAnalyzerCore {

public:

  void executeEventFromParameter(TString channelname,Event* ev);
  void executeEvent();
  void FillHistsZptWeight(TString pre,TString suf,const vector<Lepton*>& leps,const map<TString,double>& weights);

  ZptWeight();
  ~ZptWeight();

  static const int ptbinnum=42;
  static const int rapbinnum=6;
  double ptbin[ptbinnum+1]={0,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,34,36,38,40,44,48,52,56,62,70,80,100,400};
  double rapbin[rapbinnum+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
};



#endif

