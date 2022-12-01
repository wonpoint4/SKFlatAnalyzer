#ifndef FakeAnalyzer_h
#define FakeAnalyzer_h

#include "SMPAnalyzerCore.h"

class FakeAnalyzer : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void FillHists(Parameter& p);
  void FillFakeHists(Parameter& p,TString region);

  FakeAnalyzer();
  ~FakeAnalyzer();

  static const int nptbin=24;
  double ptbins[nptbin+1]={10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60,70,80,90,100,120,150,200,300,400,500,700,1000};
  static const int netabin=25;
  double etabins[netabin+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};

};



#endif

