#ifndef FakeAnalyzer_h
#define FakeAnalyzer_h

#include "SMPAnalyzerCore.h"

class FakeAnalyzer : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void FillHists(Parameter& p);
  void FillDileptonHists(Parameter& p,TString region,TString hprefix="",TString suffix="",double w=1.);
  int GetWP(const Lepton* lep);

  FakeAnalyzer();
  ~FakeAnalyzer();

  static const int nptbin=24;
  double ptbins[nptbin+1]={10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60,70,80,90,100,120,150,200,300,400,500,700,1000};
  static const int netabin=25;
  double etabins[netabin+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  static const int nmbin=47;
  double mbins[nmbin+1]={20,25,30,35,40,45,50,55,60,65,70,73,76,78,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};

};



#endif

