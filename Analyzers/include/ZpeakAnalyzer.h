#ifndef ZpeakAnalyzer_h
#define ZpeakAnalyzer_h

#include "SMPAnalyzerCore.h"

class ZpeakAnalyzer : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void FillHists(Parameter& p);

  ZpeakAnalyzer();
  ~ZpeakAnalyzer();

  static const int nptbin=16;
  double ptbins[nptbin+1]={10,15,20,25,30,35,40,45,50,55,60,80,100,150,200,400,1000};
  static const int netabin=26;
  double etabins[netabin+1]={-2.5,-2.4,-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.5};
  static const int nmassbin=98;
  double massbins[nmassbin+1]={52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150};

};



#endif

