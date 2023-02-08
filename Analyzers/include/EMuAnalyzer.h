#ifndef EMuAnalyzer_h
#define EMuAnalyzer_h

#include "SMPAnalyzerCore.h"

class EMuAnalyzer : public SMPAnalyzerCore {

public:

  virtual void executeEvent();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual void EvalWeights(Parameter& p);
  virtual void FillHists(Parameter& p);

  EMuAnalyzer();
  ~EMuAnalyzer();

  static const int mbinnum=42;
  const double mbin[mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,1500,2000,3000};

  static const int nptbin=28;
  double ptbins[nptbin+1]={0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60,70,80,90,100,120,150,200,300,400,500,700,1000};
  static const int nmbin=28;
  double mbins[nmbin+1]={52,56,60,65,70,75,80,85,90,95,100,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};

};



#endif

