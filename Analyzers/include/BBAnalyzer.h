#ifndef BBAnalyzer_h
#define BBAnalyzer_h

#include "SMPAnalyzerCore.h"

class BBAnalyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual void EvalDefaultWeight(Parameter& p);
  virtual Variations MakeVariations(const Parameter& p);
  virtual void EvalVariationsBcharge(const Parameter& p,Variations& v);
  virtual void FillHists(Parameter& p);

  bool IsDileptonSkim;

  BBAnalyzer();
  ~BBAnalyzer();

  static const int nptbin=17;
  double ptbins[nptbin+1]={0,10,15,20,25,30,35,40,45,50,55,60,80,100,150,200,400,1000};
  static const int netabin=26;
  double etabins[netabin+1]={-2.5,-2.4,-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.5};
  static const int nmassbin=40;
  double massbins[nmassbin+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};


};



#endif

