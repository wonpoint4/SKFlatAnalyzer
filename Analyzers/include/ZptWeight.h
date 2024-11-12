#ifndef ZptWeight_h
#define ZptWeight_h

#include "SMPAnalyzerCore.h"

class ZptWeight : public SMPAnalyzerCore {

public:

  virtual void executeEvent();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual Variations MakeVariations(const Parameter& p);
  virtual void FillHists(Parameter& p);

  ZptWeight();
  ~ZptWeight();

  static const int massbinnum=13;
  static constexpr const double massbin[massbinnum+1]={52,60,70,77,86,96,106,120,150,200,280,400,800,3000};
  static const int ybinnum=14;
  static constexpr const double ybin[ybinnum+1]={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.8,4.8};
  static const int ptbinnum=52;
  static constexpr const double ptbin[ptbinnum+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,62,70,80,100,120,140,160,180,200,240,280,320,360,400,450,500,550,650};
};
#endif

