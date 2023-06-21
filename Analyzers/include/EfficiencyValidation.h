#ifndef EfficiencyValidation_h
#define EfficiencyValidation_h

#include "SMPAnalyzerCore.h"

class EfficiencyValidation : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual void EvalWeights(Parameter& p);
  virtual void FillHists(Parameter& p);
  virtual void FillHistsEfficiency(Parameter& p,TString region);
  virtual bool PassSelection(Parameter& p);

  virtual double GetLeptonTriggerORSF_old(const Parameter& p,int iset=0,int imem=0);
  virtual double GetCosThetaCS(const Particle *p0,const Particle *p1,int direction=0) const;

  EfficiencyValidation();
  ~EfficiencyValidation();

  static const int mbinnum=42;
  const double mbin[mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,1500,2000,3000};
  static const int rochester_nmbin=11;
  const double rochester_mbins[rochester_nmbin+1]={54,66,76,82,86,89.5,92.7,96,100,106,116,150};
  static const int rochester_nybin=6;
  const double rochester_ybins[rochester_nybin+1]={0,0.4,0.8,1.2,1.6,2.0,2.4};
  static const int rochester_ncostbin=10;
  const double rochester_costbins[rochester_ncostbin+1]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

  static const int netabin_muonID=48;
  const double etabins_muonID[netabin_muonID+1]={-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
  static const int nptbin_muonID=10;
  const double ptbins_muonID[nptbin_muonID+1]={10, 15, 20, 25, 30, 35, 40, 45, 50, 70, 200};
  //static const int netabin_electronID=48;
  //const double etabins_electronID[netabin_electronID+1]={-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.566,-1.4442,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4442,1.566,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  static const int netabin_electronID=50;
  const double etabins_electronID[netabin_electronID+1]={-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  static const int nptbin_electronID=12;
  const double ptbins_electronID[nptbin_electronID+1]={10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 500};
};



#endif

