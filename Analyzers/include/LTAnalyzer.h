#ifndef LTAnalyzer_h
#define LTAnalyzer_h

#include "SMPAnalyzerCore.h"

class LTAnalyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual void executeEventWithParameter(Parameter& p);
  virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;executeEventWithParameter(pp);}
  virtual pair<double,double> GetCostAndPhiCS(Particle* l0,Particle* l1);
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual void EvalWeights(Parameter& p);
  virtual void ResetRecoWeights(Parameter& p);
  virtual void FillHists(Parameter& p);
  static int GetUnfoldBin(int njet,double mass,double pt,double cost,double phi);

  LTAnalyzer();
  ~LTAnalyzer();

  static const int mbinnum=42;
  static constexpr const double mbin[mbinnum+1]={52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,1500,2000,3000};
  static const int njetbin=2;
  static const int nptbin=6;
  static constexpr const double ptbins[nptbin+1]={0,20,30,50,70,100,200};
  static const int ncostbin=10;
  static constexpr const double costbins[ncostbin+1]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  static const int nphibin=10;
  static constexpr const double phibins[nphibin+1]={0.0*TMath::Pi(),0.05*TMath::Pi(),0.1*TMath::Pi(),0.15*TMath::Pi(),0.2*TMath::Pi(),0.25*TMath::Pi(),0.3*TMath::Pi(),0.35*TMath::Pi(),0.4*TMath::Pi(),0.45*TMath::Pi(),0.5*TMath::Pi()};
  static const int nresponsebin=2*nptbin*ncostbin*nphibin;

  
};




#endif

