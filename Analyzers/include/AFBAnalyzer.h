#ifndef AFBAnalyzer_h
#define AFBAnalyzer_h

#include "SMPAnalyzerCore.h"

class AFBAnalyzer : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(TString channelname,Event* ev);

  AFBAnalyzer();
  ~AFBAnalyzer();

  //double GetCosThetaCS(const vector<Lepton*>& leps);
  double GetCosThetaCS(const Particle *p0,const Particle *p1);
  double GetCosTheta(const vector<Lepton*>& leps,const vector<Jet>& jets,TString option,double fcut);
  void FillHists(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight);
  void FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w);
  void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
  TRandom3* random;
  TString hardprefix;

  static const int afb_mbinnum=39;
  const double afb_mbin[afb_mbinnum+1]={60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,600,800,1400,2000,3000,4000,5000};
  static const int afb_ybinnum=7;
  const double afb_ybin[afb_ybinnum+1]={0,0.4,0.8,1.2,1.6,2.0,2.4,2.8};
  static const int afb_ptbinnum=33;
  const double afb_ptbin[afb_ptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,180,200,300,400,600,800,1000};

  static const int grid_mbinnum=5;
  const double grid_mbin[grid_mbinnum+1]={60,80,100,120,500,5000};
  static const int grid_ybinnum=2;
  const double grid_ybin[grid_ybinnum+1]={0,1.2,2.4};
  static const int grid_ptbinnum=4;
  const double grid_ptbin[grid_ptbinnum+1]={0,20,50,100,1000};

  static const int fine_mbinnum=39;
  const double fine_mbin[fine_mbinnum+1]={60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,600,800,1400,2000,3000,4000,5000};
  static const int fine_ptbinnum=33;
  const double fine_ptbin[fine_ptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,180,200,300,400,600,800,1000};
  
  static const int lptbinnum=56;
  const double lptbin[lptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,64,68,72,76,80,85,90,95,100,110,120,130,140,150,160,180,200,250,300,350,400,500,600,700,800,900,1000};
};



#endif

