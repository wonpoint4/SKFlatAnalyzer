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
  void FillHistAll(TString channelname,TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight);
  void FillAFBSystematicHists(TString pre,TString suf,Particle* l0,Particle* l1,map<TString,double> map_weight);
  void FillAFBHists(TString pre,TString suf,Particle* l0,Particle* l1,double w);
  //static const int massbinnum=12;
  //const double massrange[massbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120};
  void FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w);
  void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
  TRandom3* random;
  TString hardprefix;

  static const int mbinnum=16;
  const double mbin[mbinnum+1]={40,60,70,80,88,91,94,100,110,120,150,200,280,400,800,2000,5000};
  static const int mbinnum_fine=32;
  const double mbin_fine[mbinnum_fine+1]={40,50,60,65,70,75,80,84,88,91,94,97,100,105,110,115,120,130,140,150,175,200,240,280,340,400,600,800,1400,2000,3000,4000,5000};
  static const int ybinnum=7;
  const double ybin[ybinnum+1]={0,0.4,0.8,1.2,1.6,2.0,2.4,2.8};
  static const int ptbinnum=8;
  const double ptbin[ptbinnum+1]={0,10,20,40,60,100,200,400,1000};
  static const int ptbinnum_fine=33;
  const double ptbin_fine[ptbinnum_fine+1]={0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,180,200,300,400,600,800,1000};
  static const int lptbinnum=56;
  const double lptbin[lptbinnum+1]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,64,68,72,76,80,85,90,95,100,110,120,130,140,150,160,180,200,250,300,350,400,500,600,700,800,900,1000};
};



#endif

