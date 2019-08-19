#ifndef AFBAnalyzer_h
#define AFBAnalyzer_h

#include "SMPAnalyzerCore.h"

class AFBAnalyzer : public SMPAnalyzerCore {

public:

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
  static const int massbinnum=16;
  const double massrange[massbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120,150,200,280,400};
  void FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w);
  void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
  TRandom3* random;
  TString hardprefix;

};



#endif

