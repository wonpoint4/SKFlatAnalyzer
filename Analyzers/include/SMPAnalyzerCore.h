#ifndef SMPAnalyzerCore_h
#define SMPAnalyzerCore_h

#include "AnalyzerCore.h"

class SMPAnalyzerCore : public AnalyzerCore {

public:

  virtual void initializeAnalyzer();
  void FillGenHists(TString pre,TString suf,TLorentzVector genl0,TLorentzVector genl1,TLorentzVector genfsr,double w);
  void FillBasicHists(TString pre,TString suf,const vector<Lepton*>& leps,double weight);
  void FillSystematicHists(TString pre,TString suf,const vector<Lepton*>& leps,map<TString,double> map_systematic);
  void SetupZPtWeight();
  double GetZPtWeight(double zpt,double zrap,Lepton::Flavour flavour);
  double DileptonTrigger_SF(TString SFhistkey0,TString SFhistkey1,const vector<Lepton*>& leps,int sys);
  void PrintGens(const vector<Gen>& gens);
  double GetBinContentUser(TH2* hist,double valx,double valy,int sys);
  void GetGenIndex(const vector<Gen>& gens,int& parton0,int& parton1,int& hardl0,int& hardl1,int& l0,int& l1,vector<int>& photons);
  std::vector<Electron> SMPGetElectrons(TString id, double ptmin, double fetamax);

  static const int zptcor_nptbin=46;
  const double zptcor_ptbin[zptcor_nptbin+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,60,70,80,90,100,120,140,160,180,200,250,400};
  static const int zptcor_nybin=6;
  const double zptcor_ybin[zptcor_nybin+1]={0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH2D *hzpt_muon,*hzpt_electron,*hzpt_norm_muon,*hzpt_norm_electron;
  TString tauprefix;
  double zptcor;

  SMPAnalyzerCore();
  ~SMPAnalyzerCore();

};
#endif

