#ifndef SMPAnalyzerCore_h
#define SMPAnalyzerCore_h

#include <tuple>
#include "AnalyzerCore.h"
#include "TRegexp.h"
#include "RoccoR.h"
#include "RocelecoR.h"
#include "TH4D.h"

class SMPAnalyzerCore : public AnalyzerCore {

public:  
  virtual void initializeAnalyzer();
  
  std::map< TString, TH4D* > maphist_TH4D;
  TH4D* GetHist4D(TString histname);
  void FillHist(TString histname,
                Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                Double_t weight,
                Int_t n_binx, Double_t x_min, Double_t x_max,
                Int_t n_biny, Double_t y_min, Double_t y_max,
                Int_t n_binz, Double_t z_min, Double_t z_max,
                Int_t n_binu, Double_t u_min, Double_t u_max);
  void FillHist(TString histname,
                Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                Double_t weight,
                Int_t n_binx, Double_t *xbins,
                Int_t n_biny, Double_t *ybins,
                Int_t n_binz, Double_t *zbins,
                Int_t n_binu, Double_t *ubins);
  void FillHist(TString histname,
                Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                Double_t weight,
                Int_t n_binx, Double_t *xbins,
                Int_t n_biny, Double_t *ybins,
                Int_t n_binz, Double_t *zbins,
                Int_t n_binu, Double_t u_min, Double_t u_max);
  virtual void WriteHist();

  using AnalyzerCore::FillHist;
  void FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double *xbins);
  void FillHist(TString histname,
		double value_x, double value_y,
		map<TString,double> weights,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max);
  void FillHist(TString histname,
		double value_x, double value_y,
		map<TString,double> weights,
		int n_binx, double *xbins,
		int n_biny, double *ybins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		map<TString,double> weights,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max,
		int n_binz, double z_min, double z_max);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		map<TString,double> weights,
		int n_binx, double *xbins,
		int n_biny, double *ybins,
		int n_binz, double *zbins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z, double value_u,
		map<TString,double> weights,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max,
		int n_binz, double z_min, double z_max,
                int n_binu, double u_min, double u_max);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z, double value_u,
		map<TString,double> weights,
		int n_binx, double *xbins,
		int n_biny, double *ybins,
		int n_binz, double *zbins,
                int n_binu, double *ubins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z, double value_u,
		map<TString,double> weights,
		int n_binx, double *xbins,
		int n_biny, double *ybins,
		int n_binz, double *zbins,
                int n_binu, double u_min, double u_max);

  void FillGenHists(TString pre,TString suf,TLorentzVector genl0,TLorentzVector genl1,TLorentzVector genfsr,double w);
  void FillDileptonHists(TString pre,TString suf,Particle* l0,Particle* l1,double w);
  void SetupZptWeight();
  void SetupZ0Weight();
  void SetupRoccoR();
  double GetZptWeight(double zpt,double zrap,Lepton::Flavour flavour);
  double GetZ0Weight(double z0);
  double Lepton_SF(TString histkey,const Lepton* lep,int sys);
  double LeptonTrigger_SF(TString triggerSF_key,const vector<Lepton*>& leps,int sys);
  double DileptonTrigger_SF(TString SFhistkey0,TString SFhistkey1,const vector<Lepton*>& leps,int sys);
  void PrintGens(const vector<Gen>& gens);
  double GetBinContentUser(TH2* hist,double valx,double valy,int sys);
  void GetDYLHEParticles(const vector<LHE>& lhes,LHE& l0,LHE& l1);
  void GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,bool dressed);
  std::vector<Electron> SMPGetElectrons(TString id, double ptmin, double fetamax);
  std::vector<Muon> SMPGetMuons(TString id,double ptmin,double fetamax);
  void FillCutflow(TString histname,TString label,double weight);
  TString Replace(TString str,TRegexp reg,TString repl);
  inline map<TString,double> Multiply(map<TString,double> a,double b){
    for(auto& [name,value]:a) value*=b;
    return a;
  }
  
  static const int zptcor_nptbin=46;
  const double zptcor_ptbin[zptcor_nptbin+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,60,70,80,90,100,120,140,160,180,200,250,400};
  static const int zptcor_nybin=6;
  const double zptcor_ybin[zptcor_nybin+1]={0,0.4,0.8,1.2,1.6,2.0,2.4};

  map<TString,TH2D*> map_hist_zpt;
  TH1D *hz0;
  TString tauprefix;
  double zptcor;
  bool IsDYSample=false;

  RoccoR* roc;
  RocelecoR* rocele;
  std::vector<Muon> MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set=0,int member=0);
  std::vector<Electron> ElectronEnergyCorrection(const vector<Electron>& electrons,int set=0,int member=0);

  SMPAnalyzerCore();
  ~SMPAnalyzerCore();

};
#endif

