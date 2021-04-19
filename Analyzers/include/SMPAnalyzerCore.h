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
  enum{
    NominalWeight=1<<0,
    SystematicWeight=1<<1,
    PDFWeight=1<<2,
  };

  class Parameter{
  public:
    TString channel;
    TString prefix,hprefix,suffix;
    vector<TString> triggers;
    vector<Gen> gens;
    vector<Muon> muons;
    vector<Electron> electrons;
    vector<Muon> amuons;
    vector<Electron> aelectrons;
    vector<Lepton*> leptons;
    Lepton* lepton0=NULL;
    Lepton* lepton1=NULL;
    Gen truth_lepton0;
    Gen truth_lepton1;
    std::map<TString,double> weightmap;
    std::map<TString,double> doublemap;
    std::map<TString,int> intmap;
    int weightbit=NominalWeight;
    TString option;
    struct Key{
      TString electronIDSF,muonIDSF,muonISOSF;
      vector<TString> triggerSF;
    };
    struct Weight{
      double lumiweight=1;
      double PUweight=1,PUweight_up=1,PUweight_down=1;
      double prefireweight=1,prefireweight_up=1,prefireweight_down=1;
      double z0weight=1;
      double zptweight=1;
      double RECOSF=1,RECOSF_up=1,RECOSF_down=1;
      double IDSF=1,IDSF_up=1,IDSF_down=1;
      double ISOSF=1,ISOSF_up=1,ISOSF_down=1;
      double triggerSF=1,triggerSF_up=1,triggerSF_down=1;
      double CFSF=1,CFSF_up=1,CFSF_down=1;
    };
    struct Cut{
      double lepton0pt=-1,lepton1pt=-1;
      double muon0pt=-1,muon1pt=-1;
      double electron0pt=-1,electron1pt=-1;
      double amuon0pt=-1,amuon1pt=-1;
      double aelectron0pt=-1,aelectron1pt=-1;
    };
    Key k;
    Weight w;
    Cut c;

    Parameter();
    ~Parameter();
    void SetChannel(TString ch);
    void SetElectronKeys(TString elID,vector<TString> trig);
    void SetMuonKeys(TString muID,TString muISO,vector<TString> trig);
    void SetLeptonPtCut(double l0pt,double l1pt);
    void SetLeptons();
    void SetGens(vector<Gen> gs);
    void SetElectrons(vector<Electron> els);
    void SetMuons(vector<Muon> mus);
    void SetAElectrons(vector<Electron> els);
    void SetAMuons(vector<Muon> mus);
  };

  virtual void initializeAnalyzer();
  virtual void beginEvent();
  virtual void executeEventWithParameter(Parameter p);
  virtual void EvalIDSF(Parameter& p);
  virtual void EvalTriggerSF(Parameter& p);
  virtual bool PassSelection(Parameter& p);
  virtual Parameter MakeParameter(TString channel,TString option="");

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
  virtual void FillHists(Parameter& p);

  void FillGenHists(TString pre,TString suf,TLorentzVector genl0,TLorentzVector genl1,TLorentzVector genfsr,double w);
  void FillDileptonHists(TString pre,TString suf,Particle* l0,Particle* l1,double w);
  static double GetPtThreshold(TString path);
  static bool IsExists(TString filepath);
  void SetupZptWeight();
  void SetupZ0Weight();
  void SetupRoccoR();
  double GetBTaggingReweight_1a_2WP(const vector<Jet>& jets, JetTagging::Parameters jtpT, JetTagging::Parameters jtpL, string Syst);
  void SetupPUJetWeight(TString ID="Medium");
  double GetPUJetWeight(const vector<Jet>& jets, int sys);
  bool isGenMatchedJet(const Jet& jet, const vector<Gen>& gens);
  double GetZptWeight(double zpt,double zrap,Lepton::Flavour flavour);
  double GetZ0Weight(double z0);

  void SetupCFRate();
  double GetCFSF(const Lepton* l,int sys=0);
  double GetCFSF(const Parameter& p,int sys=0);
  void DeleteCFRate();
  TH2* hcfrate_data=NULL;
  TH2* hcfrate_mc=NULL;
  TH2* hcfsf=NULL;

  double Lepton_SF(TString histkey,const Lepton* lep,int sys);
  double LeptonTrigger_SF(TString triggerSF_key,const vector<Lepton*>& leps,int sys);
  double LeptonTriggerOR_SF(TString triggerSF_key0,TString triggerSF_key1,const vector<Lepton*>& leps,int sys);
  double DileptonTrigger_SF(TString SFhistkey0,TString SFhistkey1,const vector<Lepton*>& leps,int sys);
  void PrintGens(const vector<Gen>& gens);
  void PrintLHEs(vector<LHE>& lhes);
  double GetBinContentUser(TH2* hist,double valx,double valy,int sys);
  double GetBinContentUser(TH3* hist,double valx,double valy,double valz,int sys);
  void GetDYLHEParticles(const vector<LHE>& lhes,LHE& l0,LHE& l1);
  void GetDYLHEParticles(const vector<LHE>& lhes,LHE& l0,LHE& l1,LHE& j0);
  void GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode);
  void GetDYGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& j0,int mode);
  static Gen SMPGetGenMatchedLepton(const Lepton& lep, const std::vector<Gen>& gens, int mode=0);
  std::vector<Electron> SMPGetElectrons(TString id, double ptmin, double fetamax);
  std::vector<Muon> SMPGetMuons(TString id,double ptmin,double fetamax);
  void FillCutflow(TString histname,TString label,double weight);
  TString Replace(TString str,TRegexp reg,TString repl);
  inline map<TString,double> Multiply(map<TString,double> a,double b){
    for(auto& iter:a) iter.second*=b;
    return a;
  }
  
  static const int zptcor_nptbin=46;
  const double zptcor_ptbin[zptcor_nptbin+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,60,70,80,90,100,120,140,160,180,200,250,400};
  static const int zptcor_nybin=6;
  const double zptcor_ybin[zptcor_nybin+1]={0,0.4,0.8,1.2,1.6,2.0,2.4};

  map<TString,TH2D*> map_hist_zpt;

  TF1 *hz0_data=NULL, *hz0_mc=NULL;
  TH2F *heff_data=NULL, *hmistag_data=NULL, *heff_mc=NULL, *hmistag_mc=NULL;
  bool IsDYSample=false;
  Event _event;
  double reductionweight=1;
  double pujetweight=1;
  double btagweight=1;

  vector<LHE> lhes;
  LHE lhe_l0,lhe_l1,lhe_j0;
  vector<Gen> gens;
  Gen gen_p0,gen_p1,gen_l0,gen_l1,gen_j0,gen_l0_dressed,gen_l1_dressed,gen_l0_bare,gen_l1_bare;

  RoccoR* roc=NULL;
  RocelecoR* rocele=NULL;

  std::vector<Muon> MuonMomentumCorrection(const vector<Muon>& muons,int sys,int set=0,int member=0);
  std::vector<Electron> ElectronEnergyCorrection(const vector<Electron>& electrons,int set=0,int member=0);

  SMPAnalyzerCore();
  ~SMPAnalyzerCore();

};
#endif

