#ifndef SMPAnalyzerCore_h
#define SMPAnalyzerCore_h

#include <tuple>
#include "AnalyzerCore.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TProfile.h"
#include "RoccoR.h"
#include "Aepcor.h"
#include "TH4D.h"
#include "EfficiencyTool.h"
#include "RocPFProb.h"
#include "Weight.h"
#include "ZptCorrection.h"

class SMPAnalyzerCore : public AnalyzerCore {

public:  
  enum{
    NominalWeight=1<<0,
    SystematicWeight=1<<1,
    PDFWeight=1<<2,
    EfficiencyWeight=1<<3,
    LeptonCorrection=1<<4,
  };

  class Variation{
    virtual TString ClassName()=0;
  };
  class VariationWeight: public Variation{
  public:
  VariationWeight(double w):weight(w){}
    double weight=1.;
    virtual TString ClassName(){return "VariationWeight";}
  };
  class VariationMuonMomentum: public Variation{
  public:
  VariationMuonMomentum(int s,int m):set(s),mem(m){}
    int set=0;
    int mem=0;
    virtual TString ClassName(){return "VariationMuonMomentum";}
  };
  class VariationElectronEnergy: public Variation{
  public:
  VariationElectronEnergy(int s,int m):set(s),mem(m){}
    int set=0;
    int mem=0;
    virtual TString ClassName(){return "VariationElectronEnergy";}
  };
  class VariationJES: public Variation{
  public:
  VariationJES(int d):direction(d){}
    int direction=0;
    virtual TString ClassName(){return "VariationJES";}
  };
  class VariationJER: public Variation{
  public:
  VariationJER(int d):direction(d){}
    int direction=0;
    virtual TString ClassName(){return "VariationJER";}
  };
  typedef std::map<TString,std::unique_ptr<SMPAnalyzerCore::Variation>> Variations;
  void AddVariationWeight(Variations& v,TString suffix,double weight);
  void AddVariationMuonMomentum(Variations& v,TString suffix,int set,int mem);
  void AddVariationElectronEnergy(Variations& v,TString suffix,int set,int mem);
  void AddVariationJES(Variations& v,TString suffix,int direction);
  void AddVariationJER(Variations& v,TString suffix,int direction);

  class Parameter{
  public:
    TString channel;
    TString prefix,hprefix,suffix,vsuffix;
    vector<TString> triggers;
    vector<Gen> gens;
    vector<Jet> jets;
    vector<Jet> bjets;
    vector<Muon> muons;
    vector<Electron> electrons;
    vector<Muon> amuons;
    vector<Electron> aelectrons;
    vector<Lepton*> leptons;
    Lepton* lepton0=NULL;
    Lepton* lepton1=NULL;
    Gen truth_lepton0;
    Gen truth_lepton1;
    std::map<TString,double> doublemap;
    std::map<TString,int> intmap;
    int variationbits=NominalWeight;
    Weight default_weight;
    Weight weight;
    int muonmomentum_set=0;
    int default_muonmomentum_set=0;
    int muonmomentum_mem=0;
    int default_muonmomentum_mem=0;
    int electronenergy_set=0;
    int default_electronenergy_set=0;
    int electronenergy_mem=0;
    int default_electronenergy_mem=0;
    int JES_direction=0;
    int default_JES_direction=0;
    int JER_direction=0;
    int default_JER_direction=0;

    TString option;
    struct Key{
      TString electronRECOSF,electronIDSF,electronIDSF2,muonTrackingSF,muonRECOSF,muonIDSF,muonISOSF,DZSF;
      vector<TString> triggerSF;
    };
    struct Weights{
      Weight lumiweight;
      Weight PUweight,PUweight_up,PUweight_down;
      Weight prefireweight,prefireweight_up,prefireweight_down;
      Weight z0weight;
      Weight zptweight,zptweight_g,zptweight_gy,zptweight_gym;
      Weight topptweight;
      Weight weakweight;
      Weight electronRECOSF;
      vector<vector<Weight>> electronRECOSF_sys;
      Weight electronIDSF;
      vector<vector<Weight>> electronIDSF_sys;
      Weight muonTrackingSF;
      vector<vector<Weight>> muonTrackingSF_sys;
      Weight muonRECOSF;
      vector<vector<Weight>> muonRECOSF_sys;
      Weight muonIDSF;
      vector<vector<Weight>> muonIDSF_sys;
      Weight muonISOSF;
      vector<vector<Weight>> muonISOSF_sys;
      Weight triggerSF,triggerSF_up,triggerSF_down,triggerSF_mode1,triggerSF_interpolation;
      vector<vector<Weight>> triggerSF_sys;
      Weight CFSF,CFSF_up,CFSF_down;
      Weight btagSF,btagSF_hup,btagSF_hdown,btagSF_lup,btagSF_ldown,btagSF_hcorr,btagSF_huncorr,btagSF_lcorr,btagSF_luncorr;
      Weight bchargeSF,bchargeSF_s0m0,bchargeSF_s0m1;
    };
    struct Cut{
      double lepton0pt=-1,lepton1pt=-1;
      double muon0pt=-1,muon1pt=-1;
      double electron0pt=-1,electron1pt=-1;
      double amuon0pt=-1,amuon1pt=-1;
      double aelectron0pt=-1,aelectron1pt=-1;
      double jetpt=40.;
      int nelectronmax=-1,nmuonmax=-1;
      int nleptonmin=2;
    };
    Key k;
    Weights w;
    Cut c;

    Parameter();
    ~Parameter();
    Parameter Clone(){ Parameter p=*this; p.SetLeptons(); return p; }
    void SetChannel(TString ch);
    void SetElectronKeys(TString elID,vector<TString> trig);
    void SetElectronKeys(TString elID,TString elID2,vector<TString> trig);
    void SetMuonKeys(TString muID,TString muISO,vector<TString> trig);
    void SetLeptonPtCut(double l0pt,double l1pt);
    void SetLeptons();
    void SetGens(vector<Gen> gs);
    void SetElectrons(vector<Electron> els);
    void SetMuons(vector<Muon> mus);
    void SetAElectrons(vector<Electron> els);
    void SetAMuons(vector<Muon> mus);
  };

  virtual void Apply(Parameter& p,TString vsuf,unique_ptr<Variation>& v);
  virtual void initializeAnalyzer();
  virtual void beginEvent();
  virtual void executeEventWithParameter(Parameter& p);
  virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p.Clone();executeEventWithParameter(pp);}
  virtual void EvalIDSF(Parameter& p);
  virtual void EvalTriggerSF(Parameter& p);
  virtual void EvalDefaultWeight(Parameter& p);
  virtual Variations MakeVariations(const Parameter& p);
  virtual void EvalVariationsPUweight(const Parameter& p,Variations& variations);
  virtual void EvalVariationsPrefireweight(const Parameter& p,Variations& variations);
  virtual void EvalVariationsCF(const Parameter& p,Variations& variations);
  virtual void EvalVariationsBtag(const Parameter& p,Variations& variations);
  virtual void EvalVariationsBcharge(const Parameter& p,Variations& v);
  virtual void EvalVariationsEtc(const Parameter& p,Variations& variations);
  virtual void EvalVariationsEfficiency(const Parameter& p,Variations& variations);
  virtual void EvalVariationsPDF(const Parameter& p,Variations& variations);
  virtual void EvalVariationsMuonMomentum(const Parameter& p,Variations& variations);
  virtual void EvalVariationsElectronEnergy(const Parameter& p,Variations& variations);
  virtual void EvalVariationsJetCorrection(const Parameter& p,Variations& variations);
  virtual bool PassSelection(Parameter& p,bool cutflow=false);
  virtual Parameter MakeParameter(TString channel,TString option="");

  std::map< TString, TH4D* > maphist_TH4D;
  TH4D* GetHist4D(TString histname);
  void FillProfile(TString histname,
		   Double_t value_x, Double_t value_y, Double_t weight,
		   Int_t n_binx, Double_t x_min, Double_t x_max);
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
                Int_t n_binx, const Double_t *xbins,
                Int_t n_biny, const Double_t *ybins,
                Int_t n_binz, const Double_t *zbins,
                Int_t n_binu, const Double_t *ubins);
  void FillHist(TString histname,
                Double_t value_x, Double_t value_y, Double_t value_z, Double_t value_u,
                Double_t weight,
                Int_t n_binx, const Double_t *xbins,
                Int_t n_biny, const Double_t *ybins,
                Int_t n_binz, const Double_t *zbins,
                Int_t n_binu, Double_t u_min, Double_t u_max);
  virtual void WriteHist();

  using AnalyzerCore::FillHist;
  void FillHist(TString histname, double value, map<TString,double> weights, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value, map<TString,double> weights, int n_bin, const double *xbins);
  void FillHist(TString histname,
		double value_x, double value_y,
		map<TString,double> weights,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max);
  void FillHist(TString histname,
		double value_x, double value_y,
		map<TString,double> weights,
		int n_binx, const double *xbins,
		int n_biny, const double *ybins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		map<TString,double> weights,
		int n_binx, double x_min, double x_max,
		int n_biny, double y_min, double y_max,
		int n_binz, double z_min, double z_max);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z,
		map<TString,double> weights,
		int n_binx, const double *xbins,
		int n_biny, const double *ybins,
		int n_binz, const double *zbins);
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
		int n_binx, const double *xbins,
		int n_biny, const double *ybins,
		int n_binz, const double *zbins,
                int n_binu, const double *ubins);
  void FillHist(TString histname,
		double value_x, double value_y, double value_z, double value_u,
		map<TString,double> weights,
		int n_binx, const double *xbins,
		int n_biny, const double *ybins,
		int n_binz, const double *zbins,
                int n_binu, double u_min, double u_max);
  virtual void FillHists(Parameter& p);
  virtual void FillHistsSyst(Parameter p,Variations& v);

  void FillGenHists(TString pre,TString suf,TLorentzVector genl0,TLorentzVector genl1,TLorentzVector genfsr,double w);
  void FillDileptonHists(TString pre,TString suf,Particle* l0,Particle* l1,double w);
  static double GetPtThreshold(TString path);
  static bool IsExists(TString filepath);
  static vector<TString> Split(TString s,TString del);

  void SetupRoccoR();
  double GetZ0Weight(double z0);

  void SetupCFRate();
  double GetCFSF(const Lepton* l,int sys=0);
  double GetCFData(const Lepton* l,int sys=0);
  double GetCFSim(const Lepton* l,int sys=0);
  double GetCFSF(const Parameter& p,int sys=0);
  void DeleteCFRate();
  TH2* hcfrate_data=NULL;
  TH2* hcfrate_mc=NULL;
  TH2* hcfsf=NULL;
  TH2* hcfscale=NULL;

  void SetupMuonTrackingSF();
  double GetMuonTrackingSF(double eta,int sys=0);
  void DeleteMuonTrackingSF();
  TH1* fMuonTrackingSF=NULL;
  bool jSetupMuonTrackingSF=false;

  double GetDYWeakWeight(double mass);

  void SetupFakeRate();
  double GetFakeTF(Parameter& p,TString option="",int sys=0);
  double GetFakeRate(const Lepton *lep);
  double GetFakeRate(Lepton::Flavour flavour,double eta,double pt);
  void DeleteFakeRate();
  TH2* fFakeRate_electron=NULL;
  TH2* fFakeRate_muon=NULL;
  map<TString,TH2*> fFakeTF;

  EfficiencyTool* fEff=NULL;
  void SetupEfficiency();
  void DeleteEfficiency();
  double GetLeptonTriggerSF(TString triggerSF_key,const vector<Lepton*>& leps,int set,int mem,TString option="");
  double GetLeptonTriggerORSF(const Parameter& p,const vector<Lepton*>& leps,int set,int mem,TString option="");
  double GetDileptonTriggerSF(TString SFhistkey0,TString SFhistkey1,TString DZSFhistkey,const vector<Lepton*>& leps,int set,int mem,TString option="");

  void PrintGens(const vector<Gen>& gens);
  static double GetBinContentUser(TH1* hist,double valx,int sys);
  static double GetBinContentUser(TH2* hist,double valx,double valy,int sys);
  static double GetBinContentUser(TH3* hist,double valx,double valy,double valz,int sys);
  void GetAFBLHEParticles(const vector<LHE>& lhes,LHE& p0,LHE& p1,LHE& l0,LHE& l1,LHE& j0);
  void GetAFBGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode);
  static Gen SMPGetGenMatchedLepton(const Lepton& lep, const std::vector<Gen>& gens, int mode=0);
  std::vector<Electron> SMPGetElectrons(TString id, double ptmin, double fetamax);
  std::vector<Muon> SMPGetMuons(TString id,double ptmin,double fetamax);
  void FillCutflow(TString histname,TString label,double weight);
  static TString Replace(TString str,TRegexp reg,TString repl);
  static map<TString,double> SelectWeights(map<TString,double> origin,vector<TString> keys);
  inline map<TString,double> Multiply(map<TString,double> a,double b){
    for(auto& iter:a) iter.second*=b;
    return a;
  }
  
  // Top pt weight
  double GetTopPtReweight2(const std::vector<Gen>& gens);

  // ZptWeight
  void SetupZptWeight();
  void DeleteZptWeight();
  ZptCorrection* fZptCorrection=NULL;

  // L1PrefiringWeight
  virtual void SetupL1PrefiringWeight();
  virtual void DeleteL1PrefiringWeight();
  virtual double getPrefiringRateEcal(double eta, double pt, TH2* h_prefmap, int sys, int mode=0) const;
  virtual double getPrefiringRatePhoton(double eta, double pt, int sys, int mode=0) const;
  virtual double getPrefiringRateJet(double eta, double pt, int sys, int mode=0) const;
  virtual double getPrefiringRateMuon(double eta, double phi, double pt, int sys) const;
  virtual double GetL1PrefiringWeight(int mode=0) const;
  RocPFProb* rocpfprob=NULL;
  TH2* fFGPP=NULL;
  TH2* fFGPM=NULL;
  TH2* fL1Prefiring_photon=NULL;
  TH2* fL1Prefiring_jet=NULL;
  TF1* fL1Prefiring_muon[12]={};
 


  bool IsDYSample=false;
  bool IsTTSample=false;
  bool IsTTLLSample=false;
  Event _event;
  double reductionweight=1;
  vector<LHE> lhes;
  LHE lhe_p0,lhe_p1,lhe_l0,lhe_l1,lhe_j0;
  vector<Gen> gens;
  Gen gen_p0,gen_p1,gen_l0,gen_l1,gen_l0_dressed,gen_l1_dressed,gen_l0_bare,gen_l1_bare;

  RoccoR* roc=NULL;
  Aepcor* rocele=NULL;
  map<TString,TH1*> fRoccorResidual;

  virtual double MuonMomentumCorrection(const Muon& muon,int set=0,int member=0);
  virtual std::vector<Muon> MuonMomentumCorrection(const vector<Muon>& muons,int set=0,int member=0,bool sort=true);
  virtual double ElectronEnergyCorrection(const Electron& electron,int set=0,int member=0);
  virtual std::vector<Electron> ElectronEnergyCorrection(const vector<Electron>& electrons,int set=0,int member=0,bool sort=true);

  double GetPFMET_T1Smear() const;
  TString GetSkimName() const;

  bool PassSLT1(const Lepton* lep) const;
  bool PassSLT2(const Lepton* lep) const;
  bool PassDLT1(const Lepton* lep) const;
  bool PassDLT2(const Lepton* lep) const;

  static vector<vector<Weight>> Make2DWeights(const vector<int>& structure);

  virtual double GetBchargeSF(const Jet& bjet,int set=-1,int mem=-1) const;
  virtual double GetBchargeSF(const Parameter& p,int set=-1,int mem=-1) const;

  SMPAnalyzerCore();
  ~SMPAnalyzerCore();

};
#endif

