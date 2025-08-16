#ifndef dybAnalyzer_h
#define dybAnalyzer_h

#include "TKey.h"
#include "SMPAnalyzerCore.h"
#include "AFBAnalyzer.h"
#include "LHAPDF/Reweighting.h"

class dybAnalyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEventWithParameter(TString channel, TString option="");
  virtual void executeEventGen();
  virtual void executeEvent();
  virtual bool IsFiredTriggers(TString channel);
  virtual bool HasDileptons(TString channel);
  virtual double jetCharge(const Jet& jet);
  virtual double GetBTaggingReweight_1a_2WP(const vector<Jet>& jets, JetTagging::Parameters jtpT, JetTagging::Parameters jtpL, string Syst="central");
  double GetbChargeSFWeight(const vector<Jet>& jets, unsigned int mode, int sys=0, TString bChargeBins="012345");

  // PUJetID, SF
  void SetupPUJetWeight();
  TH2F *heff_sf = NULL;
  TH2F *heff_sf_unc = NULL;
  bool PUJetIDPass(Jet jet, TString ID);
  double GetPUJetWeight(const vector<Jet>& jets, TString ID, int sys);

  dybAnalyzer();
  ~dybAnalyzer();

  virtual double GetCosThetaCS(const Particle *p0, const Particle *p1, int direction=0);
  virtual double GetCosThetaRecoil(const Particle *p0, const Particle *p1, Particle *b, int mode=0);
  vector<vector<double>> Make2DWeights(const vector<int>& structure);

  TString prefix, gprefix, hprefix, suffix;
  vector<Gen> gens;
  vector<Jet> jets_raw;
  vector<Jet> jets;
  vector<Jet> bjets;
  vector<Lepton*> leptons;
  vector<Muon> muons;
  vector<Muon> muons_raw;
  vector<Electron> electrons;
  vector<Electron> electrons_raw;
  Lepton* lepton0 = NULL;
  Lepton* lepton1 = NULL;
  Jet* jet0 = NULL;
  double bcharge = 0;
  std::map<TString,double> map_weight;

  double muonTrackingSF =1.;
  double muonRECOSF =1.;
  double muonIDSF =1.;
  double muonTriggerSF =1.;
  double electronRECOSF =1.;
  double electronIDSF =1.;
  double electronTriggerSF =1.;
  vector<vector<double>> muonTrackingSF_sys;
  vector<vector<double>> muonRECOSF_sys;
  vector<vector<double>> muonIDSF_sys;
  vector<vector<double>> muonTriggerSF_sys;
  vector<vector<double>> electronRECOSF_sys;
  vector<vector<double>> electronIDSF_sys;
  vector<vector<double>> electronTriggerSF_sys;

  double lumiweight = 1.;
  double PUweight = 1.;
  double prefireweight = 1.;
  double zptweight = 1., zptweight_gym = 1.;
  double weakweight = 1.;
  double btagSF = 1.;
  double topptweight = 1.;
  double pujetSF = 1.;
  double bchargeSF = 1.;

  bool IsNominalRun = true;
  bool IsSkimmed = false;
  bool IsNominalLike = true;

  LHAPDF::PDF* PDFbase = NULL;
  LHAPDF::PDF* PDFnf4 = NULL;

  // Old binnings or HS's binning
  static const int afb_mbinnum_original = 40;
  static constexpr const double afb_mbin_original[afb_mbinnum_original+1] = {52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};
  static const int afb_mbinnum_original2 = 26;
  static constexpr const double afb_mbin_original2[afb_mbinnum_original2+1] = {52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,200};
  static const int afb_mbinnum_HS = 8;
  static constexpr const double afb_mbin_HS[afb_mbinnum_HS+1] = {52,65,77,106,140,200,280,400,3000};

  static const int afb_mbinnum = 12;
  static constexpr const double afb_mbin[afb_mbinnum+1] = {52,66,76,82,86,89.5,92.7,96,100,106,116,150,200};
  static const int afb_chbinnum = 6;
  static constexpr const double afb_chbin[afb_chbinnum+1] = {0., 0.1, 0.2, 0.6, 1.0, 3.0, 5.0};
  //static const int afb_chbinnum = 12;
  //static constexpr const double afb_chbin[afb_chbinnum+1] = {-5.0, -3.0, -1.0, -0.6, -0.2, -0.1, 0., 0.1, 0.2, 0.6, 1.0, 3.0, 5.0};
  static const int afb_ybinnum = 12;
  static constexpr const double afb_ybin[afb_ybinnum+1] = {-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
  static const int afb_ptbinnum = 30;
  static constexpr const double afb_ptbin[afb_ptbinnum+1] = {0,2,4,6,8,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,120,140,160,190,250,400,650};
};

#endif

