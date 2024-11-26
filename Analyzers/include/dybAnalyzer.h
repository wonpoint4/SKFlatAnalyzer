#ifndef dybAnalyzer_h
#define dybAnalyzer_h

#include "TKey.h"
#include "SMPAnalyzerCore.h"
#include "AFBAnalyzer.h"
#include "LHAPDF/Reweighting.h"

class dybAnalyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEventWithParameter(TString channel);
  virtual void executeEventGen();
  virtual void executeEvent();
  virtual bool IsFiredTriggers(TString channel);
  virtual bool HasDileptons(TString channel);
  virtual double jetCharge(const Jet& jet);
  virtual double GetBTaggingReweight_1a_2WP(const vector<Jet>& jets, JetTagging::Parameters jtpT, JetTagging::Parameters jtpL, string Syst);
  double GetbChargeSFWeight(const vector<Jet>& jets, int mode, int sys);

  // PUJetID, SF
  void SetupPUJetWeight();
  TH2F *heff_data=NULL;
  TH2F *hmistag_data=NULL;
  TH2F *heff_mc=NULL;
  TH2F *hmistag_mc=NULL;
  bool PUJetIDPass(Jet jet, TString ID);
  double GetPUJetWeight(const vector<Jet>& jets, TString ID, int sys);

  dybAnalyzer();
  ~dybAnalyzer();

  virtual double GetCosThetaCS(const Particle *p0, const Particle *p1, int direction=0);
  virtual double GetCosThetaRecoil(const Particle *p0, const Particle *p1, Particle *b, int mode=0);

  TString prefix, gprefix, hprefix, suffix;
  vector<Gen> gens;
  vector<Jet> jets;
  vector<Jet> bjets;
  vector<Lepton*> leptons;
  vector<Muon> muons;
  vector<Electron> electrons;
  Lepton* lepton0 = NULL;
  Lepton* lepton1 = NULL;
  Jet* jet0 = NULL;
  double bcharge = 0;
  std::map<TString,double> map_weight;

  double leptonTrackingSF =1.;
  double leptonRECOSF =1.;
  double leptonIDSF =1.;
  double leptonTriggerSF =1.;

  double lumiweight = 1.;
  double PUweight = 1.;
  double prefireweight = 1.;
  double zptweight =1.;
  double weakweight = 1.;
  double btagSF = 1.;
  double topptweight = 1.;
  double pujetSF = 1.;

  bool IsNominalRun = true;
  bool IsSkimmed = false;

  LHAPDF::PDF* PDFbase = NULL;
  LHAPDF::PDF* PDFnf4 = NULL;

  static const int afb_mbinnum = 40;
  static constexpr const double afb_mbin[afb_mbinnum+1] = {52,56,60,65,70,74,77,80,82,84,86,88,89,90,91,92,93,94,96,98,100,103,106,110,115,120,130,140,150,175,200,240,280,340,400,500,600,700,800,1000,3000};
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

