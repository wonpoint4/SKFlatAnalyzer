#ifndef MCCorrection_h
#define MCCorrection_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "Lepton.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "Jet.h"
#include "Gen.h"

#include "TDirectoryHelper.h"
#include "GenFinderForDY.h"

class MCCorrection{

public:

  MCCorrection();
  ~MCCorrection();

  TDirectory *histDir;
  void ReadHistograms();

  TString MCSample;
  void SetMCSample(TString s);

  int DataYear;
  void SetDataYear(int i);

  bool IsFastSim;
  void SetIsFastSim(bool b);

  bool IgnoreNoHist;

  double MuonReco_SF(TString key, double eta, double p, int sys=0);
  double MuonID_SF(TString ID, double eta, double pt, int sys=0);
  double MuonISO_SF(TString ID, double eta, double pt, int sys=0);
  double MuonTrigger_SF(TString ID, double eta, double pt, int sys=0);
  double MuonTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon>& muons, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon *>& muons, int sys=0);

  std::map< TString, TH2F* > map_hist_Muon;

  double ElectronReco_SF(double sceta, double pt, int sys=0);
  double ElectronID_SF(TString ID, double sceta, double pt, int sys=0);
  double ElectronTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron>& electrons, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron *>& electrons, int sys=0);
  std::map< TString, TH2F* > map_hist_Electron;
  std::map< TString, TGraphAsymmErrors* > map_graph_Electron;

  std::map< TString, TH2F* > map_hist_prefire;
  double GetPrefireWeight(const std::vector<Photon>& photons, const std::vector<Jet>& jets, int sys);

  std::map< TString, TH1D* > map_hist_pileup;
  double GetPileUpWeightBySampleName(int N_pileup, int syst);
  double GetPileUpWeight(int N_pileup, int syst);

  double GetTopPtReweight(const std::vector<Gen>& gens);

  TH2D *hist_DYPtReweight_2D;
  GenFinderForDY *genFinderDY;
  double GetOfficialDYReweight(const vector<Gen>& gens, int sys=0);
};

#endif
