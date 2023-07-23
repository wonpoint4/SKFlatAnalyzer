#ifndef SkimTree_MuonTnP_h
#define SkimTree_MuonTnP_h

#include "SMPAnalyzerCore.h"

class SkimTree_MuonTnP : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  bool PassSelection(Parameter& p);
  void FillHists(Parameter& p);
  void WriteHist();

  SkimTree_MuonTnP();
  ~SkimTree_MuonTnP();

  TTree* newtree;
  float weight;
  float PUweight,PUweight_up,PUweight_down;
  float prefireweight,prefireweight_up,prefireweight_down;
  float zptweight;
  float z0weight;
  
  bool probe_isTracker;
  bool probe_isGlobal;
  bool probe_isSA;
  bool probe_isSA_unique;
  bool probe_isTight;
  bool probe_isMedium;
  bool probe_isMedium2016a;
  bool probe_TkIsoLoose;
  bool probe_PFIsoTight;
  bool probe_IsoMu24;
  bool probe_IsoMu27;
  bool probe_Mu17Leg1;
  bool probe_Mu8Leg2;

  float probe_pt;
  float probe_pt_cor;
  float probe_eta;
  float probe_phi;
  int probe_q;

  bool tag_IsoMu24;
  bool tag_IsoMu27;
  bool tag_isTight;
  bool tag_isMedium;
  bool tag_isMedium2016a;
  bool tag_TkIsoLoose;
  bool tag_PFIsoTight;

  float tag_pt;
  float tag_pt_cor;
  float tag_eta;
  float tag_phi;
  float tag_q;

  float pair_mass;
  float pair_mass_cor;
  float pair_pt;
  float pair_pt_cor;
  bool pair_EMTF;

  bool pair_gen_matched;
  float pair_gen_mass;
  float probe_gen_pt;
  float probe_gen_eta;
  float probe_gen_phi;
  float probe_gen_dR;
  float probe_gen_reldpt;

};



#endif

