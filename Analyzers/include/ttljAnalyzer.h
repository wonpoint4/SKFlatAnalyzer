#ifndef ttljAnalyzer_h
#define ttljAnalyzer_h

#include "SMPAnalyzerCore.h"
#include "dybAnalyzer.h"

class ttljAnalyzer : public dybAnalyzer {

 public:

  virtual void initializeAnalyzer();
  virtual void GetTTLJGenParticles(const vector<Gen>& gens, Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,Gen& b0,Gen& b1,Gen& j0,Gen& j1, int mode);
  virtual void executeEventWithParameter(TString channel);
  virtual void executeEvent();
  virtual void executeEventGen();
  virtual bool Hasleptons(TString channel);
  void FillingLikelihood(TString channel, vector<Jet> bjets, vector<Jet> ajets, double weight, unsigned int mode1=0, unsigned int mode2=0);
  void SetupLikelihoods(unsigned int mode1=0, unsigned int mode2=0);
  vector<unsigned int> Finding_bbjj_byLikelihood(TString channel, vector<Jet> bjets, vector<Jet> ajets, vector<double>& Likelihood_ratios);
  bool Kinematic_Cut(vector<Jet> bjets, vector<Jet> ajets);

  Gen gen_parton0,gen_parton1, gen_b0,gen_b1, gen_l0,gen_l1, gen_j0,gen_j1;
  Gen* gen_lepb = NULL;
  Gen* gen_hadb = NULL;
  Lepton* lepton0 = NULL;
  Jet* bjet0 = NULL;
  Jet* bjet1 = NULL;
  Jet* ajet0 = NULL;
  Jet* ajet1 = NULL;
  double bcharge0 = 0;
  double bcharge1 = 0;
  double acharge0 = 0;
  double acharge1 = 0;
  TLorentzVector met;

  // Likelihoods
  TH1* hMbl_correct_E = NULL;
  TH1* hMbl_wrong_E = NULL;
  TH1* hMblMET_correct_E = NULL;
  TH1* hMblMET_wrong_E = NULL;
  TH1* hMbjj_correct_E = NULL;
  TH1* hMbjj_wrong_E = NULL;
  TH1* hMjj_correct_E = NULL;
  TH1* hMjj_wrong_E = NULL;
  TH1* hdRtt_correct_E = NULL;
  TH1* hdRtt_wrong_E = NULL;
  TH1* hdPhitt_correct_E = NULL;
  TH1* hdPhitt_wrong_E = NULL;

  TH1* hMbl_correct_m = NULL;
  TH1* hMbl_wrong_m = NULL;
  TH1* hMblMET_correct_m = NULL;
  TH1* hMblMET_wrong_m = NULL;
  TH1* hMbjj_correct_m = NULL;
  TH1* hMbjj_wrong_m = NULL;
  TH1* hMjj_correct_m = NULL;
  TH1* hMjj_wrong_m = NULL;
  TH1* hdRtt_correct_m = NULL;
  TH1* hdRtt_wrong_m = NULL;
  TH1* hdPhitt_correct_m = NULL;
  TH1* hdPhitt_wrong_m = NULL;

  ttljAnalyzer();
  ~ttljAnalyzer();

 private:

};

#endif
