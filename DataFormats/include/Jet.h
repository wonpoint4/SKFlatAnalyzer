#ifndef Jet_h
#define Jet_h

#include "Particle.h"

class Jet: public Particle {
public:

  Jet();
  ~Jet();

  void SetArea(double area);
  void SetGenFlavours(double pf, double hf);
  inline int partonFlavour() const { return j_partonFlavour; }
  inline int hadronFlavour() const { return j_hadronFlavour; }
  void SetTaggerResults(std::vector<double> ds);
  void SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE);
  void SetMultiplicities(double cM, double nM);
  void SetPileupJetId(double v);

  void SetEnShift(double en_up, double en_down);
  inline double EnShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_En_up;
    else return j_En_down;
  }
  void SetResShift(double res_up, double res_down);
  inline double ResShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_Res_up;
    else return j_Res_down;
  }

  void SetTightJetID(double b);
  void SetTightLepVetoJetID(double b);
  inline bool Pass_tightJetID() const { return j_tightJetID; }
  inline bool Pass_tightLepVetoJetID() const { return j_tightLepVetoJetID; }

  bool PassID(TString ID) const;

  enum Tagger{
    CSVv2,
    DeepCSV,DeepJet, DeepCvsL, DeepCvsB,
    CvsL, CvsB,
    DeepFlavour_b, DeepFlavour_bb, DeepFlavour_lepb, DeepFlavour_c, DeepFlavour_uds, DeepFlavour_g,
  };
  
  inline TString TaggerString(Jet::Tagger t) const {
    if(t == Jet::CSVv2) return "CSVv2"; 
    if(t == Jet::DeepCSV) return "DeepCSV";    
    if(t == Jet::DeepJet) return "DeepJet";    
    if(t == Jet::DeepCvsL) return "DeepCvsL";    
    if(t == Jet::CvsL) return "CvsL";    
    if(t == Jet::CvsB) return "CvsB";    
    if(t == Jet::DeepFlavour_b) return "DeepFlavour_b";    
    if(t == Jet::DeepFlavour_bb) return "DeepFlavour_bb";    
    if(t == Jet::DeepFlavour_lepb) return "DeepFlavour_lepb";    
    if(t == Jet::DeepFlavour_c) return "DeepFlavour_c";    
    if(t == Jet::DeepFlavour_uds) return "DeepFlavour_uds";    
    if(t == Jet::DeepFlavour_g) return "DeepFlavour_g";    
    
    return "None";
  }

  enum WP{
    Loose, Medium, Tight
  };
  
  inline TString WPString(Jet::WP wp) const {
    if(wp == Jet::Loose) return "Loose";
    if(wp == Jet::Medium) return "Medium";
    if(wp == Jet::Tight) return "Tight";
    return "None";
  }

  double GetTaggerResult(Tagger tg);

private:

  double  j_area;
  int j_partonFlavour;
  int j_hadronFlavour;
  double  j_CSVv2;
  double  j_DeepCSV;
  double  j_DeepFlavour_b;
  double  j_DeepFlavour_bb;
  double  j_DeepFlavour_lepb;
  double  j_DeepFlavour_c;
  double  j_DeepFlavour_uds;
  double  j_DeepFlavour_g;
  double  j_CvsL;
  double  j_CvsB;
  double  j_DeepCvsL;
  double  j_DeepCvsB;
  double  j_chargedHadronEnergyFraction;
  double  j_neutralHadronEnergyFraction;
  double  j_neutralEmEnergyFraction;
  double  j_chargedEmEnergyFraction;
  double j_muonEnergyFraction;
  int j_chargedMultiplicity;
  int j_neutralMultiplicity;
  double j_PileupJetId;
  double j_En_up;
  double j_En_down;
  double j_Res_up;
  double j_Res_down;
  bool j_tightJetID, j_tightLepVetoJetID;

  ClassDef(Jet,1)
};

#endif
