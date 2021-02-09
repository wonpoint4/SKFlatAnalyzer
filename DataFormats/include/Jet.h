#ifndef Jet_h
#define Jet_h

#include "Particle.h"
#include "JetTaggingParameters.h"

class Jet: public Particle {
public:

  Jet();
  ~Jet();

  void SetArea(double area);
  void SetGenFlavours(double pf, double hf, double pdg);
  inline int partonFlavour() const { return j_partonFlavour; }
  inline int hadronFlavour() const { return j_hadronFlavour; }
  inline int partonPdgId() const { return j_partonPdgId; }
  inline double PileupJetId() const { return j_PileupJetId; }
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

  double GetTaggerResult(JetTagging::Tagger tg) const;

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
  int j_partonPdgId;
  double j_PileupJetId;
  double j_En_up;
  double j_En_down;
  double j_Res_up;
  double j_Res_down;
  bool j_tightJetID, j_tightLepVetoJetID;

  ClassDef(Jet,1)
};

#endif
