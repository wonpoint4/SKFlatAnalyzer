#ifndef DZAnalyzer_h
#define DZAnalyzer_h

#include "SMPAnalyzerCore.h"

class DZAnalyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual Variations MakeVariations(const Parameter& p);
  virtual void FillHists(Parameter& p);
  virtual void FillHistsDZ(Parameter& p,TString suffix);
  virtual bool PassSelection(Parameter& p);

  DZAnalyzer();
  ~DZAnalyzer();

  virtual void SetupDZSF();
  virtual double GetDZSF(Lepton* lep);
  virtual double GetDZSF_DZ(Lepton* lep1,Lepton* lep2);
  TH2 *fDZSF_muon=NULL;
  TH2 *fDZSF_electron=NULL;
  TH1 *fDZSF_muon_DZ=NULL;
  TH1 *fDZSF_electron_DZ=NULL;

  static const int ptbinnum=2;
  const double ptbin[ptbinnum+1]={10,25,200};
  static const int etabinnum=4;
  const double etabin[etabinnum+1]={0,1.3,1.5,2.1,2.5};

  static const int fineptbinnum=9;
  const double fineptbin[fineptbinnum+1]={10,15,20,25,30,35,40,50,60,200};
  static const int fineetabinnum=50;
  const double fineetabin[fineetabinnum+1]={-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
};



#endif

