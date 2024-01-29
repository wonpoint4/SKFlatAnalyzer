#ifndef Chi2Analyzer_h
#define Chi2Analyzer_h

#include "SMPAnalyzerCore.h"

class Chi2Analyzer : public SMPAnalyzerCore {

public:

  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual Parameter MakeParameter(TString key,TString option="");
  virtual void EvalDefaultWeight(Parameter& p);
  virtual Variations MakeVariations(const Parameter& p);
  virtual void FillHists(Parameter& p);

  Chi2Analyzer();
  ~Chi2Analyzer();

  static const int rochester_nmbin=11;
  const double rochester_mbins[rochester_nmbin+1]={54,66,76,82,86,89.5,92.7,96,100,106,116,150};
  static const int rochester_nybin=6;
  const double rochester_ybins[rochester_nybin+1]={0,0.4,0.8,1.2,1.6,2.0,2.4};

  static const int netabin_muonID=48;
  const double etabins_muonID[netabin_muonID+1]={-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
  static const int nptbin_muonID=10;
  const double ptbins_muonID[nptbin_muonID+1]={10, 15, 20, 25, 30, 35, 40, 45, 50, 70, 200};

  static const int netabin_electronID=50;
  const double etabins_electronID[netabin_electronID+1]={-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  static const int nptbin_electronID=12;
  const double ptbins_electronID[nptbin_electronID+1]={10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 500};
};



#endif

