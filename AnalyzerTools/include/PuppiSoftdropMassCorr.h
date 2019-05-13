#ifndef PuppiSoftdropMassCorr_h
#define PuppiSoftdropMassCorr_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TString.h"
#include "TF1.h"

#include "FatJet.h"

#include "TDirectoryHelper.h"

class PuppiSoftdropMassCorr{

public:

  PuppiSoftdropMassCorr();
  ~PuppiSoftdropMassCorr();

  TDirectory *histDir;
  void ReadHistograms();

  int DataYear;
  void SetDataYear(int i);

  TF1 *func_gen;
  TF1 *func_reco_0eta1v3;
  TF1 *func_reco_1v3eta2v5;

  std::vector<FatJet> Correct(std::vector<FatJet> jets);

};

#endif
