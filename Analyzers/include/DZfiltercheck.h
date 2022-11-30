#ifndef DZfiltercheck_h
#define DZfiltercheck_h

#include "SMPAnalyzerCore.h"

class DZfiltercheck : public SMPAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  DZfiltercheck();
  ~DZfiltercheck();

  vector<Muon> muons;

  vector<double> vec_etabins = {-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
  vector<double> vec_ptbins = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 200.};

  const int NEtaBin = vec_etabins.size()-1;
  const int NPtBin = vec_ptbins.size()-1;
};



#endif

