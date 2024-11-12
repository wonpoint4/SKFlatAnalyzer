#ifndef GenTest_h
#define GenTest_h

#include "AFBAnalyzerSyst.h"

class GenTest : public AFBAnalyzerSyst {

public:
  virtual void executeEvent();
  virtual void GetAFBGenParticles(const vector<Gen>& gens,Gen& parton0,Gen& parton1,Gen& l0,Gen& l1,int mode);
  GenTest();
  ~GenTest();

};



#endif

