#ifndef testMuon_h
#define testMuon_h

#include "AnalyzerCore.h"

class testMuon : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  testMuon();
  ~testMuon();

};



#endif

