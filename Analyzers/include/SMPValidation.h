#ifndef SMPValidation_h
#define SMPValidation_h

#include "SMPAnalyzerCore.h"

class SMPValidation : public SMPAnalyzerCore {

public:

  void executeEventFromParameter(TString channelname,Event* ev);
  void executeEvent();

  SMPValidation();
  ~SMPValidation();

};



#endif

