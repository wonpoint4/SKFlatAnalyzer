#ifndef MeasureJetTaggingEfficiency_h
#define MeasureJetTaggingEfficiency_h

#include "AnalyzerCore.h"


class MeasureJetTaggingEfficiency : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEvent();

  MeasureJetTaggingEfficiency();
  ~MeasureJetTaggingEfficiency();

  //==== Read what to measrue from data/Run2Legacy_v4/<Year>/BTag/taggermap.txt
  vector<string> Taggers;
  vector<string> WPs;
  vector<double> CutValues;

};



#endif

