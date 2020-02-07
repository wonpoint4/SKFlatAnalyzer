#include "BTagCalibrationStandalone.h"
#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include <sstream>

using namespace std;

class BTagSFUtil{

 public:
    
  // constructors
  BTagSFUtil(string MeasurementType, string BTagAlgorithmBF, TString OperatingPoint, int year, bool period,  int SystematicIndex);

  // destructor
  ~BTagSFUtil();


  // Function used in analyzer to determine if jet is tagged
  bool IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int seed);
  bool IsUncorrectedTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta);


  // Function to setup rand numb generator
  void RandInit(int seed);


  // Functions used in IsTagged to determine tag rate with syst                                                                                                                                                                                                      
  float GetJetSF(int JetFlavor, float JetPt, float JetEta);
  float GetJetSFPeriodDependant(int JetFlavor, float JetPt, float JetEta, TString iperiod, bool isperiodDep, TString tagger);
  float JetTagEfficiency(int JetFlavor, float JetPt, float JetEta);


  // Function to warn user of incorrect setup
  void ErrorLoadingCSV(TString tag_key, bool isperiodDep, TString TaggerName);

  // New functions for SKAnalzyer
  void SetMCSample(TString s);
  void SetDataYear(int i);
  void SetPeriodDependancy(bool b);


 private:


  // map to store BTagCalibrationReader objects, determined by setup/year
  std::map <TString, BTagCalibrationReader*> ReaderMap;
  
  void GetBTagPayload(TString BTagAlgorithm, TString DataPeriod);

  // SF functions
  float ScaleFactorB(float JetPt, int SystematicFlag);
  float ScaleFactorLight(float JetPt, float JetEta, int SystematicFlag);
  float ScaleFactorJet(int JetFlavor, float JetPt, float JetEta, int SystematicFlag);

  // Efficiency functions
  float TagEfficiencyB_2016(float JetPt, float JetEta);
  float TagEfficiencyC_2016(float JetPt, float JetEta);
  float TagEfficiencyLight_2016(float JetPt, float JetEta);

  float TagEfficiencyB_2017(float JetPt, float JetEta);
  float TagEfficiencyC_2017(float JetPt, float JetEta);
  float TagEfficiencyLight_2017(float JetPt, float JetEta);

  float TagEfficiencyB_2018(float JetPt, float JetEta);
  float TagEfficiencyC_2018(float JetPt, float JetEta);
  float TagEfficiencyLight_2018(float JetPt, float JetEta);

  
  TF1 *funSFb, *funSFlight[4][3];

  TString TaggerName, TaggerOP;
  TString MCSample;

  float TaggerCut;
  float BTagPtBinEdge[50];
  float SFb_error[50];
  float BTagEtaBinEdge[50];
  float FastSimPtBinEdge[50], FastSimEtaBinEdge[50][3];
  float FastSimCF_error[50][2][3]; 
  float FastSimCF[50][2][3];

  int nBTagEtaBins;
  int nBTagPtBins;
  int nFastSimPtBins, nFastSimEtaBins[3];
  int FastSimSystematic;

  int DataYear;

 
  bool period_dependancy;
};
