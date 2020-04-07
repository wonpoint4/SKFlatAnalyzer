#include "BTagSFUtil.h"
#include "BTagCalibrationStandalone.cc"
#include "BTagEfficiencies/BTagEfficienciesTTbar16.C" // Change this if MC are updated
#include "BTagEfficiencies/BTagEfficienciesTTbar17.C" // Change this if MC are updated  
#include "BTagEfficiencies/BTagEfficienciesTTbar18.C" // Change this if MC are updated  

//#################################################
//=== To do list
//#################################################
//--- TODO add fast sim efficiency once one is available  
//--- TODO add AK8 code
//--- TODO add charm tagging code
//#################################################

void BTagSFUtil::SetMCSample(TString s){
  MCSample = s;
}
void BTagSFUtil::SetDataYear(int i){
  DataYear = i;
}

void BTagSFUtil::SetPeriodDependancy(bool b){
  period_dependancy = b;
}

BTagSFUtil::BTagSFUtil(string MeasurementType, string BTagAlgorithm,  TString OperatingPoint, int year, bool pd,int SystematicIndex) {

  //=== initialise global variables
  MCSample = "";
  DataYear = year;
  period_dependancy = pd;
  
  //=== set paths 
  TString datapath = getenv("DATA_DIR");
  TString btagpath = datapath+"/"+TString::Itoa(DataYear,10)+"/BTag/";


  //=== set systematic flags
  string SystematicFlagBC = "central";
  string SystematicFlagL  = "central";
  if (abs(SystematicIndex)<10) {
    if (SystematicIndex==-1 || SystematicIndex==-2) SystematicFlagBC = "down";
    else if (SystematicIndex==+1 || SystematicIndex==+2) SystematicFlagBC = "up";
    else if (SystematicIndex==-3) SystematicFlagL = "down";
    else if (SystematicIndex==+3) SystematicFlagL = "up";
  }


  //=== Initialise TaggerCut, set by taggermap.txt
  TaggerCut = -1;
  ifstream in_tagger(btagpath+"/taggermap.txt");
  string btaggerline;
  while(getline(in_tagger,btaggerline)){
    std::istringstream is_tag( btaggerline );
    TString tstring_taggerline = btaggerline;
    if(tstring_taggerline.Contains("#")) continue;
    int a;
    TString b,c;
    float cutvalue; 
    
    is_tag >> a; // YEAR
    is_tag >> b; // TAGGER
    is_tag >> c; // WP
    is_tag >> cutvalue;  // cut value
   if(a == DataYear && b == BTagAlgorithm && c == OperatingPoint ) {
      TaggerCut = cutvalue; 
      TString flavour_string = "Heavy Flavour"; 
      TString systematic_string = SystematicFlagBC;
      if(MeasurementType=="incl") { flavour_string = "Light Flavour"; systematic_string = SystematicFlagL;}
      cout << "[BTagSFUtil::initializeAnalyzer] " << flavour_string << " : Year = "  << DataYear << "  : OperatingPoint = " << OperatingPoint << " : TaggerCut = " << TaggerCut << " : Systematic " << systematic_string <<  endl;
      break;
    }
  }// end of taggermap loop
  

  //=== Set OperatingPoint passed by var OperatingPoint
  BTagEntry::OperatingPoint op = BTagEntry::OP_LOOSE;
  if (OperatingPoint=="Loose")        {op = BTagEntry::OP_LOOSE;        TaggerOP += "L";}
  if (OperatingPoint=="Medium")       {op = BTagEntry::OP_MEDIUM;        TaggerOP += "M";}
  if (OperatingPoint=="Tight")        {op = BTagEntry::OP_TIGHT;        TaggerOP += "T";}
    

  //=== open histmap file to tell code which file to use depending on year/tagger
  ifstream in(btagpath+"/histmap.txt");
  string btagline;
  if(period_dependancy) cout << "[BTagSFUtil::initializeAnalyzer] setting up Period dependant SF" << endl;
  else  cout << "[BTagSFUtil::initializeAnalyzer] setting up non-period dependant SF" << endl;

  while(getline(in,btagline)){
    std::istringstream is( btagline );
    
    TString tstring_btagline = btagline;
    if(tstring_btagline.Contains("#")) continue;

    string a,f;
    TString b,c,d,e;

    is >> a; // YEAR
    is >> b; // TAGGER
    is >> c; // Period_dep
    is >> d; // Run_start
    is >> e; // Run_end
    is >> f; // csv file                                                                                                                                                                           

    istringstream buffer(a);
    int dy;
    buffer >> dy;

    if( dy != DataYear) continue;
    if( b != BTagAlgorithm) continue;
    
    if( (period_dependancy&& c=="PeriodDep") || (!period_dependancy && c!="PeriodDep")){   
      BTagCalibration calib(BTagAlgorithm, string(btagpath)+f);
      ReaderMap[TString(a)+"_"+b+"_"+d+"_"+e+"_bc"] = new BTagCalibrationReader(&calib, op, MeasurementType, SystematicFlagBC); 
      ReaderMap[TString(a)+"_"+b+"_"+d+"_"+e+"_l"] = new BTagCalibrationReader(&calib, op, MeasurementType, SystematicFlagL);
      
    }
  }
  
  TaggerName = BTagAlgorithm;
  TaggerOP = TaggerName;

  if (TaggerCut==-1) 
    cout << " " << TaggerName << " not supported for " << OperatingPoint << " WP" << endl;

  FastSimSystematic = 0;
  if (abs(SystematicIndex)>10) FastSimSystematic = SystematicIndex%10;

  if (TaggerCut==-1) 
    cout << "BTagSFUtil: " << BTagAlgorithm << " not a supported b-tagging algorithm" << endl;
  
}

BTagSFUtil::~BTagSFUtil() {


  for(map<TString, BTagCalibrationReader*>::iterator it = ReaderMap.begin(); it!= ReaderMap.end(); it++){
    delete it->second;
  }

}

float BTagSFUtil::JetTagEfficiency(int JetFlavor, float JetPt, float JetEta) {

  
  if(DataYear == 2016){
    if (abs(JetFlavor)==5) return TagEfficiencyB_2016(JetPt, JetEta);
    else if (abs(JetFlavor)==4) return TagEfficiencyC_2016(JetPt, JetEta);
    else return TagEfficiencyLight_2016(JetPt, JetEta);
  }
  else  if(DataYear == 2017){
    if (abs(JetFlavor)==5) return TagEfficiencyB_2017(JetPt, JetEta);
    else if (abs(JetFlavor)==4) return TagEfficiencyC_2017(JetPt, JetEta);
    else return TagEfficiencyLight_2017(JetPt, JetEta);
  }
  else  if(DataYear == 2018){
    if (abs(JetFlavor)==5) return TagEfficiencyB_2018(JetPt, JetEta);
    else if (abs(JetFlavor)==4) return TagEfficiencyC_2018(JetPt, JetEta);
    else return TagEfficiencyLight_2018(JetPt, JetEta);
  }

  return -999999.;
}



float BTagSFUtil::GetJetSF(int JetFlavor, float JetPt, float JetEta) {
  float Btag_SF(-999.);

  //=== If chosen period dependant SF, then use lumi to weight

  //=== Currently only CVSv2/DeepCSV 2017 have period dependancy
  //=== If period dependancy set true for other year give warning and set false 

  if (period_dependancy && DataYear != 2017) {
    cout << "[BTagSFUtil::GetJetSF] period dependancy set true for year other than 2017.... This is not available currently...." << endl;
  }
  

  if (period_dependancy) {
    
    if(DataYear == 2016){
      float Btag_SF_BF = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2016_" +TaggerName + "_272007_278808", true, TaggerName);
      float Btag_SF_GH = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2016_" +TaggerName + "_278820_284044", true, TaggerName);
      
      double lumi_periodB = 5.750490644;
      double lumi_periodC = 2.572903489;
      double lumi_periodD = 4.242291557;
      double lumi_periodE = 4.025228137;
      double lumi_periodF = 3.104509132;
      double lumi_periodG = 7.575824256;
      double lumi_periodH = 8.650628380;
      double total_lumi =  (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);
      
      double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
      double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
      Btag_SF = Btag_SF_BF*WeightBtoF + Btag_SF_GH*WeightGtoH;
      return Btag_SF;
    }
    
    else if(DataYear == 2017){

      float Btag_SF_B    = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2017_" +TaggerName + "_297046_299329", true, TaggerName);
      float Btag_SF_CtoE = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2017_" +TaggerName + "_297020_304671", true, TaggerName);
      float Btag_SF_EtoF = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2017_" +TaggerName + "_304671_306462", true, TaggerName);
      
      /// lumi taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2017Analysis and brilcal 
      double lumi_periodB = 4.7939;
      double lumi_periodCtoE = 9.6327 + 4.2478 + 7.8380;
      double lumi_periodEtoF = 1.4765 + 13.5401;

      double total_lumi = (lumi_periodB+lumi_periodCtoE + lumi_periodEtoF);
      
      double WeightB    = lumi_periodB   / total_lumi;
      double WeightCtoE = lumi_periodCtoE/ total_lumi;
      double WeightEtoF = lumi_periodEtoF/ total_lumi;
      Btag_SF = Btag_SF_B*WeightB + Btag_SF_CtoE*WeightCtoE +  Btag_SF_EtoF*WeightEtoF;

      return Btag_SF;

    }
    else if(DataYear == 2018){
      //=== https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X TBA 
    }
  }
  else{
    
    //=== Access SF using Full data year, no periods

    if(DataYear == 2016){

      float Btag_SF_BH = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2016_" +TaggerName + "_272007_284044", false, TaggerName);
      return Btag_SF_BH;
    }
    else if(DataYear == 2017){
      float Btag_SF_BF = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2017_" +TaggerName + "_297046_306462",false,TaggerName);
      return Btag_SF_BF;
    }
    else if(DataYear == 2018){
      float Btag_SF_AD = GetJetSFPeriodDependant(JetFlavor, JetPt, JetEta, "2018_" +TaggerName + "_315252_325175",false,TaggerName);
      return Btag_SF_AD;
    }

  }


  return Btag_SF;
}

float BTagSFUtil::GetJetSFPeriodDependant(int JetFlavor, float JetPt, float JetEta, TString tag_key, bool isPeriodDependant, TString TaggerName) {

  /// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco for pt range and systematic correlations
  float Btag_SF=1.;

  float ThisJetPt = JetPt;

  //=== SF available for 20 < pt  1000 GeV
  if     (JetPt > 999.99) ThisJetPt = 999.99;
  else if(JetPt < 20.   ) ThisJetPt = 20.;
 
  //=== declare iterator to find BTagCalibrationReader in map
  map <TString, BTagCalibrationReader*>::iterator reader_bf_it ;


  if (abs(JetFlavor)==5) {

    //=== access b flavour reader
    reader_bf_it =    ReaderMap.find(tag_key+"_bc");
    if(reader_bf_it == ReaderMap.end()) ErrorLoadingCSV(tag_key, isPeriodDependant, TaggerName);
    else Btag_SF = reader_bf_it->second->eval(BTagEntry::FLAV_B, JetEta, ThisJetPt);
  }
  else if (abs(JetFlavor)==4){ 

    //=== access c flavour reader    
    reader_bf_it = ReaderMap.find(tag_key+"_bc");
    if(reader_bf_it == ReaderMap.end()) ErrorLoadingCSV(tag_key, isPeriodDependant, TaggerName);
    else Btag_SF = reader_bf_it->second->eval(BTagEntry::FLAV_C, JetEta, ThisJetPt);
  }
  else {

    //=== access light flavour reader
    reader_bf_it = ReaderMap.find(tag_key+"_l");
    if(reader_bf_it == ReaderMap.end()) ErrorLoadingCSV(tag_key, isPeriodDependant, TaggerName);
    else Btag_SF = reader_bf_it->second->eval(BTagEntry::FLAV_UDSG, JetEta, ThisJetPt);

  }
    
  return Btag_SF;

}

void BTagSFUtil::ErrorLoadingCSV(TString tag_key, bool isperiodDep, TString TaggerName){

  cout << "[BTagSFUtil::ErrorLoadingCSV] no SF found for key " <<  tag_key << endl;
  cout << "" << endl;
  
  TString year="2016";
  if(tag_key.Contains("2017")) year="2017";
  if(tag_key.Contains("2018")) year="2018";

  //=== load histmap to check possible cvs files and check why error occured

  TString datapath = getenv("DATA_DIR");
  TString btagpath = datapath+"/"+TString::Itoa(DataYear,10)+"/BTag/";
  ifstream in(btagpath+"/histmap.txt");
  string btagline;

  //=== let user know what options have been used in setup
  if(period_dependancy) cout << "[BTagSFUtil::ErrorLoadingCSV] User requested period dependant SF" << endl;
  else  cout << "[BTagSFUtil::ErrorLoadingCSV] User requested non-period dependant SF" << endl;

  cout << "[BTagSFUtil::ErrorLoadingCSV] User running " << year << " year" << endl;
  cout << "[BTagSFUtil::ErrorLoadingCSV] Using tagger: "<< TaggerName << endl;
  cout << "[BTagSFUtil::ErrorLoadingCSV] available SF CSV files from " << btagpath << "/histmap.txt are;"  << endl;
  

  //=== Loop over histmap and output what is available
  int nlines(0);
  bool TaggerFound=false;
  bool isperiodDepAvailable=false;
  while(getline(in,btagline)){
    std::istringstream is( btagline );

    TString tstring_btagline = btagline;
    if(tstring_btagline.Contains("#")) continue;
    nlines++;

    string a,f;
    TString b,c,d,e;

    is >> a; // YEAR                                                                                                                                                                                                                                                            
    is >> b; // TAGGER                                                                                                                                                                                                                                                          
    is >> c; // Period_dep                                                                                                                                                                                                                                                      
    is >> d; // Run_start                                                                                                                                                                                                                                                       
    is >> e; // Run_end                                                                                                                                                                                                                                                         
    is >> f; // csv file                                                                                                                                                                                                                                                        
    if( b == TaggerName) TaggerFound=true;
    if( b == TaggerName && isperiodDep && c == "PeriodDep") isperiodDepAvailable=true;
    cout << "[BTagSFUtil::ErrorLoadingCSV] available CVS files;  Year = "  << a << "  tagger = " << b << " csv file = " << f << endl;

  }
  

  //=== cout some possible errors to help user
  if(nlines == 0) cout << "[BTagSFUtil::ErrorLoadingCSV] error likely that " << btagpath << "/histmap.txt is empty" << endl;
  if(!TaggerFound) cout << "[BTagSFUtil::ErrorLoadingCSV] error likely that tagger " << TaggerName << " is not found in " << btagpath << "/histmap.txt " << endl;
  if(isperiodDep && !isperiodDepAvailable) cout << "[BTagSFUtil::ErrorLoadingCSV] error in setup; user asked for period dependant SF's, which are not available for this Tagger/Year, in Analysis code change 4th argument in SetupBTagger() to false, i.e., SetupBTagger(vtaggers,v_wps, sys, false);" <<endl;
  
  //=== exit if no CSV file is found for user setup
  exit(EXIT_FAILURE);

}


bool BTagSFUtil::IsUncorrectedTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta) {
  /// return false if year is not set                                                                                                                                                                       
  if (DataYear == 0) return false;


  bool isBTagged = JetDiscriminant>TaggerCut;

  return isBTagged;

}


bool BTagSFUtil::IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int seed) {

  //=== set up rand number gen.
  
  TRandom3 rand_(seed);

  /// return false if year is not set
  if (DataYear == 0) return false;


  bool isBTagged = JetDiscriminant>TaggerCut;
  
  //=== Data: no correction needed
  if (JetFlavor==-999999) return isBTagged; 

  
  bool newBTag = isBTagged;

  //=== Get SF 
  JetEta= fabs(JetEta);
  float Btag_SF = GetJetSF(JetFlavor, JetPt, JetEta);

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //=== throw random number to apply correction
  float coin = rand_.Uniform(1.);    
  if(Btag_SF > 1){ 
    //=== use this if SF>1

    if( !isBTagged ) {

      float Btag_eff = JetTagEfficiency(JetFlavor, JetPt, fabs(JetEta));
      //=== fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1./Btag_eff) );
      
      //=== upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  
    //=== use this if SF<1
      
    //=== downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;

}

