#include"Plotter.cc"

class SherpadayPlotter:public Plotter{
public:
  void SetupEntries();
  void SetupSystematics();
  int Setup(int mode_=0);
  int mode;
  SherpadayPlotter();

  double GetChi2(TH1* h1,TH1* h2=NULL);
};
SherpadayPlotter::SherpadayPlotter(){
  vector<TString> files=Split(gSystem->GetFromPipe("find $GENERATORTOOLS_BASE/Hist/ -maxdepth 2 -name '*.root' -type f"),"\n");
  for(const auto& file:files){
    vector<TString> names=Split(file,"/");
    TString key=Replace(names.back(),".root$","");

    samples[key]=Sample(key,Sample::Type::UNDEF,1);
    samples[key].files.push_back(make_tuple(file,1.,"",""));
  }
}

int SherpadayPlotter::Setup(int mode_){
  entries.clear();
  systematics.clear();
  mode=mode_;
  
  SetupEntries();
  SetupSystematics();
  SetupPlots(Form("plot_Sherpaday/mode%d/SherpaPlotter.dat",mode));

  if(DEBUG) cout<<"[Setup] nentry: "<<entries.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}
void SherpadayPlotter::SetupEntries(){
  if(DEBUG)  cout<<"[SherpadayPlotter::SetupEntries]"<<endl;
  if(mode==0){
    entries.push_back(Sample("Sherpa NLO",Sample::Type::SIGNAL,kRed)+"Sherpa_ttW_NLO0");
    entries.push_back(Sample("MG NLO",Sample::Type::B,kGreen)+"MG_ttW_NLO0");
  }else if(mode==1){
    entries.push_back(Sample("Sherpa NLO CP5 (leptonic)",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0");
    entries.push_back(Sample("MG NLO CP5 (leptonic)",Sample::Type::B,kGreen)+"ttV_leptonic_MG_ttW_leptonic_NLO0");
    entries.push_back(Sample("MG NLO CUETP8M1 (leptonic)",Sample::Type::B,kBlue)+"ttV_leptonic_MG_ttW_leptonic_NLO0_CUETP8M1");
  }else if(mode==2){
    entries.push_back(Sample("Sherpa NLO",Sample::Type::SIGNAL,kRed)+"Sherpa_ttW_NLO0");
    entries.push_back((Sample("Sherpa NLO EW",Sample::Type::A,kBlue)+"Sherpa_ttW_NLO0")%"_weight116");
  }else if(mode==3){
    entries.push_back(Sample("Sherpa NLO (leptonic)",Sample::Type::SIGNAL,kRed)+"Sherpa_ttW_leptonic_NLO0");
    entries.push_back((Sample("Sherpa NLO EW (leptonic)",Sample::Type::A,kBlue)+"Sherpa_ttW_leptonic_NLO0")%"_weight116");
  }else if(mode==11){
    entries.push_back(Sample("Sherpa NLO CP5 (leptonic)",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0");
    entries.push_back((Sample("weight4",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight4");
    entries.push_back((Sample("weight5",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight5");
    entries.push_back((Sample("weight6",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight6");
    entries.push_back((Sample("weight7",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight7");
    entries.push_back((Sample("weight8",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight8");
    entries.push_back((Sample("weight9",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight9");
    entries.push_back((Sample("weight10",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight10");
    entries.push_back((Sample("weight11",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight11");
    entries.push_back((Sample("weight12",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight12");
    entries.push_back((Sample("weight13",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight13");
    entries.push_back((Sample("weight14",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight14");
    entries.push_back((Sample("weight15",Sample::Type::SIGNAL,kRed)+"ttV_leptonic_Sherpa_ttW_leptonic_NLO0")%"_weight15");
    entries.push_back(Sample("MG NLO CP5 (leptonic)",Sample::Type::B,kGreen)+"ttV_leptonic_MG_ttW_leptonic_NLO0");
    entries.push_back(Sample("MG NLO CUETP8M1 (leptonic)",Sample::Type::B,kBlue)+"ttV_leptonic_MG_ttW_leptonic_NLO0_CUETP8M1");
  }else if(mode==12){
    entries.push_back(Sample("Sherpa WW LO (leptonic)",Sample::Type::SIGNAL,kRed)+"WW_leptonic_Sherpa_WW_leptonic_LO0");
    entries.push_back(Sample("MG WW LO (leptonic)",Sample::Type::B,kGreen)+"WW_leptonic_MG_WW_leptonic_LO0");
  }else if(mode==13){
    entries.push_back(Sample("MG WW LO (leptonic)",Sample::Type::B,kGreen)+"WW_leptonic_MG_WW_leptonic_LO0");
    entries.push_back(Sample("MG WW LO wH (leptonic)",Sample::Type::B,kBlue)+"WW_leptonic_MG_WW_leptonic_wH_LO0");
    entries.push_back(Sample("Sherpa WW LO (leptonic)",Sample::Type::SIGNAL,kRed)+"WW_leptonic_Sherpa_WW_leptonic_LO0");
    entries.push_back(Sample("Sherpa WW LO noQED (leptonic)",Sample::Type::A,kOrange)+"WW_leptonic_Sherpa_WW_leptonic_noQED_LO0");
    //entries.push_back(Sample("Sherpa WW LO simple (leptonic)",Sample::Type::A,kOrange)+"WW_leptonic_Sherpa_WW_leptonic_simple_LO0");
    //entries.push_back(Sample("Sherpa WW LO YFS0 (leptonic)",Sample::Type::A,kBlue)+"WW_leptonic_Sherpa_WW_leptonic_YFS0_LO0");
    entries.push_back(Sample("Sherpa WW LO YFS1 (leptonic)",Sample::Type::A,kCyan)+"WW_leptonic_Sherpa_WW_leptonic_YFS1_LO0");
    //entries.push_back(Sample("Sherpa WW LO YFS2 (leptonic)",Sample::Type::A,kYellow)+"WW_leptonic_Sherpa_WW_leptonic_YFS2_LO0");
    //entries.push_back(Sample("Sherpa WW LO noHMS (leptonic)",Sample::Type::A,kMagenta)+"WW_leptonic_Sherpa_WW_leptonic_noHMS_LO0");
    entries.push_back(Sample("Sherpa WW LO EWS1 (leptonic)",Sample::Type::A,kMagenta)+"WW_leptonic_Sherpa_WW_leptonic_EWS1_LO0");
  }    
  
  /*
  if(mode==0){
    samples["1sherpa_NLO1_LO2"]=MakeSample("sherpa 0,1j NLO + 2j LO",Sample::Type::SUM,kRed,make_tuple("ttW",1/1.14942e+08));
    samples["2sherpa_NLO0"]=MakeSample("sherpa NLO",Sample::Type::SUM,kBlue,make_tuple("ttW_NLO0",1/2.9276e+06));
    samples["3sherpa_LO0"]=MakeSample("sherpa LO",Sample::Type::SUM,kGreen,make_tuple("ttW_LO0",1/248615.));
  }else if(mode==1){
    samples["1sherpa_NLO0"]=MakeSample("sherpa NLO",Sample::Type::SUM,kBlue,make_tuple("ttW_NLO0_ext",1/1.16661e+08));
    samplefrags["ttW_NLO0_EW"].type=SampleFrag::Type::BG;
    samples["2sherpa_NLO0_EW"]=MakeSample("sherpa NLO + EWcor",Sample::Type::SUM,kRed,make_tuple("ttW_NLO0_EW",1/1.38038e+08));
  }
  */
    

  //samples["1sherpa_NLO1_LO2"]=MakeSample("sherpa 0,1j NLO + 2j LO",Sample::Type::SUM,kRed,make_tuple("ttW",1/1.14942e+08));
  //samples["2sherpa_NLO1"]=MakeSample("sherpa 0,1j NLO",Sample::Type::SUM,kYellow,make_tuple("ttW_NLO1",1/7.45606e+07));
  //samples["3sherpa_NLO0"]=MakeSample("sherpa 0j NLO",Sample::Type::SUM,kMagenta,make_tuple("ttW_NLO0",1/2.9276e+06));
  //samples["4MG_NLO1"]=MakeSample("MG 0,1j NLO",Sample::Type::SUM,kBlue,make_tuple("TTW",1/91676.1));
  //samples["5MG_NLO0"]=MakeSample("MG 0j NLO",Sample::Type::SUM,kGreen,make_tuple("TTW_NLO0",1/83238.));
  //samples["6MG_LO2"]=MakeSample("MG 0,1,2j LO",Sample::Type::SUM,kBlack,make_tuple("TTW_LO2",1/53641.));
  //samples["7MG_old"]=MakeSample("MG 0,1j NLO old",Sample::Type::SUM,kGreen,make_tuple("TTW_old",1/1259300.3));
}
void SherpadayPlotter::SetupSystematics(){
  if(DEBUG)  cout<<"[SetupSystematics]"<<endl;
  systematics["scale"]=MakeSystematic("scale",Systematic::Type::ENVELOPE,1<<Sample::Type::SIGNAL,"_weight4 _weight5 _weight6 _weight7 _weight8 _weight9 _weight10");
  vector<TString> suffixes;
  for(int i=11;i<111;i++) suffixes.push_back(Form("_weight%d",i));
  systematics["pdf"]=MakeSystematic("pdf",Systematic::Type::GAUSSIAN,1<<Sample::Type::SIGNAL,suffixes);

}
double SherpadayPlotter::GetChi2(TH1* h1,TH1* h2){
  double chi2=0;
  for(int i=h1->GetXaxis()->GetFirst();i<h1->GetXaxis()->GetLast()+1;i++){
    double x1=h1->GetBinContent(i);
    double ex1=h1->GetBinError(i);
    double x2=h2?h2->GetBinContent(i):0.;
    double ex2=h2?h2->GetBinError(i):0.;
    chi2+=pow((x1-x2)/(ex1-ex2),2);
  }
  chi2/=h1->GetXaxis()->GetLast()-h1->GetXaxis()->GetFirst()+1;
  return chi2;
}
    
  
