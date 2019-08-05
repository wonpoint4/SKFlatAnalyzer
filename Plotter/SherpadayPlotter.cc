#include"Plotter.cc"

class SherpadayPlotter:public Plotter{
public:
  void SetupSamples();
  void SetupSystematics();
  int Setup(int mode_=0);
  int mode;
  SherpadayPlotter();

  double GetChi2(TH1* h1,TH1* h2=NULL);
};
SherpadayPlotter::SherpadayPlotter(){
  vector<TString> files=Split(gSystem->GetFromPipe("find $GENERATORTOOLS_BASE/*/Event -maxdepth 2 -name hists.root -type f"),"\n");
  for(const auto& file:files){
    vector<TString> names=Split(file,"/");
    TString key=names.at(names.size()-4)+"_"+names.at(names.size()-2);
    samplefrags[key]=MakeSampleFrag(key,SampleFrag::Type::SIGNAL,2,make_tuple(file,1.,"",""));
    samplefrags[key+"_EW"]=MakeSampleFrag(key,SampleFrag::Type::SIGNAL,2,make_tuple(file,1.,"","_weight116"));
  }  
  samplefrags["TTW_old"]=MakeSampleFrag("TTW_old",SampleFrag::Type::BG,2,make_tuple("/data9/Users/hsseo/mg/ttW/hists.root",1.,"",""));
}

int SherpadayPlotter::Setup(int mode_){
  samples.clear();
  systematics.clear();
  mode=mode_;
  
  LoadPlots(Form("SherpadayPlotter_mode%d.dat",mode));
  plotdir=Form("SherpadayPlotter_mode%d/",mode);

  SetupSamples();
  SetupSystematics();

  if(DEBUG) cout<<"[Setup] nsample: "<<samples.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;

  return 1;
}

void SherpadayPlotter::SetupSamples(){
  if(DEBUG)  cout<<"[SherpadayPlotter::SetupSamples]"<<endl;
  if(mode==0){
    samples["1Sherpa_ttW_NLO0"]=MakeSample("Sherpa NLO",Sample::Type::SUM,kRed,make_tuple("Sherpa_ttW_NLO0",1.));
    samples["2Sherpa_ttW_NLO0_EW"]=MakeSample("Sherpa NLO + EW",Sample::Type::SUM,kBlue,make_tuple("Sherpa_ttW_NLO0_EW",1.));
    samples["3MG_ttW_NLO0"]=MakeSample("MG NLO",Sample::Type::SUM,kGreen,make_tuple("MG_ttW_NLO0",1.));
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
  systematics["scale"]=MakeSystematic("scale",SystematicType::ENVELOPE,"_weight4 _weight5 _weight6 _weight7 _weight8 _weight9 _weight10",2);
  vector<TString> suffixes;
  for(int i=11;i<111;i++) suffixes.push_back(Form("_weight%d",i));
  systematics["pdf"]=MakeSystematic("pdf",SystematicType::GAUSSIAN,suffixes,2);

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
    
  
