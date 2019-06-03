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
  TH1::SetDefaultSumw2();
  vector<TString> files=Split(gSystem->GetFromPipe("find /data9/Users/hsseo/SHERPADAY2/sherpa/ -maxdepth 2 -name hists.root -type f"),"\n");
  for(const auto& file:files){
    auto names=Split(file,"/");
    TString key=names.at(names.size()-2);
    samplefrags[key]=MakeSampleFrag(key,SampleFrag::Type::SIGNAL,2,make_tuple(file,"",1.));
  }  
  files=Split(gSystem->GetFromPipe("find /data9/Users/hsseo/SHERPADAY2/MG/ -maxdepth 2 -name hists.root -type f"),"\n");
  for(const auto& file:files){
    auto names=Split(file,"/");
    TString key=names.at(names.size()-2);
    samplefrags[key]=MakeSampleFrag(key,SampleFrag::Type::BG,2,make_tuple(file,"",1.));
  }  
  samplefrags["TTW_old"]=MakeSampleFrag("TTW_old",SampleFrag::Type::BG,2,make_tuple("/data9/Users/hsseo/mg/ttW/hists.root","",1.));
}

int SherpadayPlotter::Setup(int mode_){
  samples.clear();
  systematics.clear();
  mode=mode_;

  SetupSamples();
  SetupSystematics();

  if(DEBUG) cout<<"[Setup] nsample: "<<samples.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;

  return 1;
}

void SherpadayPlotter::SetupSamples(){
  if(DEBUG)  cout<<"[SherpadayPlotter::SetupSamples]"<<endl;
  /*
sherpa/ttH_LO2/ Xsec: 0.0619172 +- 0.00526608 +- 0.0264722 negative fraction: 0.00062 sumw: 7.11313e+06
sherpa/ttH_NLO0/ Xsec: 0.10657 +- 0.00860063 +- 0.0130632 negative fraction: 0.17614 sumw: 386906
sherpa/ttH_NLO1/ still running?
sherpa/ttW/ Xsec: 0.612778 +- 0.0217875 +- 0.0807978 negative fraction: 0.1589 sumw: 1.14942e+08
sherpa/ttW_LO0/ Xsec: 0.365112 +- 0.0123152 +- 0.0857067 negative fraction: 0.00229 sumw: 248615
sherpa/ttW_LO1/ Xsec: 0.438108 +- 0.0140153 +- 0.145301 negative fraction: 0.00284 sumw: 5.93487e+06
sherpa/ttW_LO2/ Xsec: 0.435854 +- 0.0140032 +- 0.16438 negative fraction: 0.00284 sumw: 3.56626e+07
sherpa/ttW_NLO0/ Xsec: 0.573284 +- 0.0183521 +- 0.0736886 negative fraction: 0.11499 sumw: 2.9276e+06
sherpa/ttW_NLO1/ Xsec: 0.637182 +- 0.0211456 +- 0.072207 negative fraction: 0.16202 sumw: 7.45606e+07
sherpa/ttZtoQQorNuNu/ still running?
sherpa/ttZtoQQorNuNu_LO0/ Xsec: 0.457656 +- 0.0141548 +- 0.137168 negative fraction: 0.00104 sumw: 984768
sherpa/ttZtoQQorNuNu_LO1/ Xsec: 0.439962 +- 0.0140459 +- 0.178978 negative fraction: 0.00158 sumw: 1.91396e+07
sherpa/ttZtoQQorNuNu_LO2/ Xsec: 0.437166 +- 0.0140174 +- 0.201324 negative fraction: 0.00206 sumw: 6.25985e+07
sherpa/ttZtoQQorNuNu_NLO0/ Xsec: 0.698832 +- 0.0214304 +- 0.0881039 negative fraction: 0.15223 sumw: 3.35572e+06
sherpa/ttZtoQQorNuNu_NLO1/ Xsec: 0.703369 +- 0.0235521 +- 0.0987803 negative fraction: 0.18576 sumw: 9.27691e+07
MG/TTW/ Xsec: 0.60303 +- 0.0206511 +- 0 negative fraction: 0.243027 sumw: 91676.1
MG/TTW_LO1/ Xsec: 0.39244 +- 0.0130632 +- 0 negative fraction: 0 sumw: 0
MG/TTW_LO2/ Xsec: 0.372087 +- 0.0153835 +- 0 negative fraction: 0 sumw: 53641
MG/TTW_NLO0/ Xsec: 0.5589 +- 0.00487476 +- 0 negative fraction: 0.1683 sumw: 83238
MG/TTZtoQQorNuNu/ Xsec: 0.696217 +- 0.0307786 +- 0 negative fraction: 0.340323 sumw: 0
MG/TTZtoQQorNuNu_LO1/ Xsec: 0.51888 +- 0.0175046 +- 0 negative fraction: 0 sumw: 0
MG/TTZtoQQorNuNu_LO2/ Xsec: 0.450477 +- 0.0198374 +- 0 negative fraction: 0 sumw: 0
MG/TTZtoQQorNuNu_NLO0/ Xsec: 0.6917 +- 0.00574601 +- 0 negative fraction: 0.256767 sumw: 0
  */
  if(mode==0){
    samples["1sherpa_NLO1_LO2"]=MakeSample("sherpa 0,1j NLO + 2j LO",Sample::Type::SUM,kRed,make_tuple("ttW",1/1.14942e+08));
    samples["2sherpa_NLO0"]=MakeSample("sherpa NLO",Sample::Type::SUM,kBlue,make_tuple("ttW_NLO0",1/2.9276e+06));
    samples["3sherpa_LO0"]=MakeSample("sherpa LO",Sample::Type::SUM,kGreen,make_tuple("ttW_LO0",1/248615.));
  }else if(mode==1){
    samples["1sherpa_NLO0"]=MakeSample("sherpa NLO",Sample::Type::SUM,kBlue,make_tuple("ttW_NLO0_ext",1/1.16661e+08));
    samplefrags["ttW_NLO0_EW"].type=SampleFrag::Type::BG;
    samples["2sherpa_NLO0_EW"]=MakeSample("sherpa NLO + EWcor",Sample::Type::SUM,kRed,make_tuple("ttW_NLO0_EW",1/1.38038e+08));
  }
    

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
    
  
