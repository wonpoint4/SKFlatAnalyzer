#include"Plotter.cc"

class EfficiencyPlotter:public Plotter{
public:
  void SetupSamples();
  void SetupSystematics();
  int Setup(int channel_,int year_,TString mode_);
  TString mode;
  int channel,year;
  TString schannel,syear;
  TString analyzer;
  EfficiencyPlotter();

  double GetChi2(TH1* h1,TH1* h2=NULL);
};
EfficiencyPlotter::EfficiencyPlotter(){
  analyzer="EfficiencyValidation";
  TObjArray* arr=gSystem->GetFromPipe("find $SKFlatOutputDir$SKFlatV/"+analyzer+"/ -type f").Tokenize("\n");
  for(int i=0;i<arr->GetEntries();i++){
    TString file=((TObjString*)arr->At(i))->String();
    TString year=file("/"+analyzer+"/[0-9]+/"); year=year(analyzer.Length()+2,4);
    TString key=file(analyzer+"_.*\\.root");
    key=key(analyzer.Length()+1,key.Index(".root")-analyzer.Length()-1)+"_"+year;
    SampleFrag frag;
    frag.files.push_back(make_tuple(file,1.,"",""));
    samplefrags[key]=frag;
  } 
  
  vector<TString> syears={"2016","2017","2018"};
  for(const auto& syear:syears){
    samplefrags["muon"+syear]=MakeSampleFrag("muon"+syear,SampleFrag::Type::DATA,kBlack,TRegexp("SkimTree_SMP_DoubleMuon_[A-Z].*_"+syear));
    samplefrags["electron"+syear]=MakeSampleFrag("electron"+syear,SampleFrag::Type::DATA,kBlack,TRegexp("SkimTree_SMP_.*EG.*_[A-Z].*_"+syear));

    samplefrags["amc"+syear]=MakeSampleFrag("amc"+syear,SampleFrag::Type::SIGNAL,kRed,TRegexp("SkimTree_SMP_DYJets_"+syear));
    samplefrags["vv"+syear]=MakeSampleFrag("Diboson",SampleFrag::Type::BG,kBlue,TRegexp("SkimTree_SMP_[W-Z][W-Z]_pythia_"+syear));
    samplefrags["wjets"+syear]=MakeSampleFrag("W",SampleFrag::Type::BG,kYellow,"SkimTree_SMP_WJets_MG_"+syear);
    samplefrags["tt"+syear]=MakeSampleFrag("t#bar{t}",SampleFrag::Type::BG,kMagenta,"SkimTree_SMP_TTLL_powheg_"+syear);
  }

  vector<TString> ids={"POGTight_PFIsoTight","POGTight_TrkIsoLoose","MediumID","MediumID_selective","MediumID_Q","MediumID_selective_Q"};
  for(const auto& element:samplefrags){
    const SampleFrag& frag=element.second;
    if(frag.type!=SampleFrag::Type::UNDEFINE&&element.first.Contains(TRegexp("201[0-9]$"))){
      cout<<"ss"<<endl;
      for(const auto& id:ids){
	samplefrags[element.first+"_"+id]=MakeSampleFrag(frag.title,frag.type,frag.linecolor,element.first,1.,"","_"+id);
      }
    }
  }
}

int EfficiencyPlotter::Setup(int channel_,int year_,TString mode_=""){
  samples.clear();
  systematics.clear();
  plots.clear();
  channel=channel_;
  year=year_;
  mode=mode_;
  schannel=GetStringChannel((Channel)channel);
  syear=Form("%d",year);

  SetupSamples();
  SetupSystematics();
  //SetupPlots();

  if(DEBUG) cout<<"[Setup] nsample: "<<samples.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}

void EfficiencyPlotter::SetupSamples(){
  if(DEBUG)  cout<<"[EfficiencyPlotter::SetupSamples]"<<endl;
  vector<tuple<int,int,TString>> availables={make_tuple(0,2016,"POGTight_PFIsoTight"),
					     make_tuple(0,2017,"POGTight_PFIsoTight"),
					     make_tuple(0,2017,"POGTight_TrkIsoLoose_Q"),
					     make_tuple(0,2018,"POGTight_PFIsoTight"),
					     make_tuple(0,2018,"POGTight_TrkIsoLoose"),
					     make_tuple(1,2016,"MediumID"),
					     make_tuple(1,2016,"MediumID_selective_Q"),
					     make_tuple(1,2017,"MediumID"),
					     make_tuple(1,2017,"MediumID_Q"),
					     make_tuple(1,2017,"MediumID_selective_Q"),
					     make_tuple(1,2018,"MediumID"),
					     make_tuple(1,2018,"MediumID_selective"),};
  
  for(const auto& avail:availables){
    if(make_tuple(channel,year,mode)==avail){
      samples["data"]=MakeSample("data",Sample::Type::SUM,kBlack,make_tuple(schannel+syear+"_"+mode,1.));
      samples["data"].markerstyle=20;
      samples["data"].markersize=0.7;
      if(year==2018) samples["sim"]=MakeSample("sim",Sample::Type::STACK,kBlue,make_tuple("mg"+syear+"_"+mode,1.),make_tuple("mgtt"+syear+"_"+mode,1.),make_tuple("vv"+syear+"_"+mode,1.),make_tuple("wjets"+syear+"_"+mode,1.),make_tuple("tt"+syear+"_"+mode,1.));
      else samples["sim"]=MakeSample("sim",Sample::Type::STACK,kBlue,make_tuple("amc"+syear+"_"+mode,1.),make_tuple("amctt"+syear+"_"+mode,1.),make_tuple("vv"+syear+"_"+mode,1.),make_tuple("wjets"+syear+"_"+mode,1.),make_tuple("tt"+syear+"_"+mode,1.));
      get<0>(samples["sim"].frags[0]).title=(channel==Channel::MUON)?"#gamma*/Z#rightarrow#mu#mu":"#gamma*/Z#rightarrowee";
      return;
    }
  }
  cout<<"###ERROR#### [EfficiencyPlotter::SetupSamples] Not available configuration"<<endl;
  return;
}
void EfficiencyPlotter::SetupSystematics(){
  if(DEBUG)  cout<<"[SetupSystematics]"<<endl;
  if(channel==Channel::ELECTRON) systematics["RECOSF"]=MakeSystematic("RECOSF",SystematicType::ENVELOPE,"_RECOSF_up _RECOSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["IDSF"]=MakeSystematic("IDSF",SystematicType::ENVELOPE,"_IDSF_up _IDSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::MUON) systematics["ISOSF"]=MakeSystematic("ISOSF",SystematicType::ENVELOPE,"_ISOSF_up _ISOSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["triggerSF"]=MakeSystematic("triggerSF",SystematicType::ENVELOPE,"_triggerSF_up _triggerSF_down",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));

  if(channel==Channel::ELECTRON) systematics["noRECOSF"]=MakeSystematic("noRECOSF",SystematicType::ENVELOPE,"_noRECOSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["noIDSF"]=MakeSystematic("noIDSF",SystematicType::ENVELOPE,"_noIDSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  if(channel==Channel::MUON) systematics["noISOSF"]=MakeSystematic("noISOSF",SystematicType::ENVELOPE,"_noISOSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["notriggerSF"]=MakeSystematic("notriggerSF",SystematicType::ENVELOPE,"_notriggerSF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));

  systematics["nozptcor"]=MakeSystematic("nozptcor",SystematicType::ENVELOPE,"_nozptcor",(1<<SampleFrag::Type::SIGNAL));
  systematics["noefficiencySF"]=MakeSystematic("noefficiencySF",SystematicType::ENVELOPE,"_noefficiencySF",(1<<SampleFrag::Type::SIGNAL)+(1<<SampleFrag::Type::BG));
  systematics["efficiencySF"]=MakeSystematic("efficiencySF",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF",0);

}

double EfficiencyPlotter::GetChi2(TH1* h1,TH1* h2){
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
    
  
