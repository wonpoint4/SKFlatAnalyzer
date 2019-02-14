#include<iostream>
#include<vector>
#include"TH2D.h"
#include"TFile.h"
#include"TKey.h"
#include"TString.h"
#include"THStack.h"
#include"TLegend.h"
#include"TCanvas.h"
#include"TLine.h"
using namespace std;
/////////////////////////////////////////////////////////////////////////////
///////////////////////////// struct and enum ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////
struct Sample{
  TString name;
  int type;
  int color;
  vector<TString> files;
};
enum SuffixType{NONE,UPDOWN,NUMBER,COSTOM};
TString GetStringSuffixType(SuffixType type){
  switch(type){
  case NONE: return "NONE";
  case UPDOWN: return "UPDOWN";
  case NUMBER: return "NUMBER";
  case COSTOM: return "CUSTON";
  default: return "###WARNING### Bad SuffixType";
  }
}
enum SystematicType{ENVELOPE,GAUSSIAN,HESSIAN,MULTI};
TString GetStringSystematicType(SystematicType type){
  switch(type){
  case ENVELOPE: return "ENVELOPE";
  case GAUSSIAN: return "GAUSSIAN";
  case HESSIAN: return "HESSIAN";
  case MULTI: return "MULTI";
  default: return "###WARNING### Bad SystematicType";
  }
}
struct Systematic{
  TString name;
  SystematicType type;
  SuffixType suffixtype;
  vector<TString> suffixes;
  bool exist_data;
  bool exist_bg;
  int sysbit;
};
struct Histogram{
  TString name;
  int rebin;
  double xmin;
  double xmax;
  TString option;
};
struct Directory{
  TString name;
  vector<Histogram> histograms;
};
enum SampleType{DATA,DY,TAU,BG};
TString GetStringSampleType(SampleType type){
  switch(type){
  case DATA: return "DATA";
  case DY: return "DY";
  case TAU: return "TAU";
  case BG: return "BG";
  default: return "###WARNING### Bad SampleType";
  }
}
TString GetStringEColor(EColor color){
  switch(color){
  case kBlack: return "kBlack";
  case kRed: return "kRed";
  case kGreen: return "kGreen";
  case kBlue: return "kBlue";
  case kYellow: return "kYellow";
  case kMagenta: return "kMagenta";
  default: return "###WARNING### Bad EColor";
  }
}
enum Channel{MUON,ELECTRON};
TString GetStringChannel(Channel channel){
  switch(channel){
  case MUON: return "muon";
  case ELECTRON: return "electron";
  default: return "###WARNING### Bad Channel";
  }
}

/////////////////////////////////////////////////////////////////////////////
///////////////////////////// global variables ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////
vector<Sample> samples;
vector<Systematic> systematics;
vector<Directory> directories;
bool DEBUG=false;

/////////////////////////////////////////////////////////////////////////////
//////////////////// Add functions for global variables//////////////////////
/////////////////////////////////////////////////////////////////////////////
void AddSample(TString name_,SampleType type_,EColor color_,vector<TString> files_){
  Sample sample;
  sample.name=name_;
  sample.type=type_;
  sample.color=color_;
  sample.files=files_;
  cout<<" [AddSample] "<<sample.name<<" "<<GetStringSampleType((SampleType)sample.type)<<" "<<GetStringEColor((EColor)sample.color)<<endl;
  for(int i=0;i<sample.files.size();i++)
    cout<<"   "<<sample.files.at(i)<<endl;    
  samples.push_back(sample);
}
void AddSample(TString name_,SampleType type_,EColor color_,TString file1,TString file2="",TString file3="",TString file4="",TString file5="",TString file6="",TString file7=""){
  vector<TString> files;
  if(file1!="") files.push_back(file1);
  if(file2!="") files.push_back(file2);
  if(file3!="") files.push_back(file3);
  if(file4!="") files.push_back(file4);
  if(file5!="") files.push_back(file5);
  if(file6!="") files.push_back(file6);
  if(file7!="") files.push_back(file7);
  AddSample(name_,type_,color_,files);
}
void AddSystematic(TString name_,SystematicType type_,vector<TString> includes,bool exist_data_,bool exist_bg_){
  Systematic systematic;
  systematic.name=name_;
  systematic.type=type_;
  if(systematic.type==SystematicType::MULTI){
    systematic.sysbit=0;
    for(int i=0;i<systematics.size();i++)
      for(int j=0;j<includes.size();j++)
	if(systematics[i].name==includes[j]){
	  systematic.sysbit|=systematics[i].sysbit;
	  break;
	}
  }else{
    systematic.sysbit=1<<systematics.size();
    systematic.suffixes=includes;
  }    
  systematic.exist_data=exist_data_;
  systematic.exist_bg=exist_bg_;
  cout<<" [AddSystematic] "<<systematic.name<<" "<<GetStringSystematicType(systematic.type)<<" sysbit:"<<systematic.sysbit<<" data:"<<(systematic.exist_data?"exist":"not_exist")<<" bg:"<<(systematic.exist_bg?"exist":"not_exist")<<endl;
  cout<<"  INCLUDE=";for(int i=0;i<includes.size();i++) cout<<includes[i]<<" ";
  cout<<endl;
  systematics.push_back(systematic);
}
void AddSystematic(TString name_,SystematicType type_,TString includes_,bool exist_data_,bool exist_bg_){
  TObjArray* arr=includes_.Tokenize(" ");
  vector<TString> includes;
  for(int i=0;i<arr->GetEntries();i++){
    includes.push_back(((TObjString*)arr->At(i))->String());
  }
  AddSystematic(name_,type_,includes,exist_data_,exist_bg_);
}

void AddHistogram(Directory& directory,TString name_,int rebin_,double xmin_,double xmax_,TString option_=""){
  Histogram histogram;
  histogram.name=name_;
  histogram.rebin=rebin_;
  histogram.xmin=xmin_;
  histogram.xmax=xmax_;
  histogram.option=option_;
  if(DEBUG) std::cout<<" [AddHistogram] to "<<histogram.name<<" "<<histogram.rebin<<" "<<histogram.xmin<<" "<<histogram.xmax<<endl;
  directory.histograms.push_back(histogram);
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Setup functions///////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void SetupSamples(int channel,int year,TString analyzer,TString skim){
  TString syear=Form("%d",year);
  if(skim!="") skim="_"+skim;
  cout<<"[SetupSamples] "<<analyzer<<" "<<GetStringChannel((Channel)channel)<<year<<endl;
  TString filedir=TString(getenv("SKFlatOutputDir"))+getenv("SKFlatV")+"/"+analyzer+"/";
  if(year==2017){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleMuon_F.root");
      AddSample("#gamma*/Z#rightarrow#mu#mu",SampleType::DY,EColor::kRed,filedir+"2017/"+analyzer+skim+"_DYJets.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_B.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2017/DATA/"+analyzer+skim+"_DoubleEG_F.root");
      AddSample("#gamma*/Z#rightarrowee",SampleType::DY,EColor::kRed,filedir+"2017/"+analyzer+skim+"_DYJets.root");
    }else return "";
  }else if(year==2016){
    if(channel==Channel::MUON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleMuon_H.root");
      AddSample("#gamma*/Z#rightarrow#mu#mu",SampleType::DY,(EColor)(EColor::kRed),filedir+"2016/"+analyzer+skim+"_DYJets.root");
    }else if(channel==Channel::ELECTRON){
      AddSample("data",SampleType::DATA,EColor::kBlack,filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_B_ver2.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_C.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_D.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_E.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_F.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_G.root",filedir+"2016/DATA/"+analyzer+skim+"_DoubleEG_H.root");
      AddSample("#gamma*/Z#rightarrowee",SampleType::DY,EColor::kRed,filedir+"2016/"+analyzer+skim+"_DYJets.root");
    }else return "";
  }else return "";
  AddSample("#gamma*/Z#rightarrow#tau#tau",SampleType::TAU,EColor::kGreen,filedir+syear+"/"+analyzer+skim+"_DYJets.root");
  AddSample("Diboson",SampleType::BG,EColor::kBlue,filedir+syear+"/"+analyzer+skim+"_WW_pythia.root",filedir+syear+"/"+analyzer+skim+"_WZ_pythia.root",filedir+syear+"/"+analyzer+skim+"_ZZ_pythia.root");
  AddSample("W",SampleType::BG,EColor::kYellow,filedir+syear+"/"+analyzer+skim+"_WJets_MG.root");
  AddSample("t#bar{t}",SampleType::BG,EColor::kMagenta,year==2017?filedir+"2017/"+analyzer+skim+"_TTLL_powheg.root":filedir+"2016/"+analyzer+skim+"_TT_powheg.root");
}
void SetupSystematics(int channel,int year,TString analyzer){
  cout<<"[SetupSystematics]"<<endl;
  AddSystematic("RECOSF",SystematicType::ENVELOPE,"_RECOSF_up _RECOSF_down",0,1);
  AddSystematic("IDSF",SystematicType::ENVELOPE,"_IDSF_up _IDSF_down",0,1);
  AddSystematic("ISOSF",SystematicType::ENVELOPE,"_ISOSF_up _ISOSF_down",0,1);
  AddSystematic("triggerSF",SystematicType::ENVELOPE,"_triggerSF_up _triggerSF_down",0,1);
  AddSystematic("PUreweight",SystematicType::ENVELOPE,"_PUreweight_up _PUreweight_down",0,1);
  AddSystematic("prefireweight",SystematicType::ENVELOPE,"_prefireweight_up _prefireweight_down",0,1);
  AddSystematic("alphaS",SystematicType::ENVELOPE,"_alphaS_up _alphaS_down",0,0);
  AddSystematic("scalevariation",SystematicType::ENVELOPE,"_scalevariation0 _scalevariation1 _scalevariation2 _scalevariation3 _scalevariation5 _scalevariation6 _scalevariation7",0,0);
  
  vector<TString> prefixes;
  for(int i=0;i<100;i++) prefixes.push_back(Form("_pdf%d",i));
  if(year==2017) AddSystematic("pdf",SystematicType::HESSIAN,prefixes,0,0);
  else if(year==2016) AddSystematic("pdf",SystematicType::GAUSSIAN,prefixes,0,0);
  else cout<<"###WARNING### [SetupSystematics] wrong year"<<endl;

  AddSystematic("noRECOSF",SystematicType::ENVELOPE,"_noRECOSF",0,1);
  AddSystematic("noIDSF",SystematicType::ENVELOPE,"_noIDSF",0,1);
  AddSystematic("noISOSF",SystematicType::ENVELOPE,"_noISOSF",0,1);
  AddSystematic("notriggerSF",SystematicType::ENVELOPE,"_notriggerSF",0,1);
  AddSystematic("noPUreweight",SystematicType::ENVELOPE,"_noPUreweight",0,1);
  AddSystematic("noprefireweight",SystematicType::ENVELOPE,"_noprefireweight",0,1);
  AddSystematic("nozptcor",SystematicType::ENVELOPE,"_nozptcor",0,0);
  AddSystematic("noefficiencySF",SystematicType::ENVELOPE,"_noefficiencySF",0,1);
  AddSystematic("IDSF_POG",SystematicType::ENVELOPE,"_IDSF_POG",0,1);
  AddSystematic("efficiencySF",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF",0,0);
  AddSystematic("totalsys",SystematicType::MULTI,"RECOSF IDSF ISOSF triggerSF PUreweight prefireweight alphaS scalevariation pdf nozptcor",0,0);
}
void SetupDirectoriesSMPValidation(int channel,int year){
  cout<<"[SetupDirectoriesSMPValidation] SMPValidation "<<GetStringChannel((Channel)channel)<<year<<endl;
  TString region[]={"OS","OS_Z","OS_Z_y0.0to0.4","OS_Z_y0.4to0.8","OS_Z_y0.8to1.2","OS_Z_y1.2to1.6","OS_Z_y1.6to2.0","OS_Z_y2.0to2.4","SS"};
  for(int i=0;i<sizeof(region)/sizeof(TString);i++){
    Directory directory;
    directory.name=GetStringChannel((Channel)channel)+Form("%d",year)+"/"+region[i]+"/";
    cout<<directory.name<<" ";
    AddHistogram(directory,"dimass",0,80,100);
    AddHistogram(directory,"dipt",4,0,200);
    AddHistogram(directory,"dirap",0,0,0);
    AddHistogram(directory,"l0pt",2,0,100);
    AddHistogram(directory,"l0eta",0,0,0);
    AddHistogram(directory,"l0riso",0,0,0);
    AddHistogram(directory,"l1pt",2,0,100);
    AddHistogram(directory,"l1eta",0,0,0);
    AddHistogram(directory,"l1riso",0,0,0);
    AddHistogram(directory,"lldelR",0,0,0);
    AddHistogram(directory,"lldelphi",0,0,0);
    AddHistogram(directory,"nPV",0,0,0);
    directories.push_back(directory);
  }
  cout<<endl;
  directories[0].histograms[0].rebin=4;
  directories[0].histograms[0].xmin=0;
  directories[0].histograms[0].xmax=0;
}
void SetupDirectoriesAFBAnalyzer(int channel,int year){
  cout<<"[SetupDirectoriesAFBAnalyzer] AFBAnalyzer "<<GetStringChannel((Channel)channel)<<year<<endl;
  const int ybinnum=6;
  double yrange[ybinnum+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  const int mbinnum=12;
  double massrange[mbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120};
  vector<TString> region;
  region.push_back(Form("OS_m%.0fto%.0f",massrange[0],massrange[mbinnum]));
  for(int iy=0;iy<ybinnum;iy++){
    region.push_back(Form("OS_m%.0fto%.0f_y%.1fto%.1f",massrange[0],massrange[mbinnum],yrange[iy],yrange[iy+1]));
    for(int im=0;im<mbinnum;im++){
      //region.push_back(Form("OS_m%.0fto%.0f_y%.1fto%.1f",massrange[im],massrange[im+1],yrange[iy],yrange[iy+1]));
    }
  }
  for(int i=0;i<region.size();i++){
    Directory directory;
    directory.name=GetStringChannel((Channel)channel)+Form("%d",year)+"/"+region[i]+"/";
    cout<<directory.name<<" ";
    AddHistogram(directory,"dimass",2,60,120);
    AddHistogram(directory,"dipt",4,0,200);
    AddHistogram(directory,"dirap",0,0,0);
    AddHistogram(directory,"l0pt",2,0,100);
    AddHistogram(directory,"l0eta",0,0,0);
    AddHistogram(directory,"l0riso",0,0,0);
    AddHistogram(directory,"l1pt",2,0,100);
    AddHistogram(directory,"l1eta",0,0,0);
    AddHistogram(directory,"l1riso",0,0,0);
    AddHistogram(directory,"lldelR",0,0,0);
    AddHistogram(directory,"lldelphi",0,0,0);
    AddHistogram(directory,"nPV",0,0,0);
    AddHistogram(directory,"costhetaCS",0,0,0,"1:noleg");
    AddHistogram(directory,"abscosthetaCS",0,0,0,"1:noleg");
    if(region[i].Contains("m60to120")) AddHistogram(directory,"AFB",0,0,0,"AFB");
    directories.push_back(directory);
  }
  cout<<endl;
}
TString Setup(int channel,int year,TString analyzer="SMPValidation",TString skim="SkimTree_SMP"){
  samples.clear();
  systematics.clear();
  directories.clear();

  SetupSamples(channel,year,analyzer,skim);
  SetupSystematics(channel,year,analyzer);
  if(analyzer.Contains("SMPValidation")) SetupDirectoriesSMPValidation(channel,year);
  else if(analyzer.Contains("AFBAnalyzer")) SetupDirectoriesAFBAnalyzer(channel,year);
  else cout<<"###WARNING### [Setup] unknown Analyzer="<<analyzer<<endl;

  cout<<"[Setup] nsample: "<<samples.size()<<endl;
  cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  cout<<"[Setup] ndirectories: "<<directories.size()<<endl;  

  TString schannel=GetStringChannel((Channel)channel);
  TString syear=Form("%d",year);
  return schannel+syear;
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Core functions////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TH1* GetHist(TString filename,TString histname){
  TFile f(filename);
  TH1* hist=(TH1*)f.Get(histname);
  if(hist){
    hist->SetDirectory(0);
    return hist;
  }else{
    cout<<"###WARNING### [GetHist] no "<<histname<<" in "<<filename<<endl;
    return NULL;
  }
}
vector<TH1*> GetHists(TString datahistname,TString dyhistname,TString bghistname,TString option){
  vector<TH1*> hists;
  if(option.Contains("AFB")){
    TString datahistname_forward=datahistname;
    TString dyhistname_forward=dyhistname;
    TString bghistname_forward=bghistname;
    datahistname_forward.ReplaceAll("AFB","forward");
    dyhistname_forward.ReplaceAll("AFB","forward");
    bghistname_forward.ReplaceAll("AFB","forward");
    TString datahistname_backward=datahistname;
    TString dyhistname_backward=dyhistname;
    TString bghistname_backward=bghistname;
    datahistname_backward.ReplaceAll("AFB","backward");
    dyhistname_backward.ReplaceAll("AFB","backward");
    bghistname_backward.ReplaceAll("AFB","backward");
    hists=GetHists(datahistname_forward,dyhistname_forward,bghistname_forward,"");
    vector<TH1*> hists_backward=GetHists(datahistname_backward,dyhistname_backward,bghistname_backward,"");
    hists.insert(hists.end(),hists_backward.begin(),hists_backward.end());
  }
  else{
    for(int i=0;i<samples.size();i++){
      TH1* hist=NULL;
      TString histname;
      if(samples[i].type==SampleType::DATA) histname=datahistname;
      else if(samples[i].type==SampleType::DY) histname=dyhistname;
      else if(samples[i].type==SampleType::TAU) histname=bghistname(0,bghistname.Last('/')+1)+"tau_"+bghistname(bghistname.Last('/')+1,bghistname.Length());
    else histname=bghistname;
      for(int j=0;j<samples[i].files.size();j++){
	if(!hist) hist=GetHist(samples[i].files[j],histname);
	else hist->Add(GetHist(samples[i].files[j],histname));
      }
      if(hist){
	hist->SetName(samples[i].name);
	hist->SetTitle(dyhistname);
	hist->GetXaxis()->SetTitle(dyhistname);
	hist->SetLineColor(samples[i].color);
	hist->SetFillColor(samples[i].color);
      }
      hists.push_back(hist);
    }
  }
  return hists;
}
vector<TH1*> GetHists(TString histname,TString option){
  return GetHists(histname,histname,histname,option);
}
TH1* GetHistAFB(TH1* hist_forward,TH1* hist_backward){
  TH1* hist=(TH1*)hist_forward->Clone("afb");
  hist->Reset();
  for(int i=0;i<hist->GetNbinsX()+2;i++){
    double valf=hist_forward->GetBinContent(i);
    double ef=hist_forward->GetBinError(i);
    double valb=hist_backward->GetBinContent(i);
    double eb=hist_backward->GetBinError(i);
    hist->SetBinContent(i,(valf-valb)/(valf+valb));
    hist->SetBinError(i,2*sqrt(ef*ef*valb*valb+eb*eb*valf*valf)/pow(valf+valb,2));
  }
  return hist;
}
TH1* GetHMC(const vector<TH1*>& hists,TString option=""){
  if(option.Contains("AFB")){
    TH1* hist_forward=hists.at(1);
    TH1* hist_backward=hists.at(samples.size()+1);
    TH1* hist=GetHistAFB(hist_forward,hist_backward);
    hist->SetName("DY");
    return hist;        
  }else if(option.Contains("stack")){
    THStack* hstack=new THStack;
    for(int i=(int)hists.size()-1;i>0;i--){
      if(hists.at(i)){
	if(samples.at(i).type==SampleType::DY) hists.at(i)->SetFillColor(samples.at(i).color-9);
	hstack->Add(hists.at(i),"HIST");
      }
    }
    return (TH1*)hstack;
  }else if(option.Contains("BGSub")) return (TH1*)hists.at(1)->Clone();
  else{
    TH1* hist=(TH1*)hists.at(1)->Clone("hmc");
    for(int i=2;i<(int)hists.size();i++) if(hists.at(i)) hist->Add(hists.at(i));
    return hist;
  }
}
TH1* GetHMC(THStack* hstack){
  TList* list=hstack->GetHists();
  TH1* hist=(TH1*)list->At(list->GetSize()-1)->Clone("hmc");
  for(int i=0;i<list->GetSize()-1;i++){
    hist->Add((TH1*)list->At(i));
  }
  return hist;
}
TH1* GetHData(const vector<TH1*>& hists,TString option=""){
  if(option.Contains("AFB")){
    vector<TH1*> hists_forward,hists_backward;
    hists_forward.insert(hists_forward.begin(),hists.begin(),hists.begin()+samples.size());
    hists_backward.insert(hists_backward.begin(),hists.begin()+samples.size(),hists.end());
    TH1* hist_forward=GetHData(hists_forward,"BGSub");
    TH1* hist_backward=GetHData(hists_backward,"BGSub");
    TH1* hist=GetHistAFB(hist_forward,hist_backward);
    delete hist_forward;delete hist_backward;
    hist->SetName("data");
    return hist;    
  }else if(option.Contains("BGSub")){
    TH1* hist=(TH1*)hists.at(0)->Clone("hdata");
    for(int i=2;i<(int)samples.size();i++) if(hists.at(i)) hist->Add(hists.at(i),-1.);
    return hist;
  }
  else return (TH1*)hists.at(0)->Clone();
}
TLegend* GetLegend(const vector<TH1*>& hists,TString option){
  int histssize=hists.size();
  if(option.Contains("BGSub")) histssize=2;
  double horizontalshift=0;
  if(option.Contains("leftleg")) horizontalshift=-0.53;
  TLegend* legend=new TLegend(0.67+horizontalshift,0.88-histssize*0.07,0.89+horizontalshift,0.88);
  for(int i=0;i<histssize;i++){
    TString att="f";
    if(i==0) att="lp";
    legend->AddEntry(hists.at(i),hists.at(i)->GetName(),att);
  }
  legend->SetBorderSize(0);
  return legend;
}
TLegend* GetLegend(TH1* h1,TH1* h2,TString option){
  vector<TH1*> hists;
  hists.push_back(h1);
  if(strstr(h2->ClassName(),"THStack")){
    TList* list=((THStack*)h2)->GetHists();
    for(int i=list->GetSize()-1;i>=0;i--){
      hists.push_back((TH1*)list->At(i));
    }
  }else hists.push_back(h2);
  return GetLegend(hists,option);
}
TH1* GetEnvelope(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetHMC((THStack*)central);
  TH1* syshist=(TH1*)central->Clone("sys");
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,0);
  for(int i=0;i<(int)variations.size();i++){
    for(int j=0;j<syshist->GetNbinsX()+1;j++){
      double diff=fabs(syshist->GetBinContent(j)-variations.at(i)->GetBinContent(j));
      if(diff>syshist->GetBinError(j)) syshist->SetBinError(j,diff);
    }
  }
  return syshist;
}
TH1* GetEnvelope(TH1* central,TH1* variation1,TH1* variation2=NULL,TH1* variation3=NULL,TH1* variation4=NULL,TH1* variation5=NULL,TH1* variation6=NULL,TH1* variation7=NULL,TH1* variation8=NULL,TH1* variation9=NULL){
  vector<TH1*> variations;
  if(variation1) variations.push_back(variation1);
  if(variation2) variations.push_back(variation2);
  if(variation3) variations.push_back(variation3);
  if(variation4) variations.push_back(variation4);
  if(variation5) variations.push_back(variation5);
  if(variation6) variations.push_back(variation6);
  if(variation7) variations.push_back(variation7);
  if(variation8) variations.push_back(variation8);
  if(variation9) variations.push_back(variation9);
  return GetEnvelope(central,variations);
}    
TH1* GetHessianError(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetHMC((THStack*)central);
  TH1* syshist=(TH1*)central->Clone("sys");
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,0);
  for(int i=0;i<(int)variations.size();i++){
    for(int j=0;j<syshist->GetNbinsX()+1;j++){
      double diff=fabs(syshist->GetBinContent(j)-variations.at(i)->GetBinContent(j));
      syshist->SetBinError(j,sqrt(pow(syshist->GetBinError(j),2)+pow(diff,2)));
    }
  }
  return syshist;
}  
TH1* GetRMSError(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetHMC((THStack*)central);
  TH1* syshist=(TH1*)central->Clone("sys");
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,0);
  for(int i=0;i<(int)variations.size();i++){
    for(int j=0;j<syshist->GetNbinsX()+1;j++){
      double diff=fabs(syshist->GetBinContent(j)-variations.at(i)->GetBinContent(j));
      syshist->SetBinError(j,sqrt(pow(syshist->GetBinError(j),2)+pow(diff,2)));
    }
  }
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,syshist->GetBinError(i)/sqrt(variations.size()));
  return syshist;
}  
int AddError(TH1* hist,TH1* sys){
  for(int i=1;i<hist->GetNbinsX()+1;i++){
    if(fabs(hist->GetBinContent(i)-sys->GetBinContent(i))*1000000>fabs(hist->GetBinContent(i))){
      cout<<"###WARNING### [AddError] sys hist is wrong"<<endl;
      cout.precision(20);
      cout<<i<<" "<<hist->GetBinContent(i)<<" "<<sys->GetBinContent(i)<<" "<<fabs(hist->GetBinContent(i)-sys->GetBinContent(i))<<endl;
      return -1;
    }
  }
  for(int i=1;i<hist->GetNbinsX()+1;i++){
    double err1=hist->GetBinError(i);
    double err2=sys->GetBinError(i);
    hist->SetBinError(i,sqrt(err1*err1+err2*err2));
  }
  return 1;
}
TCanvas* GetCompare(TH1* h1,TH1* h1sys,TH1* h2,TH1* h2sys,int rebin,double xmin,double xmax,TString option){
  THStack *hstack=NULL;
  if(strstr(h2->ClassName(),"THStack")!=NULL){
    hstack=(THStack*)h2;
    h2=GetHMC(hstack);
  }
  if(rebin){
    h1->Rebin(rebin);
    if(h1sys) h1sys->Rebin(rebin);
    if(hstack){
      TList* list=hstack->GetHists();
      for(int i=0;i<list->GetSize();i++){
	((TH1*)list->At(i))->Rebin(rebin);
      }
    }
    h2->Rebin(rebin);
    if(h2sys) h2sys->Rebin(rebin);
  }
  if(xmin||xmax){
    h1->GetXaxis()->SetRangeUser(xmin,xmax);
    h2->GetXaxis()->SetRangeUser(xmin,xmax);
  }

  h1->SetStats(0);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.7);
  TH1 *h1total=NULL,*h2total=NULL;
  if(h1sys){
    h1total=(TH1*)h1->Clone("h1total");
    AddError(h1total,h1sys);
  }else{
    h1total=h1;
  }
  if(h2sys){
    h2total=(TH1*)h2->Clone("h2total");
    AddError(h2total,h2sys);
  }else{
    h2total=h2;
  }

  TCanvas* c1=new TCanvas(h2->GetTitle(),h2->GetTitle(),800,800);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);

  h1total->Draw("e");
  h1total->GetXaxis()->SetLabelSize(0);
  h1total->GetXaxis()->SetTitle("");
  if(h1sys) h1->Draw("same e1");

  if(hstack){
    hstack->Draw("same");
    h2->SetFillStyle(3144);
    h2->SetFillColor(samples.at(1).color+2);
    if(h2sys){
      h2total->SetFillStyle(3244);
      h2total->SetFillColor(samples.at(1).color);
    }      
  }else if(h2sys){
    h2->SetFillStyle(3001);
    h2->SetFillColor(samples.at(1).color-9);
    h2total->SetFillStyle(3001);
    h2total->SetFillColor(samples.at(1).color+1);
  }else{
    h2->SetFillStyle(3001);
  }
  h2total->Draw("same e2");
  if(h2sys) h2->Draw("same e2");

  TLegend* legend=GetLegend(h1,hstack?(TH1*)hstack:h2,option);
  if(!option.Contains("1:noleg")) legend->Draw();

  if(option.Contains("logy")){
    gPad->SetLogy();
  }else{
    double maximum=h1total->GetMaximum()>h2total->GetMaximum()?h1total->GetMaximum():h2total->GetMaximum();
    double minimum=h1total->GetMinimum()<h2total->GetMinimum()?h1total->GetMinimum():h2total->GetMinimum();
    double range=fabs(maximum-minimum);
    h1total->GetYaxis()->SetRangeUser(minimum<0?minimum-0.1*range:0,maximum+0.1*range);
  }
  h1total->Draw("same a e");

  c1->cd(2);
  gPad->SetPad(0,0,1,0.35);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.2);
  gPad->SetGridx();gPad->SetGridy();

  TH1* ratio=(TH1*)h1->Clone("ratio");
  ratio->SetTitle("");
  ratio->SetStats(0);
  double defaultval=1.;
  if(option.Contains("diff")){
    defaultval=0.;
    ratio->Add(h2,-1);
    ratio->GetYaxis()->SetTitle("Data - Simulation");
    ratio->GetYaxis()->SetRangeUser(-0.037,0.037);
    ratio->GetYaxis()->SetLabelSize(0.06);
  }else{
    defaultval=1.;
    ratio->Divide(h2);
    ratio->GetYaxis()->SetTitle("Data/Simulation");
    ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    ratio->GetYaxis()->SetLabelSize(0.1);
  }
  ratio->GetYaxis()->SetTitleSize(0.1);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetXaxis()->SetTitle(h2->GetTitle());
  ratio->GetXaxis()->SetTitleSize(0.09);
  ratio->GetXaxis()->SetLabelSize(0.09);
  ratio->Draw();
  TLine *line=new TLine(xmin?xmin:ratio->GetXaxis()->GetXmin(),defaultval,xmax?xmax:ratio->GetXaxis()->GetXmax(),defaultval);
  line->SetLineColor(2);
  line->Draw();
  if(h1sys||h2sys){
    TH1* ratiosys=(TH1*)h2->Clone("ratiosys");
    for(int i=1;i<ratiosys->GetNbinsX()+1;i++){
      double yh1sys=h1sys?h1sys->GetBinContent(i):0;
      double eyh1sys=(yh1sys==0)?0:h1sys->GetBinError(i);
      double yh2sys=h2sys?h2sys->GetBinContent(i):0;
      double eyh2sys=yh2sys==0?0:h2sys->GetBinError(i);
      ratiosys->SetBinContent(i,defaultval);
      if(option.Contains("diff")) ratiosys->SetBinError(i,sqrt(pow(eyh1sys,2)+pow(eyh2sys,2)));  
      else ratiosys->SetBinError(i,sqrt((yh1sys?pow(eyh1sys/yh1sys,2):0)+(yh2sys?pow(eyh2sys/yh2sys,2):0)));  
    }
    ratiosys->SetFillStyle(3002);
    ratiosys->SetFillColor(samples.at(1).color);
    ratiosys->Draw("same p e2");
    TLegend* syslegend=new TLegend(0.6,0.75,0.89,0.95);
    if(option.Contains("diff")) syslegend->AddEntry(ratio,"Data - Simulation (Stat.)","lp");
    else syslegend->AddEntry(ratio,"Data/Simulation (Stat.)","lp");
    syslegend->AddEntry(ratiosys,"Syst.","f");
    if(!option.Contains("2:noleg")) syslegend->Draw();
  }
  ratio->Draw("same");

  return c1;
}
TCanvas* GetCompare(TString histname,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option=""){
  if(option.Contains("AFB")) option+=" BGSub diff leftleg";
  vector<TH1*> hists_central=GetHists(histname,histname,histname,option);
  TH1* hdata_central=GetHData(hists_central,option);
  TH1* hmc_central=GetHMC(hists_central,option);
  THStack* hstack_central=option.Contains("BGSub")?NULL:(THStack*)GetHMC(hists_central,"stack");
  vector<TH1*> hdata_syss,hmc_syss;
  for(int i=0;i<systematics.size();i++){
    if(systematics[i].type==SystematicType::MULTI) continue;
    if(sysbit&systematics[i].sysbit){
      if(DEBUG) std::cout<<"sysname="<<systematics[i].name<<" systype="<<GetStringSystematicType(systematics[i].type)<<endl;
      vector<TH1*> hdata_variations;
      vector<TH1*> hmc_variations;
      for(int j=0;j<systematics[i].suffixes.size();j++){
	TString datahistname=histname+(systematics[i].exist_data?systematics[i].suffixes[j]:"");
	TString dyhistname=histname+systematics[i].suffixes[j];
	TString bghistname=histname+(systematics[i].exist_bg?systematics[i].suffixes[j]:"");
	if(DEBUG) std::cout<<datahistname<<" "<<dyhistname<<" "<<bghistname<<endl;
	vector<TH1*> gethists=GetHists(datahistname,dyhistname,bghistname,option);
	hdata_variations.push_back(GetHData(gethists,option));
	hmc_variations.push_back(GetHMC(gethists,option));
	for(int k=0;k<(int)gethists.size();k++) delete gethists.at(k);
      }
      if(systematics[i].type==SystematicType::ENVELOPE){
	hdata_syss.push_back(GetEnvelope(hdata_central,hdata_variations));
	hmc_syss.push_back(GetEnvelope(hmc_central,hmc_variations));
      }else if(systematics[i].type==SystematicType::GAUSSIAN){
	hdata_syss.push_back(GetRMSError(hdata_central,hdata_variations));
	hmc_syss.push_back(GetRMSError(hmc_central,hmc_variations));
      }else if(systematics[i].type==SystematicType::HESSIAN){
	hdata_syss.push_back(GetHessianError(hdata_central,hdata_variations));
	hmc_syss.push_back(GetHessianError(hmc_central,hmc_variations));
      }else{
	cout<<"###WARNING### [GetCompare] Wrong SystematicType "<<systematics[i].type<<endl;
      }
      for(int j=0;j<(int)hdata_variations.size();j++){
	delete hdata_variations.at(j);
	delete hmc_variations.at(j);
      }
      if(DEBUG) std::cout<<systematics[i].name+": "<<hdata_variations.size()<<" variations"<<endl;
    }
  }
  TH1 *hdata_sys=NULL,*hmc_sys=NULL;
  if(hdata_syss.size()>0){
    hdata_sys=(TH1*)hdata_central->Clone("hdata_sys");
    hmc_sys=(TH1*)hmc_central->Clone("hmc_sys");
    for(int i=0;i<hdata_sys->GetNbinsX()+2;i++){
      hdata_sys->SetBinError(i,0);
      hmc_sys->SetBinError(i,0);
    }
    for(int i=0;i<(int)hdata_syss.size();i++){
      AddError(hdata_sys,hdata_syss.at(i));
      AddError(hmc_sys,hmc_syss.at(i));
      delete hdata_syss.at(i);delete hmc_syss.at(i);
    }
  }
  TCanvas* c1=GetCompare(hdata_central,hdata_sys,hstack_central?(TH1*)hstack_central:hmc_central,hmc_sys,rebin,xmin,xmax,option);
  return c1;
}
TH1* GetAxisParent(TPad* pad){
  TList* list=pad->GetListOfPrimitives();
  for(int i=0;i<list->GetSize();i){
    if(strstr(list->At(i)->ClassName(),"TH")!=NULL) return (TH1*)list->At(i);
  }
  return NULL;
}
TCanvas* GetCompareAFBAll(vector<TString> histnames,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option=""){
  TString canvasname=histnames[0];
  canvasname.ReplaceAll("_y0.0to0.4","");
  int nhist=histnames.size();
  TCanvas* c1=new TCanvas(canvasname,canvasname,800,800);
  c1->Divide(nhist,1);
  for(int i=0;i<nhist;i++){
    TCanvas* ctemp=GetCompare(histnames[i],sysbit,rebin,xmin,xmax,option+(i==0?" ":" 1:noleg 2:noleg"));
    TH1* hdata=GetAxisParent((TPad*)ctemp->GetPad(1));
    TH1* hratio=GetAxisParent((TPad*)ctemp->GetPad(2));
    TLegend *leg1,*leg2;
    double sf0=(0.85/nhist)/(0.1+0.85/nhist);
    double sf1=(0.85/nhist)/(0.05+0.85/nhist);
    ctemp->GetPad(1)->SetLeftMargin(i==0?1-sf0:0);
    ctemp->GetPad(2)->SetLeftMargin(i==0?1-sf0:0);
    ctemp->GetPad(1)->SetRightMargin(i==nhist-1?1-sf1:0);
    ctemp->GetPad(2)->SetRightMargin(i==nhist-1?1-sf1:0);

    hdata->GetYaxis()->SetRangeUser(-0.14,0.23);
    hdata->GetXaxis()->SetNdivisions(503);

    TString ratiotitle=hratio->GetXaxis()->GetTitle();
    ratiotitle=ratiotitle("_y.*/");
    ratiotitle=ratiotitle(1,ratiotitle.Length()-2);
    hratio->GetXaxis()->SetTitle(ratiotitle);
    hratio->GetXaxis()->CenterTitle();
    hratio->GetXaxis()->SetNdivisions(503);
    hratio->GetXaxis()->SetLabelSize(0.15);
    hratio->GetXaxis()->SetLabelOffset(-0.05);
    hratio->GetXaxis()->SetTitleOffset(0.5);
    hratio->GetXaxis()->SetTitleSize(0.2);
    if(i==0){
      leg1=(TLegend*)ctemp->GetPad(1)->GetPrimitive("TPave");
      leg1->SetX1(1-sf0+0.02);
      leg1->SetX2(0.9);
      leg1->SetTextSize(0.13);
      leg2=(TLegend*)ctemp->GetPad(2)->GetPrimitive("TPave");
      if(leg2){
	((TLegendEntry*)leg2->GetListOfPrimitives()->At(0))->SetLabel("Stat.");
	leg2->SetX1(1-sf0+0.02);
	leg2->SetX2(0.9);
	leg2->SetTextSize(0.13);
      }
      hdata->GetYaxis()->SetLabelSize(0.12);
      hdata->GetYaxis()->SetLabelOffset(0.02);
      hdata->GetYaxis()->SetTitle("A_{FB}");
      hdata->GetYaxis()->SetTitleSize(0.15);
      hdata->GetYaxis()->SetTitleOffset(1.3);

      hratio->GetYaxis()->SetLabelSize(0.1);
      hratio->GetYaxis()->SetLabelOffset(0.02);
      hratio->GetYaxis()->SetTitleOffset(1.8);
      hratio->GetYaxis()->SetTitleSize(0.12);

      hratio->GetXaxis()->SetLabelSize(hratio->GetXaxis()->GetLabelSize()*sf0);
      hratio->GetXaxis()->SetLabelOffset(0.001);
      hratio->GetXaxis()->SetTitleSize(hratio->GetXaxis()->GetTitleSize()*sf0);
      hratio->GetXaxis()->SetTitleOffset(hratio->GetXaxis()->GetTitleOffset()/sf0);
    }else if(i==nhist-1){
      hratio->GetXaxis()->SetLabelSize(hratio->GetXaxis()->GetLabelSize()*sf1);
      hratio->GetXaxis()->SetLabelOffset(0.025);
      hratio->GetXaxis()->SetTitleSize(hratio->GetXaxis()->GetTitleSize()*sf1);
      hratio->GetXaxis()->SetTitleOffset(hratio->GetXaxis()->GetTitleOffset()/sf1);
      hratio->GetXaxis()->SetLabelOffset(-0.02);
    }
    c1->cd(i+1);
    gPad->SetPad((i==0?0:0.1)+0.85*i/nhist,0.0,(i==nhist-1?0.15:0.1)+0.85*(i+1)/nhist,1);
    ctemp->GetPad(1)->SetGridx();
    ctemp->GetPad(1)->SetGridy();
    ctemp->DrawClonePad();
    delete ctemp;
  }
  c1->cd(0);
  TPad* titlepad=new TPad("titlepad","titlepad",0,0.94,1,1);
  titlepad->Draw();
  titlepad->cd();
  TPaveText* pavetitle=new TPaveText(0.1,0.1,0.9,0.9);
  pavetitle->AddText(c1->GetTitle());
  pavetitle->Draw();
  return c1;
}


void SaveAll(TString outputdir="plot"){
  int oldlevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kWarning;
  int dmax=directories.size();
  int smax=systematics.size();
  TCanvas* c=NULL;
  for(int id=0;id<dmax;id++){
    std::cout<<"mkdir -p "+outputdir+"/"+directories[id].name<<endl;
    system("mkdir -p "+outputdir+"/"+directories[id].name);
    int hmax=directories[id].histograms.size();
    for(int ih=0;ih<hmax;ih++){
      Histogram *histogram=&(directories[id].histograms[ih]);
      TString this_histname=directories[id].name+histogram->name;
      c=GetCompare(this_histname,0,histogram->rebin,histogram->xmin,histogram->xmax,histogram->option);
      c->SaveAs(outputdir+"/"+this_histname+".png");
      delete c;      
    }
    for(int is=0;is<smax;is++){
      if(DEBUG) std::cout<<"mkdir -p "+outputdir+"/"+directories[id].name+systematics[is].name<<endl;
      system("mkdir -p "+outputdir+"/"+directories[id].name+systematics[is].name);
      for(int ih=0;ih<hmax;ih++){
	Histogram *histogram=&(directories[id].histograms[ih]);
	TString this_histname=directories[id].name+histogram->name;
	c=GetCompare(this_histname,systematics[is].sysbit,histogram->rebin,histogram->xmin,histogram->xmax,histogram->option);
	c->SaveAs(outputdir+"/"+this_histname(0,this_histname.Last('/')+1)+systematics[is].name+this_histname(this_histname.Last('/'),this_histname.Length())+".png");
	delete c;
	/*
	if(systematics[is].suffixes.size()==1){
	  c=GetCompare(this_histname+"/"+systematics[is].suffixes[0],0,histogram->rebin,histogram->xmin,histogram->xmax,histogram->option);
	  c->SaveAs(outputdir+"/"+this_histname(0,this_histname.Last('/')+1)+systematics[is].name+this_histname(this_histname.Last('/'),this_histname.Length())+".png");
	}
	*/
      }
    }
  }  
  gErrorIgnoreLevel=oldlevel;
}
void SaveAFB(TString outputdir="plot"){
  int oldlevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kWarning;
  int smax=systematics.size();
  TString dirname=directories[0].name(0,directories[0].name.Index('/'))+"/summary/";
  TCanvas* c=NULL;
  std::cout<<"mkdir -p "+outputdir+"/"+dirname<<endl;
  system("mkdir -p "+outputdir+"/"+dirname);
  vector<TString> histnames;
  for(int i=0;i<directories.size();i++){
    if(directories[i].name.Contains(TRegexp("_y.*to.*/"))) 
      histnames.push_back(directories[i].name+"AFB");
  }
  c=GetCompareAFBAll(histnames,0,0,0,0,"AFB");
  c->SaveAs(outputdir+"/"+dirname+"AFB.png");
  delete c;      
  for(int is=0;is<smax;is++){
    c=GetCompareAFBAll(histnames,systematics[is].sysbit,0,0,0,"AFB");
    c->SaveAs(outputdir+"/"+dirname+"AFB_"+systematics[is].name+".png");
    delete c;
  }
  gErrorIgnoreLevel=oldlevel;
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// etc.//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void PrintKeys(TList* keys){
  for(int i=0;i<keys->GetSize();i++){
    TKey* key=(TKey*)keys->At(i);
    if(strcmp(key->GetClassName(),"TDirectoryFile")==0) PrintKeys(((TDirectoryFile*)key->ReadObj())->GetListOfKeys());   
    else cout<<key->GetMotherDir()->GetPath()<<" "<<key->GetName()<<"\n";
  }
}
void PrintHistList(int samplenum=1){
  TFile f(samples.at(samplenum).files.at(0));
  cout<<f.GetName()<<" "<<f.GetTitle()<<endl;
  PrintKeys(f.GetListOfKeys());
}
