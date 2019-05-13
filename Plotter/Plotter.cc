#include"Global.h"
#include"Sample.cc"

/////////////////////////////////////////////////////////////////////////////
///////////////////////////// struct and enum ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TString GetStringEColor(EColor color){
  switch(color){
  case kBlack: return "kBlack";
  case kRed: return "kRed";
  case kGreen: return "kGreen";
  case kBlue: return "kBlue";
  case kYellow: return "kYellow";
  case kMagenta: return "kMagenta";
  case kOrange: return "kOrange";
  default: return "###WARNING### Bad EColor";
  }
}

enum SystematicType{ENVELOPE,GAUSSIAN,HESSIAN,MULTI};
TString GetStringSystematicType(SystematicType type){
  switch(type){
  case ENVELOPE: return "ENVELOPE";
  case GAUSSIAN: return "GAUSSIAN";
  case HESSIAN: return "HESSIAN";
  case MULTI: return "MULTI";
  default: return "###ERROR### Bad SystematicType";
  }
}
struct Systematic{
  TString name;
  SystematicType type;
  vector<TString> suffixes;
  int sysbit;
  int varibit;
};
struct Plot{
  enum Type{Compare,Ratio,CompareAndRatio,CompareAndDiff};
  Type type;
  TString name;
  TString histname;
  TString sysname;
  int rebin;
  double xmin;
  double xmax;
  TString option;
  int npad;
  
};
enum Channel{MUON,ELECTRON};
TString GetStringChannel(Channel channel){
  switch(channel){
  case MUON: return "muon";
  case ELECTRON: return "electron";
  default: return "###ERROR### Bad Channel";
  }
}

/////////////////////////////////////////////////////////////////////////////
//////////////////// Add functions for global variables//////////////////////
/////////////////////////////////////////////////////////////////////////////
class Plotter{
public:
  Plotter();
  ~Plotter();
  map<TString,Sample> samples;
  map<TString,Systematic> systematics;
  map<TString,Plot> plots;

  //AddSystematic
  Systematic MakeSystematic(TString name_,SystematicType type_,vector<TString> includes,int varibit);
  Systematic MakeSystematic(TString name_,SystematicType type_,TString includes_,int varibit);

  //AddPlot
  void AddPlot(TString name_,int rebin_,double xmin_,double xmax_,TString option_="");
  
  //Hist
  TH1* GetHistRaw(TString filename,TString histname);  
  TH1* GetHistSampleFrag(const SampleFrag& frag,TString histname);
  TH1* GetHist(TString samplekey,TString histname,TString suffix="",int varibit=0);
  tuple<TH1*,TH1*> GetHistWithTotal(TString samplekey,TString histname,TString sysname);
  //  vector<TH1*> GetHists(TString datahistname,TString signalhistname="",TString bghistname="",int rebin=0,double xmin=0,double xmax=0,TString option="");
  TH1* GetHistWeightedAFB(TH1* hist_forward_num,TH1* hist_forward_den,TH1* hist_backward_num,TH1* hist_backward_den);
  TH1* GetHistAFB(TH1* hist_forward,TH1* hist_backward);
  TH1* GetHMC(const vector<TH1*>& hists,TString option="");
  static TH1* GetTH1(TH1* hstack);
  TH1* GetHData(const vector<TH1*>& hists,TString option="");
  static TLegend* GetLegend(const vector<TH1*>& hists,TString option);
  static TLegend* GetLegend(const vector<tuple<TH1*,TH1*>>& hists,TString option);
  TH1* GetEnvelope(TH1* central,const vector<TH1*>& variations);
  TH1* GetEnvelope(TH1* central,TH1* variation1,TH1* variation2=NULL,TH1* variation3=NULL,TH1* variation4=NULL,TH1* variation5=NULL,TH1* variation6=NULL,TH1* variation7=NULL,TH1* variation8=NULL,TH1* variation9=NULL);
  TH1* GetHessianError(TH1* central,const vector<TH1*>& variations);
  TH1* GetRMSError(TH1* central,const vector<TH1*>& variations);
  int AddError(TH1* hist,TH1* sys);
  bool CheckHists(vector<TH1*> hists);
  void PrintSampleFrags();
  void PrintSamples(bool detail=false);
  static TCanvas* GetCompare(vector<tuple<TH1*,TH1*>> hists,TString option);
  static TCanvas* GetRatio(vector<tuple<TH1*,TH1*>> hists,TString option);
  static TCanvas* GetCompareAndRatio(vector<tuple<TH1*,TH1*>> hists,TString option);
  TCanvas* GetCanvas(TCanvas* (*func)(vector<tuple<TH1*,TH1*>>,TString),TString histname,TString sysname,int rebin,double xmin,double xmax,TString option);
  TCanvas* GetCanvas(TString plotkey);
  static TH1* GetAxisParent(TVirtualPad* pad);
  void RebinXminXmax(TH1* hist,int rebin,double xmin,double xmax);
  /*
  TCanvas* GetCompare(TH1* h1,TH1* h1sys,TH1* h2,TH1* h2sys,TString option);
  TCanvas* GetCompare(TString datahistname,TString signalhistname,TString bghistname,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option="");
  TCanvas* GetCompare(TString histname,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option="");
  void SavePlots(TString outputdir="plot",int njob=1,int ijob=0);
  
  set<TString> GetHistKeys(TList* keys,TRegexp regexp=".*");
  set<TString> GetHistKeys(int samplenum=1,TRegexp regexp=".*");
  void PrintHistKeys(int samplenum=1,TRegexp regexp=".*");
  set<TString> ParseHistKeys(set<TString> histkeys_raw,set<TString> prefixes,set<TString> suffixes,set<TString> excludes);
  void AddPlotsAuto(set<TString> excludes);
  void PrintPlots(TRegexp reg);
  */
  //TCanvas* GetCompareAFBAll(vector<TString> histnames,int sysbit=0,TString option="");
  //TCanvas* GetCompareAFBAll(TRegexp regexp,int sysbit=0,TString option="");
  TString GetHistNameWithPrefix(TString prefix,TString histname);
};
Plotter::Plotter(){
  TH1::SetDefaultSumw2(kTRUE);
}
Plotter::~Plotter(){
}
Systematic Plotter::MakeSystematic(TString name_,SystematicType type_,vector<TString> includes,int varibit_){
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
  systematic.varibit=varibit_;
  cout<<" [AddSystematic] "<<systematic.name<<" "<<GetStringSystematicType(systematic.type)<<" sysbit:"<<systematic.sysbit<<" varibit:"<<systematic.varibit<<endl;
  cout<<"  INCLUDE=";for(int i=0;i<includes.size();i++) cout<<includes[i]<<" ";
  cout<<endl;
  return systematic;
}
Systematic Plotter::MakeSystematic(TString name_,SystematicType type_,TString includes_,int varibit_){
  TObjArray* arr=includes_.Tokenize(" ");
  vector<TString> includes;
  for(int i=0;i<arr->GetEntries();i++){
    includes.push_back(((TObjString*)arr->At(i))->String());
  }
  return MakeSystematic(name_,type_,includes,varibit_);
}
/*
void Plotter::AddPlot(TString name_,int rebin_,double xmin_,double xmax_,TString option_=""){
  Plot plot;
  plot.name=name_;
  plot.rebin=rebin_;
  plot.xmin=xmin_;
  plot.xmax=xmax_;
  plot.option=option_;
  if(DEBUG) std::cout<<" [AddPlot] to "<<plot.name<<" "<<plot.rebin<<" "<<plot.xmin<<" "<<plot.xmax<<endl;
  plots["/"+name_]=plot;
}
*/
/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Core functions////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TH1* Plotter::GetHistRaw(TString filename,TString histname){
  TFile *f=new TFile(filename);
  TH1* hist=(TH1*)f->Get(histname);
  if(hist){
    hist->SetDirectory(0);
    if(DEBUG>1) cout<<"###INFO### [GetHist] get "<<histname<<" in "<<filename<<endl;
  }else{
    if(DEBUG>0) cout<<"###WARNING### [GetHist] no "<<histname<<" in "<<filename<<endl;
  }
  f->Close();
  delete f;
  return hist;
}
TH1* Plotter::GetHistSampleFrag(const SampleFrag& frag,TString histname){
  TH1* hist=NULL;
  for(unsigned int i=0;i<frag.files.size();i++){
    TString filepath=get<0>(frag.files[i]);
    TString histnamewithprefix=GetHistNameWithPrefix(get<1>(frag.files[i]),histname);
    double fileweight=get<2>(frag.files[i]);
    TH1* this_hist=GetHistRaw(filepath,histnamewithprefix);
    if(!hist){
      hist=(TH1*)this_hist->Clone();
      hist->Reset();
    }
    hist->Add(this_hist,fileweight);
    delete this_hist;
  }
  frag.ApplyHistStyle(hist);
  return hist;
} 
TH1* Plotter::GetHist(TString samplekey,TString histname,TString suffix,int varibit){
  TH1* hist=NULL;
  Sample* this_sample=&samples[samplekey];
  if(histname.Contains("weightedAFB")){
    TString histname_forward_num=histname; histname_forward_num.ReplaceAll("weightedAFB","forward_num");
    TString histname_forward_den=histname; histname_forward_num.ReplaceAll("weightedAFB","forward_den");
    TString histname_backward_num=histname; histname_backward_num.ReplaceAll("weightedAFB","backward_num");
    TString histname_backward_den=histname; histname_backward_num.ReplaceAll("weightedAFB","backward_den");
    TH1* hist_forward_num=GetHist(samplekey,histname_forward_num,suffix,varibit);
    TH1* hist_forward_den=GetHist(samplekey,histname_forward_den,suffix,varibit);
    TH1* hist_backward_num=GetHist(samplekey,histname_backward_num,suffix,varibit);
    TH1* hist_backward_den=GetHist(samplekey,histname_backward_den,suffix,varibit);
    hist=GetHistWeightedAFB(hist_forward_num,hist_forward_den,hist_backward_num,hist_backward_den);
    delete hist_forward_num; delete hist_forward_den; delete hist_backward_num; delete hist_backward_den;
    return hist;
  }else if(histname.Contains("AFB")){
    TString histname_forward=histname; histname_forward.ReplaceAll("AFB","forward");
    TString histname_backward=histname; histname_backward.ReplaceAll("AFB","backward");
    TH1* hist_forward=GetHist(samplekey,histname_forward,suffix,varibit);
    TH1* hist_backward=GetHist(samplekey,histname_backward,suffix,varibit);
    hist=GetHistAFB(hist_forward,hist_backward);
    delete hist_forward; delete hist_backward;
    return hist;
  }else{
    for(int i=this_sample->frags.size()-1;i>=0;i--){
      SampleFrag *frag=&get<0>(this_sample->frags[i]);
      double fragweight=get<1>(this_sample->frags[i]);
      TString histnamewithsuffix=histname;
      if((1<<frag->type)&varibit) histnamewithsuffix+=suffix;
      TH1* this_hist=GetHistSampleFrag(*frag,histnamewithsuffix);
      if(this_sample->type==Sample::Type::SUM){
	if(!hist){
	  hist=(TH1*)this_hist->Clone();
	  hist->Reset();
	}
	hist->Add(this_hist,fragweight);
	delete this_hist;
      }else if(this_sample->type==Sample::Type::STACK){
	if(!hist){
	  hist=(TH1*)new THStack(this_sample->title,this_sample->title);
	}
	THStack* hstack=(THStack*)hist;
	if(fragweight!=1.) if(DEBUG>0) cout<<"###WARNING### [Plotter::GetHist] STACK with weight!=1. is not possible yet (weight="<<fragweight<<")"<<endl;
	hstack->Add(this_hist,"HIST");
      }else{
	if(DEBUG>0) cout<<"###WARNING### [Plotter::GetHist] wrong Sample::Type "<<this_sample->type<<endl;
      }
    }
    this_sample->ApplyHistStyle(hist);
    hist->SetTitle(histname+suffix);
    return hist;
  }
}
tuple<TH1*,TH1*> Plotter::GetHistWithTotal(TString samplekey,TString histname,TString sysname){
  TH1* central=GetHist(samplekey,histname);
  THStack* hstack=NULL;
  if(strstr(central->ClassName(),"THStack")!=NULL){
    hstack=(THStack*)central;
    central=GetTH1((TH1*)hstack);
  }
  int sysbit=systematics[sysname].sysbit;
  vector<TH1*> syss;
  for(const auto& sysmap:systematics){
    const Systematic& systematic=sysmap.second;
    if(systematic.type==SystematicType::MULTI) continue;
    if(sysbit&systematic.sysbit){
      if(DEBUG>1) std::cout<<"###INFO### [Plotter::GetHistWithSys] sysname="<<systematic.name<<" systype="<<GetStringSystematicType(systematic.type)<<endl;
      vector<TH1*> variations;
      for(const auto& suffix:systematic.suffixes){
	TH1* this_hist=GetHist(samplekey,histname,suffix,systematic.varibit);
	variations.push_back(GetTH1(this_hist));
	delete this_hist;
      }
      if(systematic.type==SystematicType::ENVELOPE){
	syss.push_back(GetEnvelope(central,variations));
      }else if(systematic.type==SystematicType::GAUSSIAN){
	syss.push_back(GetRMSError(central,variations));
      }else if(systematic.type==SystematicType::HESSIAN){
	syss.push_back(GetHessianError(central,variations));
      }else{
	cout<<"###ERROR### [GetCompare] Wrong SystematicType "<<systematic.type<<endl;
      }
      if(DEBUG>1) std::cout<<"###INFO### [Plotter::GetHistWithSys] "<<systematic.name+": "<<variations.size()<<" variations"<<endl;
      for(unsigned int j=0;j<variations.size();j++){
	delete variations.at(j);
      }
    }
  }
  TH1 *total=NULL;
  if(syss.size()>0){
    total=(TH1*)central->Clone();
    for(int i=0;i<(int)syss.size();i++){
      AddError(total,syss.at(i));
      delete syss.at(i);
    }
  }
  if(hstack){
    delete central;
    central=(TH1*)hstack;
  }
  return make_tuple(central,total);
}

TH1* Plotter::GetHistWeightedAFB(TH1* hist_forward_num,TH1* hist_forward_den,TH1* hist_backward_num,TH1* hist_backward_den){
  TH1* hist=(TH1*)hist_forward_num->Clone();
  hist->SetBit(kCanDelete);
  hist->Reset();
  for(int i=0;i<hist->GetNbinsX()+2;i++){
    double valfn=hist_forward_num->GetBinContent(i);
    double valbn=hist_backward_num->GetBinContent(i);
    double valfd=hist_forward_den->GetBinContent(i);
    double efd=hist_forward_den->GetBinError(i);
    double valbd=hist_backward_den->GetBinContent(i);
    double ebd=hist_backward_den->GetBinError(i);
    hist->SetBinContent(i,3./8.*(valfd-valbd)/(valfn+valbn));
    hist->SetBinError(i,3./8.*(valbn*valfd+valfn*valbd)/pow(valfn+valbn,2)*sqrt(pow(efd/valfd,2)+pow(ebd/valbd,2)));
  }
  TString this_title=hist_forward_num->GetTitle();
  hist->SetTitle(this_title.ReplaceAll("forward_num","weightedAFB"));
  return hist;
}
TH1* Plotter::GetHistAFB(TH1* hist_forward,TH1* hist_backward){
  if(!hist_forward||!hist_backward){
    return NULL;
  }
  TH1* hist=(TH1*)hist_forward->Clone();
  hist->SetBit(kCanDelete);
  hist->Reset();
  for(int i=0;i<hist->GetNbinsX()+2;i++){
    double valf=hist_forward->GetBinContent(i);
    double ef=hist_forward->GetBinError(i);
    double valb=hist_backward->GetBinContent(i);
    double eb=hist_backward->GetBinError(i);
    hist->SetBinContent(i,(valf-valb)/(valf+valb));
    hist->SetBinError(i,2*sqrt(ef*ef*valb*valb+eb*eb*valf*valf)/pow(valf+valb,2));
  }
  TString this_title=hist_forward->GetTitle();
  hist->SetTitle(this_title.ReplaceAll("forward","AFB"));
  return hist;
}
/*
TH1* Plotter::GetHMC(const vector<TH1*>& hists,TString option=""){
  TH1* hist=NULL;
  if(option.Contains("weightedAFB")){
    vector<TH1*> hists_forward_num(hists.begin(),hists.begin()+samples.size());
    vector<TH1*> hists_backward_num(hists.begin()+samples.size(),hists.end()+samples.size()*2);
    vector<TH1*> hists_forward_den(hists.begin()+samples.size()*2,hists.begin()+samples.size()*3);
    vector<TH1*> hists_backward_den(hists.begin()+samples.size()*3,hists.begin()+samples.size()*4);
    TH1* hist_forward_num=GetHMC(hists_forward_num,"BGSub");
    TH1* hist_backward_num=GetHMC(hists_backward_num,"BGSub");
    TH1* hist_forward_den=GetHMC(hists_forward_den,"BGSub");
    TH1* hist_backward_den=GetHMC(hists_backward_den,"BGSub");
    hist=GetHistWeightedAFB(hist_forward_num,hist_forward_den,hist_backward_num,hist_backward_den);
    delete hist_forward_num;delete hist_backward_num;delete hist_forward_den;delete hist_backward_den;
  }else if(option.Contains("AFB")){
    vector<TH1*> hists_forward(hists.begin(),hists.begin()+samples.size());
    vector<TH1*> hists_backward(hists.begin()+samples.size(),hists.begin()+samples.size()*2);
    TH1* hist_forward=GetHMC(hists_forward,"BGSub");
    TH1* hist_backward=GetHMC(hists_backward,"BGSub");
    hist=GetHistAFB(hist_forward,hist_backward);
    delete hist_forward;delete hist_backward;
  }else if(option.Contains("stack")){
    THStack* hstack=new THStack("hstack","hstack");
    hstack->SetBit(kCanDelete);
    for(int i=(int)samples.size()-1;i>=0;i--){
      if(samples[i].type==SampleFrag::Type::DATA) continue;
      if(hists.at(i)){
	TH1* this_hist=(TH1*)hists.at(i)->Clone();
	if(samples.at(i).type==SampleFrag::Type::SIGNAL) this_hist->SetFillColor(samples.at(i).color-9);
	this_hist->SetBit(kCanDelete);
	hstack->Add(this_hist,"HIST");
      }
    }
    return (TH1*)hstack;
  }else if(option.Contains("BGSub")){
    for(int i=0;i<(int)samples.size();i++){
      if(samples[i].type==SampleFrag::Type::SIGNAL&&hists.at(i)){
	if(hist) hist->Add(hists.at(i));
	else{
	  hist=(TH1*)hists.at(i)->Clone("hmc");
	  hist->SetBit(kCanDelete);
	}
      }
    }
  }else{
    for(int i=0;i<(int)samples.size();i++){
      if(samples[i].type!=SampleFrag::Type::DATA&&hists.at(i)){
	if(hist) hist->Add(hists.at(i));
	else{
	  hist=(TH1*)hists.at(i)->Clone("hmc");
	  hist->SetBit(kCanDelete);
	}
      }
    }
  }
  return hist;
}
*/
TH1* Plotter::GetTH1(TH1* hstack){
  TH1* hist=NULL;
  if(strstr(hstack->ClassName(),"THStack")!=NULL){
    hist=(TH1*)((THStack*)hstack)->GetHists()->Last()->Clone();
    hist->Reset();
    for(const auto& this_hist:*((THStack*)hstack)->GetHists()){
      hist->Add((TH1*)this_hist);
    }
  }else{
    hist=(TH1*)hstack->Clone();
  }
  return hist;
}
TH1* Plotter::GetHData(const vector<TH1*>& hists,TString option=""){
  TH1* hist=NULL;
  if(option.Contains("weightedAFB")){
    vector<TH1*> hists_forward_num(hists.begin(),hists.begin()+samples.size());
    vector<TH1*> hists_backward_num(hists.begin()+samples.size(),hists.end()+samples.size()*2);
    vector<TH1*> hists_forward_den(hists.begin()+samples.size()*2,hists.begin()+samples.size()*3);
    vector<TH1*> hists_backward_den(hists.begin()+samples.size()*3,hists.begin()+samples.size()*4);
    TH1* hist_forward_num=GetHData(hists_forward_num,"BGSub");
    TH1* hist_backward_num=GetHData(hists_backward_num,"BGSub");
    TH1* hist_forward_den=GetHData(hists_forward_den,"BGSub");
    TH1* hist_backward_den=GetHData(hists_backward_den,"BGSub");
    hist=GetHistWeightedAFB(hist_forward_num,hist_forward_den,hist_backward_num,hist_backward_den);
    delete hist_forward_num;delete hist_backward_num;delete hist_forward_den;delete hist_backward_den;
  }else if(option.Contains("AFB")){
    vector<TH1*> hists_forward(hists.begin(),hists.begin()+samples.size());
    vector<TH1*> hists_backward(hists.begin()+samples.size(),hists.begin()+samples.size()*2);
    TH1* hist_forward=GetHData(hists_forward,"BGSub");
    TH1* hist_backward=GetHData(hists_backward,"BGSub");
    hist=GetHistAFB(hist_forward,hist_backward);
    delete hist_forward;delete hist_backward;
  }else if(option.Contains("BGSub")){
    hist=GetHData(hists,"");
    for(int i=0;i<(int)samples.size();i++) if(samples[i].type==SampleFrag::Type::BG&&hists.at(i)) hist->Add(hists.at(i),-1.);
  }else{
    for(int i=0;i<(int)samples.size();i++){
      if(samples[i].type==SampleFrag::Type::DATA&&hists.at(i)){
	if(hist) hist->Add(hists.at(i));
	else{
	  hist=(TH1*)hists.at(i)->Clone();
	  hist->SetBit(kCanDelete);
	}
      }
    }
  }
  return hist;
}
TLegend* Plotter::GetLegend(const vector<TH1*>& hists,TString option){
  double horizontalshift=0;
  if(option.Contains("leftleg")) horizontalshift=-0.53;
  TLegend* legend=new TLegend(0.67+horizontalshift,0.88-0.07,0.89+horizontalshift,0.88);
  for(const auto& hist:hists){
    if(strstr(hist->ClassName(),"THStack")){
      TIter next(((THStack*)hist)->GetHists(),kIterBackward);
      while(TObject *obj = next()){
	legend->AddEntry(obj,obj->GetName());
      }
    }else legend->AddEntry(hist,hist->GetName());
  }
  legend->SetY1(0.88-legend->GetNRows()*0.07);
  legend->SetBorderSize(0);
  return legend;
}
TLegend* Plotter::GetLegend(const vector<tuple<TH1*,TH1*>>& hists,TString option){
  vector<TH1*> hists_final;
  for(const auto& histtuple : hists){
    hists_final.push_back(get<0>(histtuple));
  }
  return GetLegend(hists_final,option);
}
TH1* Plotter::GetEnvelope(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetBit(kCanDelete);
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,0);
  for(int i=0;i<(int)variations.size();i++){
    for(int j=0;j<syshist->GetNbinsX()+1;j++){
      double diff=fabs(syshist->GetBinContent(j)-variations.at(i)->GetBinContent(j));
      if(diff>syshist->GetBinError(j)) syshist->SetBinError(j,diff);
    }
  }
  return syshist;
}
TH1* Plotter::GetEnvelope(TH1* central,TH1* variation1,TH1* variation2=NULL,TH1* variation3=NULL,TH1* variation4=NULL,TH1* variation5=NULL,TH1* variation6=NULL,TH1* variation7=NULL,TH1* variation8=NULL,TH1* variation9=NULL){
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
TH1* Plotter::GetHessianError(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetBit(kCanDelete);
  for(int i=1;i<syshist->GetNbinsX()+1;i++) syshist->SetBinError(i,0);
  for(int i=0;i<(int)variations.size();i++){
    for(int j=0;j<syshist->GetNbinsX()+1;j++){
      double diff=fabs(syshist->GetBinContent(j)-variations.at(i)->GetBinContent(j));
      syshist->SetBinError(j,sqrt(pow(syshist->GetBinError(j),2)+pow(diff,2)));
    }
  }
  return syshist;
}  
TH1* Plotter::GetRMSError(TH1* central,const vector<TH1*>& variations){
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetBit(kCanDelete);
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
int Plotter::AddError(TH1* hist,TH1* sys){
  for(int i=1;i<hist->GetNbinsX()+1;i++){
    if(fabs(hist->GetBinContent(i)-sys->GetBinContent(i))*1000000>fabs(hist->GetBinContent(i))){
      cout<<"###ERROR### [AddError] systematic hist is wrong"<<endl;
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
TCanvas* Plotter::GetCompare(vector<tuple<TH1*,TH1*>> hists,TString option){
  if(hists.size()==0){
    cout<<"zero length"<<endl;
    return NULL;
  }
  TH1* axisowner=get<1>(hists[0]);
  if(!axisowner) axisowner=get<0>(hists[0]);

  TCanvas* c=new TCanvas;
  axisowner->SetStats(0);
  axisowner->Draw();
  for(unsigned int i=0;i<hists.size();i++){
    if(strstr(get<0>(hists[i])->ClassName(),"THStack"))
      get<0>(hists[i])->Draw("same e");
  }
  for(unsigned int i=0;i<hists.size();i++){
    if(get<1>(hists[i])){
      get<1>(hists[i])->Draw("same e");
    }
    if(strstr(get<0>(hists[i])->ClassName(),"THStack")==NULL)
      get<0>(hists[i])->Draw("same e");
  }
  /*
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
*/
  TLegend* legend=GetLegend(hists,option);
  if(!option.Contains("1:noleg")) legend->Draw();

  if(option.Contains("logy")){
    c->SetLogy();
  }else{
    //double maximum=h1total->GetMaximum()>h2total->GetMaximum()?h1total->GetMaximum():h2total->GetMaximum();
    //double minimum=h1total->GetMinimum()<h2total->GetMinimum()?h1total->GetMinimum():h2total->GetMinimum();
    //double range=fabs(maximum-minimum);
    //h1total->GetYaxis()->SetRangeUser(minimum<0?minimum-0.1*range:0,maximum+0.1*range);
  }
  //h1total->Draw("same a e");
  return c;
}
TCanvas* Plotter::GetRatio(vector<tuple<TH1*,TH1*>> hists,TString option){
  vector<tuple<TH1*,TH1*>> hists_new;
  TH1* base=NULL;
  if(strstr(get<0>(hists[0])->GetName(),"data")&&!strstr(get<0>(hists[1])->GetName(),"data"))
    base=GetTH1(get<0>(hists[1]));
  else 
    base=GetTH1(get<0>(hists[0]));
  
  for(int i=0;i<base->GetNbinsX()+2;i++){
    base->SetBinError(i,0);
  }
  for(const auto& histtuple:hists){
    TH1 *hist0=NULL,*hist1=NULL;
    if(get<0>(histtuple)){
      if(strstr(get<0>(histtuple)->ClassName(),"THStack")) hist0=GetTH1(get<0>(histtuple));
      else hist0=(TH1*)get<0>(histtuple)->Clone();
      hist0->Divide(base);
    }
    if(get<1>(histtuple)){
      hist1=(TH1*)get<1>(histtuple)->Clone();
      hist1->Divide(base);
    }
    hists_new.push_back(make_tuple(hist0,hist1));
  }  
  delete base;
  TCanvas* c=GetCompare(hists_new,option);
  TH1* axisowner=GetAxisParent(c);
  if(option.Contains("widewidey")){
    axisowner->GetYaxis()->SetRangeUser(0.01,1.99);
    axisowner->GetYaxis()->SetNdivisions(506);
  }else if(option.Contains("widey")){
    axisowner->GetYaxis()->SetRangeUser(0.501,1.499);
    axisowner->GetYaxis()->SetNdivisions(506);
  }else{
    axisowner->GetYaxis()->SetRangeUser(0.801,1.199);
    axisowner->GetYaxis()->SetNdivisions(504);
  }
  return c;
}
TCanvas* Plotter::GetCompareAndRatio(vector<tuple<TH1*,TH1*>> hists,TString option){
  TCanvas* c=new TCanvas;
  c->Divide(1,2);
  TCanvas* c1temp=GetCompare(hists,option);
  c->cd(1);
  c1temp->DrawClonePad();
  delete c1temp;
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);
  TH1* axisparent=GetAxisParent(gPad);
  axisparent->GetXaxis()->SetLabelSize(0);
  axisparent->GetXaxis()->SetTitle("");


  TCanvas* c2temp=GetRatio(hists,option);
  c->cd(2);
  c2temp->DrawClonePad();
  delete c2temp;
  gPad->SetPad(0,0,1,0.365);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.2);
  gPad->SetGridx();gPad->SetGridy();
  gPad->SetFillStyle(0);

  axisparent=GetAxisParent(gPad);
  axisparent->SetTitle("");
  axisparent->SetStats(0);
  axisparent->GetYaxis()->SetLabelSize(0.1);
  /*
  double defaultval=1.;
  if(option.Contains("diff")){
    defaultval=0.;
    ratio1->Add(h2,-1);
    ratio2->GetYaxis()->SetTitle("Data - Simulation");
    ratio2->GetYaxis()->SetRangeUser(-0.067,0.067);
    ratio2->GetYaxis()->SetLabelSize(0.06);
  }else{
    ratio2->GetYaxis()->SetTitle("Data/Simulation");
    ratio2->GetYaxis()->SetLabelSize(0.1);
  }
  */
  axisparent->GetYaxis()->SetTitleSize(0.1);
  axisparent->GetYaxis()->SetTitleOffset(0.5);
  axisparent->GetXaxis()->SetTitle(axisparent->GetTitle());
  axisparent->GetXaxis()->SetTitleSize(0.09);
  axisparent->GetXaxis()->SetLabelSize(0.09);
  /*
  if(h1sys||h2sys){
    TH1* ratiosys=(TH1*)h2->Clone("ratiosys");
    ratiosys->SetBit(kCanDelete);
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
    TLegend* syslegend=new TLegend(0.6,0.75,0.89,0.95);
    if(option.Contains("diff")) syslegend->AddEntry(ratio1,"Data - Simulation (Stat.)","lp");
    else syslegend->AddEntry(ratio1,"Data/Simulation (Stat.)","lp");
    syslegend->AddEntry(ratiosys,"Syst.","f");
    if(!option.Contains("2:noleg")) syslegend->Draw();
    ratiosys->Draw("same e2");
  }
  */
  return c;
}
TCanvas* Plotter::GetCanvas(TCanvas* (*func)(vector<tuple<TH1*,TH1*>>,TString),TString histname,TString sysname,int rebin,double xmin,double xmax,TString option){
  vector<tuple<TH1*,TH1*>> hists;
  for(const auto& samplemap:samples){
    hists.push_back(GetHistWithTotal(samplemap.first,histname,sysname));
    RebinXminXmax(get<0>(hists.back()),rebin,xmin,xmax);
    RebinXminXmax(get<1>(hists.back()),rebin,xmin,xmax);
  }
  return func(hists,option);
}
void Plotter::RebinXminXmax(TH1* hist,int rebin,double xmin,double xmax){
  if(hist){
    if(strstr(hist->ClassName(),"THStack")){
      for(const auto& obj:*((THStack*)hist)->GetHists())
	RebinXminXmax((TH1*)obj,rebin,xmin,xmax);
    }else{
      if(rebin) hist->Rebin(rebin);
      if(xmin||xmax) hist->GetXaxis()->SetRangeUser(xmin,xmax);
    }
  }
}
TCanvas* Plotter::GetCanvas(TString plotkey){
  Plot& plot=plots[plotkey];
  TCanvas* c=new TCanvas(plot.name,plot.name);
  c->Divide(plot.npad);
  for(int i=1;i<=plot.npad;i++){
    c->cd(i);
    if(plot.type==Plot::Type::CompareAndRatio) GetCanvas(GetCompareAndRatio,plot.histname,plot.sysname,plot.rebin,plot.xmin,plot.xmax,plot.option)->Draw();
  }
  return c;
}  
bool Plotter::CheckHists(vector<TH1*> hists){
  bool flag_data=false,flag_signal=false;
  for(unsigned int i=0;i<samples.size();i++){
    if(!hists.at(i)) continue;
    if(samples[i].type==SampleFrag::Type::DATA) flag_data=true;
    if(samples[i].type==SampleFrag::Type::SIGNAL) flag_signal=true;
  }
  if(flag_data&&flag_signal) return true;
  else return false;
}
/*
TCanvas* Plotter::GetCompare(TString datahistname,TString signalhistname,TString bghistname,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option=""){
  if(option.Contains("AFB")) option+=" BGSub diff leftleg";
  vector<TH1*> hists_central=GetHists(datahistname,signalhistname,bghistname,rebin,xmin,xmax,option);
  if(!CheckHists(hists_central)){
    if(DEBUG) cout<<"###WARNING### [GetCompare] Not enough hists for compare ("<<datahistname<<","<<signalhistname<<","<<bghistname<<")"<<endl;
    return NULL;
  }
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
	TString datahistnamesys=datahistname+(systematics[i].vary_data?systematics[i].suffixes[j]:"");
	TString signalhistnamesys=signalhistname+(systematics[i].vary_signal?systematics[i].suffixes[j]:"");
	TString bghistnamesys=bghistname+(systematics[i].vary_bg?systematics[i].suffixes[j]:"");
	if(DEBUG) std::cout<<datahistnamesys<<" "<<signalhistnamesys<<" "<<bghistnamesys<<endl;
	vector<TH1*> gethists=GetHists(datahistnamesys,signalhistnamesys,bghistnamesys,rebin,xmin,xmax,option);
	if(!CheckHists(gethists)){
	  if(DEBUG) cout<<"###WARNING### [GetCompare] Not enough hists for compare ("<<datahistnamesys<<","<<signalhistnamesys<<","<<bghistnamesys<<")"<<endl;
	  return NULL;
	}
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
	cout<<"###ERROR### [GetCompare] Wrong SystematicType "<<systematics[i].type<<endl;
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
    hdata_sys->SetBit(kCanDelete);
    hmc_sys=(TH1*)hmc_central->Clone("hmc_sys");
    hmc_sys->SetBit(kCanDelete);
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
  for(int i=0;i<(int)hists_central.size();i++) delete hists_central.at(i);
  return GetCompare(hdata_central,hdata_sys,hstack_central?(TH1*)hstack_central:hmc_central,hmc_sys,option);
}
TCanvas* Plotter::GetCompare(TString histname,int sysbit=0,int rebin=0,double xmin=0,double xmax=0,TString option=""){
  return GetCompare(histname,histname,histname,sysbit,rebin,xmin,xmax,option);
}
*/
TH1* Plotter::GetAxisParent(TVirtualPad* pad){
  TList* list=pad->GetListOfPrimitives();
  for(int i=0;i<list->GetSize();i++){
    if(strstr(list->At(i)->ClassName(),"TH")!=NULL) return (TH1*)list->At(i);
  }
  return NULL;
}
/*
TCanvas* Plotter::GetCompareAFBAll(vector<TString> histnames,int sysbit=0,TString option=""){
  TString canvasname=histnames[0];
  canvasname.ReplaceAll("_y0.0to0.4","");
  int nhist=histnames.size();
  TCanvas* c1=new TCanvas(canvasname,canvasname,800,800);
  c1->Divide(nhist,1);
  for(int i=0;i<nhist;i++){
    TCanvas* ctemp=GetCompare(histnames[i],sysbit,0,0.,0.,option+(i==0?" ":" 1:noleg 2:noleg"));
    TH1* hdata=GetAxisParent(ctemp->GetPad(1));
    TH1* hratio=GetAxisParent(ctemp->GetPad(2));
    TLegend *leg1,*leg2;
    double sf0=(0.85/nhist)/(0.1+0.85/nhist);
    double sf1=(0.85/nhist)/(0.05+0.85/nhist);
    ctemp->GetPad(1)->SetLeftMargin(i==0?1-sf0:0);
    ctemp->GetPad(2)->SetLeftMargin(i==0?1-sf0:0);
    ctemp->GetPad(1)->SetRightMargin(i==nhist-1?1-sf1:0);
    ctemp->GetPad(2)->SetRightMargin(i==nhist-1?1-sf1:0);

    hdata->GetYaxis()->SetRangeUser(-0.19,0.33);
    hdata->GetXaxis()->SetNdivisions(503);

    TString ratiotitle=hratio->GetXaxis()->GetTitle();
*/
//    ratiotitle=ratiotitle("_y.*/");
/*
    ratiotitle=ratiotitle(1,ratiotitle.Length()-2);
    hratio->GetXaxis()->SetTitle(ratiotitle);
    hratio->GetXaxis()->CenterTitle();
    hratio->GetXaxis()->SetNdivisions(503);
    hratio->GetXaxis()->SetLabelSize(0.15);
    hratio->GetXaxis()->SetLabelOffset(-0.05);
    hratio->GetXaxis()->SetTitleOffset(0.5);
    hratio->GetXaxis()->SetTitleSize(0.15);
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

TCanvas* Plotter::GetCompareAFBAll(TRegexp regexp,int sysbit=0,TString option=""){
  vector<TString> histnames;
  for(auto it=plots.begin();it!=plots.end();it++){
    if(it->first.Contains(regexp)){
      histnames.push_back(it->second.name);
    }
  }
  return GetCompareAFBAll(histnames,sysbit,option);
}


void Plotter::SavePlots(TString outputdir="plot",int njob=1,int ijob=0){
  int oldlevel=gErrorIgnoreLevel;
  if(!DEBUG) gErrorIgnoreLevel=kWarning;

  set<TString> dirs;
  int smax=systematics.size();
  TCanvas* c=NULL;
  int nplot=(plots.size()+1)/njob;
  for(auto ip=next(plots.begin(),ijob*nplot);ip!=next(plots.begin(),(ijob+1)*nplot)&&ip!=plots.end();ip++){
    Plot *plot=&ip->second;
    TString histname=plot->name;
    TString dir=outputdir+"/"+histname(0,histname.Last('/'));
    if(dirs.find(dir)==dirs.end()){
      std::cout<<"mkdir -p "+dir<<endl;
      system("mkdir -p "+dir);
      dirs.insert(dir);
    }
    if(DEBUG>1) cout<<"###INFO### [SavePlots] save "<<outputdir+"/"+histname+".png"<<endl;
    c=GetCompare(histname,0,plot->rebin,plot->xmin,plot->xmax,plot->option);
    c->SaveAs(outputdir+"/"+histname+".png");
    THStack* hstack=(THStack*)c->GetPad(1)->GetPrimitive("hstack");
    if(hstack) hstack->GetHists()->Clear();
    delete c;
    for(int is=0;is<smax;is++){
      TString this_dir=dir+"/"+systematics[is].name;
      if(dirs.find(this_dir)==dirs.end()){
	if(DEBUG) std::cout<<"mkdir -p "+this_dir<<endl;
	system("mkdir -p "+this_dir);
	dirs.insert(this_dir);
      }
      if(DEBUG>1) cout<<"###INFO### [SavePlots] save "<<this_dir+"/"+histname(histname.Last('/'),histname.Length())+".png"<<endl;
      c=GetCompare(histname,systematics[is].sysbit,plot->rebin,plot->xmin,plot->xmax,plot->option);
      if(c){
	c->SaveAs(this_dir+"/"+histname(histname.Last('/'),histname.Length())+".png");
	THStack* hstack=(THStack*)c->GetPad(1)->GetPrimitive("hstack");
	if(hstack) hstack->GetHists()->Clear();
	delete c;
      }
      if(systematics[is].suffixes.size()==1){
	if(DEBUG>1) cout<<"###INFO### [SavePlots] save "<<this_dir+"/"+histname(histname.Last('/'),histname.Length())+"_raw.png"<<endl;
	c=GetCompare(histname+(systematics[is].vary_data?systematics[is].suffixes[0]:""),histname+(systematics[is].vary_signal?systematics[is].suffixes[0]:""),histname+(systematics[is].vary_bg?systematics[is].suffixes[0]:""),0,plot->rebin,plot->xmin,plot->xmax,plot->option);
	if(c){
	  c->SaveAs(this_dir+"/"+histname(histname.Last('/'),histname.Length())+"_raw.png");
	  THStack* hstack=(THStack*)c->GetPad(1)->GetPrimitive("hstack");
	  if(hstack) hstack->GetHists()->Clear();
	  delete c;
	}
      }
    }
  }  
  gErrorIgnoreLevel=oldlevel;
}
*/
/////////////////////////////////////////////////////////////////////////////
////////////////////////////// etc.//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*
set<TString> Plotter::GetHistKeys(TList* keys,TRegexp regexp=".*"){
  set<TString> histkeys;
  for(int i=0;i<keys->GetSize();i++){
    TKey* key=(TKey*)keys->At(i);
    if(strcmp(key->GetClassName(),"TDirectoryFile")==0){
      set<TString> this_histkeys=GetHistKeys(((TDirectoryFile*)key->ReadObj())->GetListOfKeys(),regexp);
      histkeys.insert(this_histkeys.begin(),this_histkeys.end());
    }else{
      TString path=key->GetMotherDir()->GetPath();
      path=path(path.Index(":")+1,path.Length())+"/"+key->GetName();
      if(path.Contains(regexp)) histkeys.insert(path);
    }
  }
  return histkeys;
}

set<TString> Plotter::GetHistKeys(int samplenum=1,TRegexp regexp=".*"){
  TFile f(samples.at(samplenum).files.at(0));
  cout<<f.GetName()<<"/"<<f.GetTitle()<<endl;
  TList* keys=f.GetListOfKeys();
  set<TString> histkeys=GetHistKeys(keys,regexp);
  f.Close();
  return histkeys;
}
void Plotter::PrintHistKeys(int samplenum=1,TRegexp regexp=".*"){
  set<TString> histkeys=GetHistKeys(samplenum,regexp);
  for(auto it=histkeys.begin();it!=histkeys.end();it++){
    cout<<*it<<endl;
  }
}
set<TString> Plotter::ParseHistKeys(set<TString> histkeys_raw,set<TString> prefixes,set<TString> suffixes,set<TString> excludes){
  set<TString> histkeys;
  for(auto it=histkeys_raw.begin();it!=histkeys_raw.end();it++){
    bool next_flag=false;
    TString histkey=*it;
    for(auto ite=excludes.begin();!next_flag&&ite!=excludes.end();ite++){
      if(histkey.Contains(*ite)){
	next_flag=true;
	break;
      }
    }
    for(auto its=suffixes.begin();!next_flag&&its!=suffixes.end();its++){
      if(histkey.Contains(TRegexp(*its+"$"))){
	next_flag=true;
	break;
      }
    }
    if(next_flag) continue;
    for(auto itp=prefixes.begin();itp!=prefixes.end();itp++){
      int extent,preindex=histkey.Index(TRegexp("/"+*itp),&extent);
      if(preindex!=-1) histkey.Remove(preindex+1,extent-1);
    }
    histkey.ReplaceAll("forward","AFB");
    histkey.ReplaceAll("backward","AFB");
    histkey.ReplaceAll("AFB_num","weightedAFB");
    histkey.ReplaceAll("AFB_den","weightedAFB");
    histkeys.insert(histkey);
  }
  return histkeys;
}
void Plotter::AddPlotsAuto(set<TString> excludes={"^/electron..../","/qqbar","/qbarq","/gq","/qg","/gqbar","/qbarg","/qq","/gg"}){
  set<TString> data_histkeys_raw=GetHistKeys(0), signal_histkeys_raw=GetHistKeys(1);
  set<TString> prefixes,suffixes;
  for(unsigned int i=0;i<samples.size();i++)
    for(unsigned int j=0;j<samples[i].prefixes.size();j++)
      if(samples[i].prefixes[j]!="")
	prefixes.insert(samples[i].prefixes[j]);
  for(unsigned int i=0;i<systematics.size();i++)
    for(unsigned int j=0;j<systematics[i].suffixes.size();j++)
      if(systematics[i].suffixes[j]!="")
	suffixes.insert(systematics[i].suffixes[j]);
  set<TString> data_histkeys=ParseHistKeys(data_histkeys_raw,prefixes,suffixes,excludes);
  set<TString> signal_histkeys=ParseHistKeys(signal_histkeys_raw,prefixes,suffixes,excludes);
  for(auto it=signal_histkeys.begin();it!=signal_histkeys.end();it++){
    if(data_histkeys.find(*it)!=data_histkeys.end()){
      Plot plot;
      plot.name=*it;plot.name=plot.name(1,plot.name.Length());plot.rebin=0;plot.xmin=0;plot.xmax=0;
      if((*it).Contains("weightedAFB")) plot.option="weightedAFB ";
      else if((*it).Contains("AFB")) plot.option="AFB ";
      plots[*it]=plot;
    }
  }
}
void Plotter::PrintPlots(TRegexp reg){
  for(auto it=plots.begin();it!=plots.end();it++){
    if(it->first.Contains(reg)) cout<<it->first<<endl;
  }
}
*/
TString Plotter::GetHistNameWithPrefix(TString prefix,TString histname){
  TString this_histname=histname(0,histname.Last('/')+1)+prefix+histname(histname.Last('/')+1,histname.Length());
  return this_histname;
}

void Plotter::PrintSampleFrags(){
  for(auto it=samplefrags.begin();it!=samplefrags.end();it++){
    cout<<"@ Key: "<<it->first<<" ";
    it->second.Print();
  }
}

void Plotter::PrintSamples(bool detail){
  for(auto it=samples.begin();it!=samples.end();it++){
    cout<<"@ Key: "<<it->first<<" ";
    it->second.Print(detail);
  }
}
