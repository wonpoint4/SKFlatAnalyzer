#ifndef __PLOTTER_CC__
#define __PLOTTER_CC__
#include"Global.h"
#include"Sample.cc"
#include"Plot.cc"
#include"Utils.h"
#include"Systematic.cc"
class Plotter{
public:
  Plotter();
  ~Plotter();
  vector<Sample> entries;
  map<TString,Systematic> systematics;
  map<TString,Plot> plots;
  map<TString,set<TString>> histkeys_cache;
  static TDirectory *pdir;
  TString plotfile;
  TString plotdir="fig";

  //Setup
  void AddFile(TString key,TString file);
  void ScanFiles(TString path);
  virtual void SetupEntries(TString mode);
  void Reset();
  const Sample& GetEntry(TString entrytitle) const;
  void AddEntry(TString key);

  //Print
  void PrintFiles(bool detail=false,TRegexp regexp=".*");
  void PrintSamples(bool detail=false,TRegexp regexp=".*");
  void PrintEntries(bool detail=false,TRegexp regexp=".*");
  void PrintHistKeys(const Sample& sample,TString regexp=".*");
  void PrintHistKeys(int ientry=0,TString regexp=".*");
  void PrintPlots(TRegexp reg=".*");
  void PrintSystematics();

  //Hist
  virtual TH1* GetHist(TString samplekey,TString plotkey,TString additional_option="");
  virtual TH1* GetHist(int ientry,TString plotkey,TString additional_option="");
  virtual TH1* GetHist(const Sample& sample,Plot plot,TString additional_option="");
  virtual TH1* GetHistFromSample(const Sample& sample,Plot plot,TString additional_option="");
  virtual TH1* GetHistFromFile(TString filename,TString histname);
  virtual void GetHistActionForAdditionalClass(TObject*& obj,Plot plot);
  tuple<TH1*,TH1*> GetHistPair(const Sample& sample,const Plot& plot);
  vector<tuple<TH1*,TH1*>> GetHistPairs(const Plot& plot);

  //Uncertainty
  TH1* GetEnvelope(TH1* central,const vector<TH1*>& variations);
  TH1* GetEnvelope(TH1* central,TH1* variation1,TH1* variation2=NULL,TH1* variation3=NULL,TH1* variation4=NULL,TH1* variation5=NULL,TH1* variation6=NULL,TH1* variation7=NULL,TH1* variation8=NULL,TH1* variation9=NULL);
  TH1* GetHessianError(TH1* central,const vector<TH1*>& variations);
  TH1* GetRMSError(TH1* central,const vector<TH1*>& variations);
  int AddError(TH1* hist,TH1* sys);

  //Canvas
  TCanvas* GetCompare(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  vector<tuple<TH1*,TH1*>> GetRatioHistPairs(vector<tuple<TH1*,TH1*>> histpairs,TString option);
  TCanvas* GetRatio(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetDiff(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetSig(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetCompareAndRatio(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetCompareAndDiff(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetCompareAndSig(vector<tuple<TH1*,TH1*>> hists,Plot plot);
  TCanvas* GetDoubleRatio(Plot& plot);
  virtual TCanvas* DrawPlot(Plot plotkey,TString additional_option="");
  virtual TCanvas* DrawPlot(TString plotkey,TString additional_option="");
  vector<TCanvas*> DrawPlots(TRegexp plotkey,TString additional_option="");
  void SavePlots(TRegexp plotkey);
  void SavePlot(TString plotkey,TString option="",bool delete_canvas=true);
  void SavePlotAll();

  //Canvas tools
  TString GetHistNameWithPrefixAndSuffix(TString histname,TString prefix,TString suffix);
  static TH1* GetAxisParent(TVirtualPad* pad);
  TH1* GetTH1(TH1* hstack,bool deleteorigin=false);
  bool CheckHists(vector<TH1*> hists);
  static vector<TH1*> VectorTH1(vector<tuple<TH1*,TH1*>>& hists);
  tuple<double,double> GetMinMax(const vector<TH1*>& hists);
  void RebinXminXmax(TH1* hist,int rebin,double xmin,double xmax);
  void Normalize(TH1* hists,double val=1.);
  void Normalize(vector<TH1*> hists,double val=1.);
  void WidthWeight(TH1* hist);
  void WidthWeight(vector<TH1*> hists);
  bool IsEntry(const Sample& sample);
  static TLegend* GetLegend(const vector<TH1*>& hists,TString option);
  static TLegend* GetLegend(const vector<tuple<TH1*,TH1*>>& hists,TString option);
  
  //Plot
  void SetupPlots(TString filename);
  set<TString> GetHistKeys(TDirectory* dir,TRegexp regexp=".*");
  set<TString> GetHistKeys(TString filename,TString regexp=".*");
  set<TString> GetHistKeys(const Sample& sample,TString regexp=".*");
  //set<TString> ParseHistKeys(set<TString> histkeys,set<TString> fixes,set<TString> excludes={});
  //set<TString> GetParsedHistKeys(set<TString> excludes={});
  virtual void AddPlot(TString plotkey,TString plotoption="");
  //void AddPlots(TRegexp regexp=".*",TString plotoption="");
  void RemovePlot(TString plotkey);
  void RemovePlots(TRegexp regexp=".*");
  void SetPlotRebinXminXmaxAuto(TString plotkey);
  void SetPlotsRebinXminXmaxAuto();
  void SetPlotOption(TString plotname,TString option,TString suffix="");
  void SetPlotsOption(TRegexp plotname,TString option,TString suffix="");
  void RemovePlotOption(TString plotname,TString option);
  void RemovePlotsOption(TRegexp plotname,TString option);
  void SavePlotFile();

  //MakeSystematic
  Systematic MakeSystematic(TString name,Systematic::Type type,int varibit,vector<TString> includes);
  Systematic MakeSystematic(TString name,Systematic::Type type,int varibit,TString includes);

  //etc
  double GetChi2(TH1* h1,TH1* h2) const;
  int GetAutoColor() const;
  
  //working

  //TCanvas* GetCompareAFBAll(vector<TString> histnames,int sysbit=0,TString option="");
  //TCanvas* GetCompareAFBAll(TRegexp regexp,int sysbit=0,TString option="");
};

TDirectory* Plotter::pdir=NULL;

Plotter::Plotter(){
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(700);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTickLength(0.015,"Y");
}
Plotter::~Plotter(){
  SavePlotFile();
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Setup ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void Plotter::AddFile(TString key,TString file){
  if(samples.find(key)==samples.end()){
    TFile f(file);
    if(f.IsZombie()) PError("[Plotter::AddFile] Wrong file "+file);
    else{
      char buf[64];
      int nbytes=0,objlen=0,keylen=0;
      f.GetRecordHeader(buf,100,64,nbytes,objlen,keylen);
      if(nbytes==1||nbytes>f.GetSize()){
	PWarning("[Plotter::AddFile] potentially corrupted file... recover using hadd");
	system("TEMPDIR=$(mktemp -d);hadd $TEMPDIR/temp.root "+file+";cp $TEMPDIR/temp.root "+file+";rm -r $TEMPDIR");
      }
      samples[key]=Sample(file,Sample::Type::FILE);
    }
  }else PError("[Plotter::AddFile] sample "+key+" already exists");
}
void Plotter::ScanFiles(TString path){
  vector<TString> files=Split(gSystem->GetFromPipe("find "+path+" -type f -name '*.root'"),"\n");
  if(!path.EndsWith("/")) path+="/";
  for(const auto& file:files){
    TString key=Replace(file,path,"");
    key=Replace(key,".root$","");
    AddFile(key,file);
  }
}
void Plotter::SetupEntries(TString mode){
  if(Verbosity) std::cout<<"[Plotter::SetupEntries] mode = '"<<mode<<"'"<<endl;
  entries.clear();
  vector<TString> entry_keys=Split(mode," ");
  for(auto entry_key:entry_keys){
    AddEntry(entry_key);
  }
  return;
}
void Plotter::Reset(){
  entries.clear();
  systematics.clear();
  SavePlotFile();
  plots.clear();
  plotdir="fig";
}
const Sample& Plotter::GetEntry(TString entrytitle) const {
  for(const auto& entry:entries){
    if(entry.title==entrytitle) return entry;
  }
  PError("[Plotter::GetEntry(TString entrytitle)] no entry with title "+entrytitle);
  return entries.at(0);
}
void Plotter::AddEntry(TString key){
  if(samples.find(key)!=samples.end()){
    if(samples[key].IsFile()){
      Sample sample(key,Sample::Type::UNDEF,GetAutoColor());
      sample=sample+samples[key];
      entries.push_back(sample);
    }else entries.push_back(samples[key]);
  }else{
    PError("[Plotter::AddEntry] no sample with key "+key);
  }
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Print ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void Plotter::PrintFiles(bool detail,TRegexp regexp){
  std::cout<<"--------------------------"<<endl;
  for(const auto& [key,sample]:samples){
    if(sample.type!=Sample::Type::FILE) continue;
    if(key.Contains(regexp)){
      std::cout<<"@Key: "<<key<<" "<<endl;
      if(detail) sample.Print(detail,"  ");
      std::cout<<"--------------------------"<<endl;
    }
  }
}
void Plotter::PrintSamples(bool detail,TRegexp regexp){
  std::cout<<"--------------------------"<<endl;
  for(const auto& [key,sample]:samples){
    if(sample.type==Sample::Type::FILE) continue;
    if(key.Contains(regexp)){
      std::cout<<"@Key: "<<key<<" "<<endl;
      sample.Print(detail,"  ");
      std::cout<<"--------------------------"<<endl;
    }
  }
}
void Plotter::PrintEntries(bool detail,TRegexp regexp){
  std::cout<<"--------------------------"<<endl;
  for(const auto& entry:entries){
    if(entry.title.Contains(regexp)){
      entry.Print(detail);
      std::cout<<"--------------------------"<<endl;
    }
  }
}
void Plotter::PrintHistKeys(const Sample& sample,TString regexp=".*"){
  set<TString> histkeys=GetHistKeys(sample,regexp);
  for(const auto& key:histkeys) std::cout<<key<<endl;
}
void Plotter::PrintHistKeys(int ientry,TString regexp=".*"){
  set<TString> histkeys=GetHistKeys(entries[ientry],regexp);
  for(const auto& key:histkeys) std::cout<<key<<endl;
}
void Plotter::PrintPlots(TRegexp reg){
  for(auto it=plots.begin();it!=plots.end();it++){
    if(it->first.Contains(reg)) it->second.Print();
  }
}
void Plotter::PrintSystematics(){
  for(const auto& [key,sys]:systematics){
    sys.Print();
  }
}


/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Hist /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TH1* Plotter::GetHist(TString samplekey,TString plotkey,TString additional_option){
  if(plots.find(plotkey)==plots.end()) AddPlot(plotkey);
  Plot plot=plots[plotkey];
  return GetHist(samples[samplekey],plot,additional_option);
}
TH1* Plotter::GetHist(int ientry,TString plotkey,TString additional_option){
  if(plots.find(plotkey)==plots.end()) AddPlot(plotkey);
  Plot plot=plots[plotkey];
  if(entries.size()<=ientry){
    PError("[Plotter::GetHist] entries.size()<=ientry");
    return NULL;
  }
  return GetHist(entries[ientry],plot,additional_option);
} 
TH1* Plotter::GetHist(const Sample& sample,Plot plot,TString additional_option){
  PInfo("[Plotter::GetHist] "+sample.title+" "+plot.histname+" "+plot.option+" "+additional_option);_depth++;
  TString axisstring=plot.histname(TRegexp("([x-zu,]*)$"));
  TH1* hist=GetHistFromSample(sample,plot,additional_option);
  if(hist){
    hist->SetTitle(plot.name+plot.suffix);
    if(axisstring=="(x)") RebinXminXmax(hist,plot.rebin,plot.xmin,plot.xmax);
    else if(axisstring=="(y)") RebinXminXmax(hist,plot.rebin,plot.ymin,plot.ymax);
    else if(axisstring=="(z)") RebinXminXmax(hist,plot.rebin,plot.zmin,plot.zmax);
    else if(axisstring=="(u)") RebinXminXmax(hist,plot.rebin,plot.umin,plot.umax);
    else RebinXminXmax(hist,plot.rebin,plot.xmin,plot.xmax);
    if(plot.option.Contains("widthweight")) WidthWeight(hist);
    //if(plot.option.Contains("norm")) Normalize(hist);
  }
  _depth--;
  return hist;
}
TH1* Plotter::GetHistFromSample(const Sample& sample,Plot plot,TString additional_option){
  for(auto [reg,newstr]:sample.replace) TPRegexp(reg).Substitute(plot.histname,newstr);
  PInfo("[Plotter::GetHistFromSample] "+sample.title+" "+plot.histname+" "+plot.option+" "+additional_option);_depth++;
  plot.SetOption(additional_option);
  TH1* hist=NULL;
  if(sample.type==Sample::Type::STACK){
    for(int i=sample.subs.size()-1;i>=0;i--){
      TH1* this_hist=GetHistFromSample(sample.subs[i],plot);
      if(this_hist){
	if(!hist){
	  hist=(TH1*)new THStack(sample.title,sample.title);
	}
	((THStack*)hist)->Add(this_hist,"HIST");
      }
    }
  }else if(sample.type==Sample::Type::SUM){
    for(const auto& sub:sample.subs){
      TH1* this_hist=GetHistFromSample(sub,plot);
      if(this_hist){
	if(!hist){
	  hist=(TH1*)this_hist->Clone();
	  hist->SetDirectory(pdir);
	  hist->Reset();
	}
	hist->Add(this_hist);
	delete this_hist;
      }
    }
  }else{
    TString axisstring=plot.histname(TRegexp("([x-zu,]*)$"));
    TString histname=Replace(plot.histname,axisstring,"");
    histname=histname+((1<<sample.type)&plot.varibit?plot.suffix:"");
    for(const auto& file:sample.subs){
      TH1* this_hist=GetHistFromFile(file.title,histname);
      if(this_hist){
	if(strstr(this_hist->ClassName(),"TH1")!=NULL){
	}else if(strstr(this_hist->ClassName(),"TH2")!=NULL){
	  TH2* hist2d=(TH2*)this_hist;
	  int ixmin = plot.xmin ? hist2d->GetXaxis()->FindBin(plot.xmin) : 0 ;
	  int ixmax = plot.xmax ? hist2d->GetXaxis()->FindBin(plot.xmax-0.00001) : -1 ;	
	  int iymin = plot.ymin ? hist2d->GetYaxis()->FindBin(plot.ymin) : 0 ;
	  int iymax = plot.ymax ? hist2d->GetYaxis()->FindBin(plot.ymax-0.00001) : -1 ;
	  if(axisstring=="") this_hist=(TH1*)hist2d->ProjectionY("_py",ixmin,ixmax);
	  else if(axisstring=="(x)") this_hist=(TH1*)hist2d->ProjectionX("_px",iymin,iymax);
	  else if(axisstring=="(y)") this_hist=(TH1*)hist2d->ProjectionY("_py",ixmin,ixmax);
	  else PError("[Plotter::GetHistFromSample] wrong axisstring or classname");
	  delete hist2d;
	}else if(strstr(this_hist->ClassName(),"TH3")!=NULL){
	  TH3* hist3d=(TH3*)this_hist;
	  int ixmin=0,iymin=0,izmin=0;
	  int ixmax=-1,iymax=-1,izmax=-1;
	  if(plot.xmin||plot.xmax){
	    ixmin=hist3d->GetXaxis()->FindBin(plot.xmin);
	    ixmax=hist3d->GetXaxis()->FindBin(plot.xmax-0.00001);
	  }
	  if(plot.ymin||plot.ymax){
	    iymin=hist3d->GetYaxis()->FindBin(plot.ymin);
	    iymax=hist3d->GetYaxis()->FindBin(plot.ymax-0.00001);
	  }
	  if(plot.zmin||plot.zmax){
	    izmin=hist3d->GetZaxis()->FindBin(plot.zmin);
	    izmax=hist3d->GetZaxis()->FindBin(plot.zmax-0.00001);
	  }
	  if(axisstring=="") this_hist=(TH1*)hist3d->ProjectionZ("_pz",ixmin,ixmax,iymin,iymax);
	  else if(axisstring=="(x)") this_hist=(TH1*)hist3d->ProjectionX("_px",iymin,iymax,izmin,izmax);
	  else if(axisstring=="(y)") this_hist=(TH1*)hist3d->ProjectionY("_py",ixmin,ixmax,izmin,izmax);
	  else if(axisstring=="(z)") this_hist=(TH1*)hist3d->ProjectionZ("_pz",ixmin,ixmax,iymin,iymax);
	  else PError("[Plotter::GetHistFromSample] wrong axisstring or classname");
	  delete hist3d;
	}else{
	  TObject* obj=(TObject*)this_hist;
	  GetHistActionForAdditionalClass(obj,plot);
	  this_hist=(TH1*)obj;
	}
	this_hist->SetName(sample.title);
	if(!hist){
	  hist=(TH1*)this_hist->Clone();
	  hist->SetDirectory(pdir);
	  hist->Reset();
	}
	hist->Add(this_hist,sample.weight);
	delete this_hist;
      }
    }
  }
  if(hist){
    sample.ApplyStyle(hist);
  }
  _depth--;
  return hist;
}
TH1* Plotter::GetHistFromFile(TString filename,TString histname){
  TH1* hist=NULL;
  if(histname.MaybeRegexp()){
    for(const auto& key:GetHistKeys(filename,"^"+histname+"$")){
      TH1* this_hist=GetHistFromFile(filename,key);
      if(this_hist){
	if(hist){
	  hist->Add(this_hist);
	  delete this_hist;
	}else hist=this_hist;
      }
    }
  }else{
    TFile *f=TFile::Open(filename);
    if(f->IsZombie()) PError("It is a zombie file");
    if(!f->IsOpen()) PError("Cannot open the file");
    hist=(TH1*)f->Get(histname);
    if(hist){
      hist->SetDirectory(pdir);
      if(!hist->GetSumw2N()) hist->Sumw2();
      PInfo("[Plotter::GetHistFromFile] get "+histname+" from "+filename);
    }
    f->Close();
    delete f;
  }
  if(!hist)
    PDebug("[Plotter::GetHistFromFile] no "+histname+" in "+filename);
  return hist;
}
void Plotter::GetHistActionForAdditionalClass(TObject*& obj,Plot plot){
  PError((TString)"[Plotter::GetHistActionForAdditionalClass(TObject*& obj,Plot plot)] Unsupported class "+obj->ClassName());
}
tuple<TH1*,TH1*> Plotter::GetHistPair(const Sample& sample,const Plot& plot){
  PInfo("[Plotter::GetHistPair] "+sample.title+" "+plot.name+" "+plot.histname);_depth++;

  TH1* central=GetHist(sample,plot);
  TH1 *total=NULL;
  if(plot.sysname!=""){
    THStack* hstack=NULL;
    if(strstr(central->ClassName(),"THStack")!=NULL){
      hstack=(THStack*)central;
      central=GetTH1((TH1*)hstack);
    }
    int sysbit=systematics[plot.sysname].sysbit;
    vector<TH1*> syss;
    for(const auto& sysmap:systematics){
      const Systematic& systematic=sysmap.second;
      if(systematic.type==Systematic::Type::MULTI) continue;
      if(sysbit&systematic.sysbit){
	PDebug("[Plotter::GetHistPair] sysname='"+systematic.name+"' systype='"+systematic.GetTypeString()+"'");
	vector<TH1*> variations;
	for(const auto& suffix:systematic.suffixes){
	  TH1* this_hist=GetHist(sample,plot,Form("suffix:%s varibit:%d sysname:",suffix.Data(),systematic.varibit));
	  variations.push_back(GetTH1(this_hist,true));
	}
	if(systematic.type==Systematic::Type::ENVELOPE){
	  syss.push_back(GetEnvelope(central,variations));
	}else if(systematic.type==Systematic::Type::GAUSSIAN){
	  syss.push_back(GetRMSError(central,variations));
	}else if(systematic.type==Systematic::Type::HESSIAN){
	  syss.push_back(GetHessianError(central,variations));
	}else{
	  PError((TString)"[Plotter::GetHistPair] Wrong Systematic::Type "+Form("%d",systematic.type));
	}
	PDebug("[Plotter::GetHistPair] '"+systematic.name+"': "+variations.size()+" variations");
	for(unsigned int j=0;j<variations.size();j++){
	  delete variations.at(j);
	}
      }
    }
    if(syss.size()>0){
      total=(TH1*)central->Clone();
      total->SetDirectory(pdir);
      for(int i=0;i<(int)syss.size();i++){
	AddError(total,syss.at(i));
	delete syss.at(i);
      }
    }
    if(hstack){
      delete central;
      central=(TH1*)hstack;
    }
    sample.ApplyStyle(total,true);
  }
  if(plot.option.Contains("norm")){ //FIXME: Why here?
    Normalize(total);
    Normalize(central);
  }
  _depth--;
  return make_tuple(central,total);
}
vector<tuple<TH1*,TH1*>> Plotter::GetHistPairs(const Plot& plot){
  vector<tuple<TH1*,TH1*>> histpairs;
  for(const auto& entry:entries) histpairs.push_back(GetHistPair(entry,plot));
  return histpairs;
}  


/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Uncertainty //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TH1* Plotter::GetEnvelope(TH1* central,const vector<TH1*>& variations){
  PAll("[Plotter::GetEnvelope]");
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetDirectory(pdir);
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
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetEnvelope(TH1* central,TH1* variation1,TH1* variation2=NULL,TH1* variation3=NULL,TH1* variation4=NULL,TH1* variation5=NULL,TH1* variation6=NULL,TH1* variation7=NULL,TH1* variation8=NULL,TH1* variation9=NULL)]"<<endl;
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
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetHessianError(TH1* central,const vector<TH1*>& variations)]"<<endl;
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetDirectory(pdir);
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
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetRMSError(TH1* central,const vector<TH1*>& variations)]"<<endl;
  if(strstr(central->ClassName(),"THStack")) central=GetTH1(central);
  TH1* syshist=(TH1*)central->Clone("sys");
  syshist->SetDirectory(pdir);
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
  PAll("[Plotter::AddError]");
  for(int i=1;i<hist->GetNbinsX()+1;i++){
    if(fabs(hist->GetBinContent(i)-sys->GetBinContent(i))*1000000>fabs(hist->GetBinContent(i))){
      PError("[Plotter::AddError] systematic hist is wrong");
      std::cout.precision(20);
      std::cout<<i<<" "<<hist->GetBinContent(i)<<" "<<sys->GetBinContent(i)<<" "<<fabs(hist->GetBinContent(i)-sys->GetBinContent(i))<<endl;
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


/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Canvas ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TCanvas* Plotter::GetCompare(vector<tuple<TH1*,TH1*>> histpairs,Plot plot){
  plot.Print();
  PInfo("[Plotter::GetCompare] "+plot.option);_depth++;
  vector<TH1*> hists=VectorTH1(histpairs);
  TH1* axisowner=get<1>(histpairs[0]);
  if(!axisowner) axisowner=get<0>(histpairs[0]);
  
  TCanvas* c=new TCanvas;
  if(strstr(axisowner->ClassName(),"THStack")==NULL){
    axisowner->SetStats(0);
    axisowner->Draw(axisowner->GetOption());
  }else{
    axisowner->Draw();
    axisowner=((THStack*)axisowner)->GetHistogram();
    axisowner->SetName("hframe");
    axisowner->Draw();
  }
  if(c->GetPrimitive("title")) ((TPaveText*)c->GetPrimitive("title"))->SetTextSize(0.075);

  for(const auto& hist:hists)
    if(strstr(hist->ClassName(),"THStack")){
      hist->Draw("same");
      TH1* hist_sum=GetTH1(hist);
      hist_sum->Draw(TString("same ")+hist->GetOption());
    }

  for(const auto& hist:hists)
    if(strstr(hist->ClassName(),"THStack")==NULL&&hist->GetFillColor()>0)
      hist->Draw(TString("same ")+hist->GetOption());

  for(const auto& hist:hists)
    if(strstr(hist->ClassName(),"THStack")==NULL&&hist->GetFillColor()==0)
      hist->Draw(TString("same ")+hist->GetOption());

  TLegend* legend=GetLegend(histpairs,plot.option);
  if(!plot.option.Contains("noleg")) legend->Draw();

  if(plot.option.Contains("logx")){c->SetLogx();axisowner->GetXaxis()->SetMoreLogLabels();}
  if(plot.option.Contains("logy")){
    tuple<double,double> minmax=GetMinMax(hists);
    double minimum=get<0>(minmax),maximum=get<1>(minmax);
    if(minimum<0) minimum=maximum/1000;
    axisowner->GetYaxis()->SetRangeUser(minimum/100,maximum*20);
    c->SetLogy();
  }else{
    tuple<double,double> minmax=GetMinMax(hists);
    double minimum=get<0>(minmax),maximum=get<1>(minmax);
    double range=fabs(maximum-minimum);
    axisowner->GetYaxis()->SetRangeUser(minimum/range<-0.01?minimum-0.1*range:0,maximum+0.1*range);
  }
  axisowner->GetXaxis()->SetTitle(plot.xtitle);
  axisowner->GetYaxis()->SetTitle(plot.ytitle);
  gPad->RedrawAxis();
  _depth--;
  return c;
}
vector<tuple<TH1*,TH1*>> Plotter::GetRatioHistPairs(vector<tuple<TH1*,TH1*>> histpairs,TString option){
  vector<tuple<TH1*,TH1*>> histpairs_new;
  TH1* base=NULL;
  TString baseopt=option("base:[0-9]*");
  if(baseopt!="") base=GetTH1(get<0>(histpairs[((TString)baseopt(5,100)).Atoi()]));
  else if(histpairs.size()>2) base=GetTH1(get<0>(histpairs[0]));
  else base=GetTH1(get<0>(histpairs[1]));
    
  if(base){
    for(int i=0;i<base->GetNbinsX()+2;i++){
      base->SetBinError(i,0);
    }
    for(const auto& [hist0,hist1]:histpairs){
      TH1 *hist0_new=NULL,*hist1_new=NULL;
      if(hist0){
	hist0_new=GetTH1(hist0);
	hist0_new->Divide(base);
      }
      if(hist1){
	hist1_new=GetTH1(hist1);
	hist1_new->Divide(base);
      }
      histpairs_new.push_back(make_tuple(hist0_new,hist1_new));
    }  
    delete base;
  }
  return histpairs_new;
}  
TCanvas* Plotter::GetRatio(vector<tuple<TH1*,TH1*>> histpairs,Plot plot){
  PInfo("[Plotter::GetRatio] "+plot.option);_depth++;
  vector<tuple<TH1*,TH1*>> histpairs_new=GetRatioHistPairs(histpairs,plot.option);
  TCanvas* c=GetCompare(histpairs_new,plot-"norm");
  TH1* axisowner=GetAxisParent(c);
  if(axisowner){
    axisowner->GetYaxis()->SetTitle("Ratio");
    if(plot.option.Contains("widewidey")){
      axisowner->GetYaxis()->SetRangeUser(0.01,1.99);
      axisowner->GetYaxis()->SetNdivisions(506);
    }else if(plot.option.Contains("widey")){
      axisowner->GetYaxis()->SetRangeUser(0.501,1.499);
      axisowner->GetYaxis()->SetNdivisions(506);
    }else{
      axisowner->GetYaxis()->SetRangeUser(0.801,1.199);
      axisowner->GetYaxis()->SetNdivisions(504);
    }
  }
  _depth--;
  return c;
}
TCanvas* Plotter::GetDiff(vector<tuple<TH1*,TH1*>> tu_hists,Plot plot){
  PInfo("[Plotter::GetDiff]"+plot.option);
  vector<tuple<TH1*,TH1*>> tu_hists_new;
  TH1* base=NULL;
  TString baseopt=plot.option("base:[0-9]*");
  if(baseopt!="") base=GetTH1(get<0>(tu_hists[((TString)baseopt(5,100)).Atoi()]));
  else if(tu_hists.size()>2) base=GetTH1(get<0>(tu_hists[0]));
  else base=GetTH1(get<0>(tu_hists[1]));
  
  for(int i=0;i<base->GetNbinsX()+2;i++){
    base->SetBinError(i,0);
  }
  for(const auto& histtuple:tu_hists){
    TH1 *hist0=NULL,*hist1=NULL;
    if(get<0>(histtuple)){
      hist0=GetTH1(get<0>(histtuple));
      hist0->Add(base,-1);
    }
    if(get<1>(histtuple)){
      hist1=(TH1*)get<1>(histtuple)->Clone();
      hist1->SetDirectory(pdir);
      hist1->Add(base,-1);
    }
    tu_hists_new.push_back(make_tuple(hist0,hist1));
  }  
  delete base;
  TCanvas* c=GetCompare(tu_hists_new,plot);
  TH1* axisowner=GetAxisParent(c);
  
  axisowner->GetYaxis()->SetTitle("Diff");
  axisowner->GetYaxis()->SetLabelSize(0.06);

  tuple<double,double> minmax=GetMinMax(VectorTH1(tu_hists_new));
  double minimum=get<0>(minmax),maximum=get<1>(minmax);
  double range=fabs(maximum-minimum);
  axisowner->GetYaxis()->SetRangeUser(minimum/range<-0.01?minimum-0.1*range:0,maximum+0.1*range);
  _depth--;
  return c;
}
TCanvas* Plotter::GetSig(vector<tuple<TH1*,TH1*>> tu_hists,Plot plot){
  PInfo("[Plotter::GetSig] "+plot.option);_depth++;
  vector<TH1*> hists;
  for(auto tu_hist:tu_hists){
    if(get<1>(tu_hist)) hists.push_back(get<1>(tu_hist));
    else hists.push_back(get<0>(tu_hist));
  }
  TH1* base=NULL;
  TString baseopt=plot.option("base:[0-9]*");
  if(baseopt!="") base=GetTH1(hists[((TString)baseopt(5,100)).Atoi()]);
  else if(tu_hists.size()>2) base=GetTH1(hists[0]);
  else base=GetTH1(hists[1]);
  
  vector<tuple<TH1*,TH1*>> tu_hists_new;
  for(auto hist:hists){
    if(hist){
      hist=GetTH1(hist,true);
      for(int i=0;i<hist->GetNbinsX()+2;i++){
	double error=sqrt(pow(hist->GetBinError(i),2)+pow(base->GetBinError(i),2));
	double content=(hist->GetBinContent(i)-base->GetBinContent(i))/error;
	hist->SetBinContent(i,content);
	hist->SetBinError(i,0);
      }
    }
    tu_hists_new.push_back(make_tuple(hist,(TH1*)NULL));
  }
  TCanvas* c=GetCompare(tu_hists_new,plot);
  TH1* axisowner=GetAxisParent(c);
  
  axisowner->GetYaxis()->SetTitle("Difference (#sigma)");
  axisowner->GetYaxis()->SetLabelSize(0.06);

  if(plot.option.Contains("widey")) axisowner->GetYaxis()->SetRangeUser(-4.9,4.9);
  else axisowner->GetYaxis()->SetRangeUser(-2.9,2.9);
  //tuple<double,double> minmax=GetMinMax(VectorTH1(tu_hists_new));
  //double minimum=get<0>(minmax),maximum=get<1>(minmax);
  //double range=fabs(maximum-minimum);
  //axisowner->GetYaxis()->SetRangeUser(minimum/range<-0.01?minimum-0.1*range:0,maximum+0.1*range);
  _depth--;
  return c;
}
TCanvas* Plotter::GetCompareAndRatio(vector<tuple<TH1*,TH1*>> hists,Plot plot){
  PInfo("[Plotter::GetCompareAndRatio] "+plot.option);_depth++;
    
  TCanvas* c=new TCanvas;
  c->Divide(1,2);
  TCanvas* c1temp=GetCompare(hists,plot-"2:");
  c->cd(1);
  c1temp->DrawClonePad();
  delete c1temp;
  if(c->GetPad(1)->GetPrimitive("title")) ((TPaveText*)c->GetPad(1)->GetPrimitive("title"))->SetTextSize(0.075);
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);
  TH1* axisparent=GetAxisParent(gPad);
  axisparent->GetYaxis()->SetLabelSize(0.06);
  axisparent->GetYaxis()->SetTitleSize(0.06);
  axisparent->GetYaxis()->SetTitleOffset(1.2);
  axisparent->GetXaxis()->SetLabelSize(0);
  axisparent->GetXaxis()->SetTitle("");
  gPad->Update();
  gPad->Modified();

  TCanvas* c2temp=GetRatio(hists,plot-"1:"+"noleg");
  c->cd(2);
  if(c2temp){
    c2temp->DrawClonePad();
    delete c2temp;
  }
  gPad->SetPad(0,0,1,0.365);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridx();gPad->SetGridy();
  gPad->SetFillStyle(0);
  axisparent=GetAxisParent(gPad);
  if(axisparent){
    axisparent->SetTitle("");
    axisparent->SetStats(0);
    axisparent->GetYaxis()->SetLabelSize(0.1);
    axisparent->GetYaxis()->SetTitle("Ratio");
    axisparent->GetYaxis()->SetTitleSize(0.12);
    axisparent->GetYaxis()->SetTitleOffset(0.6);
    axisparent->GetXaxis()->SetTitle(plot.xtitle);
    axisparent->GetXaxis()->SetTitleSize(0.12);
    axisparent->GetXaxis()->SetLabelSize(0.12);
  }
  gPad->Update();
  gPad->Modified();
  _depth--;
  return c;
}
TCanvas* Plotter::GetCompareAndSig(vector<tuple<TH1*,TH1*>> hists,Plot plot){
  PInfo("[Plotter::GetCompareAndSig] "+plot.option);_depth++;
    
  TCanvas* c=new TCanvas;
  c->Divide(1,2);
  TCanvas* c1temp=GetCompare(hists,plot-"2:");
  c->cd(1);
  c1temp->DrawClonePad();
  delete c1temp;
  if(c->GetPad(1)->GetPrimitive("title")) ((TPaveText*)c->GetPad(1)->GetPrimitive("title"))->SetTextSize(0.075);
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);
  TH1* axisparent=GetAxisParent(gPad);
  axisparent->GetYaxis()->SetLabelSize(0.06);
  axisparent->GetYaxis()->SetTitleSize(0.06);
  axisparent->GetYaxis()->SetTitleOffset(1.2);
  axisparent->GetXaxis()->SetLabelSize(0);
  axisparent->GetXaxis()->SetTitle("");
  gPad->Update();
  gPad->Modified();

  TCanvas* c2temp=GetSig(hists,plot-"1:"+"noleg");
  c->cd(2);
  if(c2temp){
    c2temp->DrawClonePad();
    delete c2temp;
  }
  gPad->SetPad(0,0,1,0.365);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridy();
  gPad->SetFillStyle(0);
  axisparent=GetAxisParent(gPad);
  if(axisparent){
    axisparent->SetTitle("");
    axisparent->SetStats(0);
    axisparent->GetYaxis()->SetLabelSize(0.1);
    axisparent->GetYaxis()->SetTitleSize(0.12);
    axisparent->GetYaxis()->SetTitleOffset(0.6);
    axisparent->GetXaxis()->SetTitle(axisparent->GetTitle());
    axisparent->GetXaxis()->SetTitleSize(0.12);
    axisparent->GetXaxis()->SetLabelSize(0.12);
  }
  gPad->Update();
  gPad->Modified();
  _depth--;
  return c;
}
TCanvas* Plotter::GetCompareAndDiff(vector<tuple<TH1*,TH1*>> hists,Plot plot){
  PInfo("[Plotter::GetCompareAndDiff] "+plot.option);_depth++;
    
  TCanvas* c=new TCanvas;
  c->Divide(1,2);
  TCanvas* c1temp=GetCompare(hists,plot-"2:");
  c->cd(1);
  c1temp->DrawClonePad();
  delete c1temp;
  if(c->GetPad(1)->GetPrimitive("title")) ((TPaveText*)c->GetPad(1)->GetPrimitive("title"))->SetTextSize(0.075);
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);
  TH1* axisparent=GetAxisParent(gPad);
  axisparent->GetYaxis()->SetLabelSize(0.06);
  axisparent->GetYaxis()->SetTitleSize(0.06);
  axisparent->GetYaxis()->SetTitleOffset(1.2);
  axisparent->GetXaxis()->SetLabelSize(0);
  axisparent->GetXaxis()->SetTitle("");
  gPad->Update();
  gPad->Modified();

  TCanvas* c2temp=GetDiff(hists,plot-"1:"+"noleg");
  c->cd(2);
  if(c2temp){
    c2temp->DrawClonePad();
    delete c2temp;
  }
  gPad->SetPad(0,0,1,0.365);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.3);
  gPad->SetGridx();gPad->SetGridy();
  gPad->SetFillStyle(0);

  axisparent=GetAxisParent(gPad);
  if(axisparent){
    axisparent->SetTitle("");
    axisparent->SetStats(0);
    axisparent->GetYaxis()->SetLabelSize(0.1);
    axisparent->GetYaxis()->SetTitleSize(0.12);
    axisparent->GetYaxis()->SetTitleOffset(0.6);
    axisparent->GetYaxis()->SetNdivisions(305);
    axisparent->GetXaxis()->SetTitle(axisparent->GetTitle());
    axisparent->GetXaxis()->SetTitleSize(0.12);
    axisparent->GetXaxis()->SetLabelSize(0.12);
  }
  gPad->Update();
  gPad->Modified();
  _depth--;
  return c;
}
TCanvas* Plotter::GetDoubleRatio(Plot& plot){
  vector<vector<tuple<TH1*,TH1*>>> vec_histpairs;
  for(auto& sub:plot.subplots){
    vec_histpairs.push_back(GetRatioHistPairs(GetHistPairs(sub),plot.option));
  }
  int baseindex=0;
  if(plot.subplots.size()==2){
    baseindex=1;
    TString baseopt=plot.option("base:[0-9]*");
    if(baseopt!="") baseindex=((TString)baseopt(5,100)).Atoi();
    for(int i=0;i<vec_histpairs[baseindex].size();i++){
      TH1* hist=get<0>(vec_histpairs[baseindex][i]);
      bool flag=false;
      for(int j=1;j<hist->GetNbinsX()+1;j++) if(hist->GetBinContent(j)!=1.) flag=true;
      if(!flag){
	for(int j=0;j<vec_histpairs.size();j++){
	  vec_histpairs[j].erase(vec_histpairs[j].begin()+i);
	}
      }
    }
    vector<tuple<TH1*,TH1*>> histpairs_new;
    for(int i=0;i<vec_histpairs[baseindex].size();i++){
      TH1* base=(TH1*)get<0>(vec_histpairs[baseindex][i])->Clone();
      for(int j=0;j<base->GetNbinsX()+2;j++) base->SetBinError(j,0);
      for(int j=0;j<vec_histpairs.size();j++){
	if(get<0>(vec_histpairs[j][i])) get<0>(vec_histpairs[j][i])->Divide(base);
	if(get<1>(vec_histpairs[j][i])) get<1>(vec_histpairs[j][i])->Divide(base);
	histpairs_new.push_back(make_tuple(get<0>(vec_histpairs[j][i]),get<1>(vec_histpairs[j][i])));
      }
      delete base;
    }
    TCanvas* c=GetCompare(histpairs_new,plot.option.ReplaceAll("norm",""));
    TH1* axisowner=GetAxisParent(c);
    axisowner->SetTitle(plot.name);
    c->SetGridy();
    axisowner->GetYaxis()->SetTitle("DoubleRatio");
    if(plot.option.Contains("widewidey")){
      axisowner->GetYaxis()->SetRangeUser(0.01,1.99);
      axisowner->GetYaxis()->SetNdivisions(506);
    }else if(plot.option.Contains("widey")){
      axisowner->GetYaxis()->SetRangeUser(0.501,1.499);
      axisowner->GetYaxis()->SetNdivisions(506);
    }else{
      axisowner->GetYaxis()->SetRangeUser(0.801,1.199);
      axisowner->GetYaxis()->SetNdivisions(504);
    }
    return c;
  }else{
    if(DEBUG>0) std::cout<<"###ERROR### not yet"<<endl;
    exit(1);
  }
}

TCanvas* Plotter::DrawPlot(Plot plot,TString additional_option){
  PInfo("[Plotter::DrawPlot] "+plot.name+" "+plot.histname+" "+plot.option+" "+additional_option);_depth++;
  plot.SetOption(additional_option);
  TCanvas* c=NULL;
  if(plot.type==Plot::Type::DoubleRatio){
    c=GetDoubleRatio(plot);
  }else{
    vector<tuple<TH1*,TH1*>> histpairs=GetHistPairs(plot);
    vector<TH1*> hists=VectorTH1(histpairs);
    if(hists.size()<entries.size()){
      PError("[Plotter::DrawPlot] #hists < #entries ");
      _depth--;
      return NULL;
    }
    if(plot.type==Plot::Type::Compare) c=GetCompare(histpairs,plot);
    else if(plot.type==Plot::Type::Ratio) c=GetRatio(histpairs,plot);
    else if(plot.type==Plot::Type::Diff) c=GetDiff(histpairs,plot);
    else if(plot.type==Plot::Type::Sig) c=GetSig(histpairs,plot);
    else if(plot.type==Plot::Type::CompareAndRatio) c=GetCompareAndRatio(histpairs,plot);
    else if(plot.type==Plot::Type::CompareAndDiff) c=GetCompareAndDiff(histpairs,plot);
    else if(plot.type==Plot::Type::CompareAndSig) c=GetCompareAndSig(histpairs,plot);
  }    
  c->Update();
  c->Modified();
  _depth--;
  return c;
}
TCanvas* Plotter::DrawPlot(TString plotkey,TString additional_option){
  if(plots.find(plotkey)==plots.end()) AddPlot(plotkey);
  Plot plot=plots[plotkey];
  return DrawPlot(plot,additional_option);
}
vector<TCanvas*> Plotter::DrawPlots(TRegexp plotkey,TString additional_option){
  vector<TCanvas*> canvases;
  for(const auto& [key,plot]:plots)
    if(key.Contains(plotkey)) canvases.push_back(DrawPlot(key,additional_option));
  return canvases;
}
void Plotter::SavePlot(TString plotkey,TString option,bool delete_canvas){
  PDebug("[Plotter::SavePlot(TString plotkey)]");
  //pdir=new TDirectory;
  TString format="png";
  if(option.Contains("pdf")) format="pdf";
  TCanvas* c=DrawPlot(plotkey,option);
  gSystem->Exec("mkdir -p "+plotdir+"/"+Dirname(plotkey));
  c->SaveAs(plotdir+"/"+plotkey+"."+format);
  if(delete_canvas){
    delete c;
    //pdir->Delete();
  }
}
void Plotter::SavePlots(TRegexp plotkey){
  for(const auto& [key,plot]:plots)
    if(key.Contains(plotkey)) SavePlot(key);
}  
void Plotter::SavePlotAll(){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::SavePlotAll(TString plotkey)]"<<endl;
  for(const auto& [key,plot]:plots) SavePlot(key);
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Canvas Tools /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
TString Plotter::GetHistNameWithPrefixAndSuffix(TString histname,TString prefix,TString suffix){
  PError("[Plotter::GetHistNameWithPrefixAndSuffix] "+histname+" "+prefix+" "+suffix);
  TString this_histname=histname(0,histname.Last('/')+1)+prefix+histname(histname.Last('/')+1,histname.Length())+suffix;
  return this_histname;
}
TH1* Plotter::GetAxisParent(TVirtualPad* pad){
  PAll("[Plotter::GetAxisParent(TVirtualPad* pad)]");
  if(pad){
    if(pad->FindObject("hframe")) return (TH1*)pad->FindObject("hframe");
    TList* list=pad->GetListOfPrimitives();
    for(int i=0;i<list->GetSize();i++){
      if(strstr(list->At(i)->ClassName(),"TH")!=NULL) return (TH1*)list->At(i);
    }
  }
  return NULL;
}
TH1* Plotter::GetTH1(TH1* hstack,bool deleteorigin){
  PAll("[Plotter::GetTH1(TH1* hstack)]");
  TH1* hist=NULL;
  if(hstack){
    if(strstr(hstack->ClassName(),"THStack")!=NULL){
      hist=(TH1*)((THStack*)hstack)->GetHists()->Last()->Clone();
      hist->SetDirectory(pdir);
      hist->Reset();
      for(const auto& entry:entries){
	if(entry.title==hstack->GetName()){
	  entry.ApplyStyle(hist);
	  hist->SetTitle(hstack->GetTitle());
	  break;
	}
      }
      for(const auto& this_hist:*((THStack*)hstack)->GetHists()){
	hist->Add((TH1*)this_hist);
      }
    }else{
      hist=(TH1*)hstack->Clone();
      hist->SetDirectory(pdir);
    }
    if(deleteorigin) delete hstack;
  }
  return hist;
}
bool Plotter::CheckHists(vector<TH1*> hists){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::CheckHists(vector<TH1*> hists)]"<<endl;
  bool flag_data=false,flag_signal=false;
  for(unsigned int i=0;i<entries.size();i++){
    if(!hists.at(i)) continue;
    if(entries[i].type==Sample::Type::DATA) flag_data=true;
    if(entries[i].type==Sample::Type::SIGNAL) flag_signal=true;
  }
  if(flag_data&&flag_signal) return true;
  else return false;
}
vector<TH1*> Plotter::VectorTH1(vector<tuple<TH1*,TH1*>>& hists){
  PAll("[Plotter::VectorTH1(const vector<tuple<TH1*,TH1*>>& hists)]");
  vector<TH1*> hists_out;
  for(const auto& tu:hists){
    TH1* hist0=get<0>(tu); TH1* hist1=get<1>(tu); 
    if(hist1) hists_out.push_back(hist1);
    if(hist0) hists_out.push_back(hist0);
  }
  return hists_out;
}
tuple<double,double> Plotter::GetMinMax(const vector<TH1*>& hists){
  PAll("[Plotter::GetMinMax(vector<TH1*> hists)]");
  double maximum=-999999999;
  double minimum=999999999;
  for(auto hist:hists){
    if(hist){
      if(strstr(hist->ClassName(),"THStack")) hist=GetTH1(hist);
      double this_maximum=hist->GetBinContent(hist->GetMaximumBin());
      double this_minimum=hist->GetBinContent(hist->GetMinimumBin());
      if(maximum<this_maximum) maximum=this_maximum;
      if(minimum>this_minimum) minimum=this_minimum;
    }
  }
  return make_tuple(minimum,maximum);
}
void Plotter::RebinXminXmax(TH1* hist,int rebin,double xmin,double xmax){
  PAll("[Plotter::RebinXminXmax(TH1* hist,int rebin,double xmin,double xmax)]");
  if(hist){
    if(strstr(hist->ClassName(),"THStack")){
      for(const auto& obj:*((THStack*)hist)->GetHists())
	RebinXminXmax((TH1*)obj,rebin,xmin,xmax);
    }else{
      if(rebin){
	double this_xmin=hist->GetBinLowEdge(hist->GetXaxis()->GetFirst());
	double this_xmax=hist->GetBinLowEdge(hist->GetXaxis()->GetLast()+1);
	hist->Rebin(rebin);
	hist->GetXaxis()->SetRangeUser(this_xmin,this_xmax);
      }
      if(xmin==0&&xmax==0){
	if(hist->GetXaxis()->GetFirst()==1&&hist->GetXaxis()->GetLast()==hist->GetNbinsX()){
	  int start=0,end=0;
	  for(int i=0;i<hist->GetNbinsX()+2;i++){
	    if(hist->GetBinContent(i)&&hist->GetBinError(i)){
	      start=i;
	      break;
	    }
	  }
	  for(int i=hist->GetNbinsX()+1;i>=0;i--){
	    if(hist->GetBinContent(i)&&hist->GetBinError(i)){
	      end=i;
	      break;
	    }
	  }
	  if(start<end) hist->GetXaxis()->SetRange(start,end);
	}
      }else hist->GetXaxis()->SetRangeUser(xmin,xmax);
    }
  }
}

void Plotter::Normalize(TH1* hist,double val){
  PAll("[Plotter::Normalize(TH1* hists,double val=1.)]");
  if(hist){
    if(strstr(hist->ClassName(),"THStack")){
      TH1* hsim=GetTH1(hist);
      double scale=val/hsim->Integral();
      for(auto obj:*((THStack*)hist)->GetHists()){
	TH1* hsim_sub=(TH1*)obj;
	hsim_sub->Scale(scale);
      }
    }else{
      hist->Scale(val/hist->Integral());
    }
  }
}
void Plotter::Normalize(vector<TH1*> hists,double val){
  for(auto& hist:hists) Normalize(hist,val);
}
void Plotter::WidthWeight(TH1* hist){
  PAll("[Plotter::WidthWeight(TH1* hist)]");
  if(hist){
    if(strstr(hist->ClassName(),"THStack")){
      TH1* hsim=GetTH1(hist);
      for(auto obj:*((THStack*)hist)->GetHists()){
	TH1* hsim_sub=(TH1*)obj;
	WidthWeight(hsim_sub);
      }
    }else{
      for(int i=1;i<hist->GetNbinsX()+1;i++){
	double width=hist->GetBinWidth(i);
	hist->SetBinContent(i,hist->GetBinContent(i)/width);
	hist->SetBinError(i,hist->GetBinError(i)/width);
      }
    }
  }
}
void Plotter::WidthWeight(vector<TH1*> hists){
  for(auto& hist:hists) WidthWeight(hist);
}
bool Plotter::IsEntry(const Sample& sample){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::IsEntry(const Sample& sample)]"<<endl;  
  for(const auto& entry:entries){
    if(&entry==&sample) return true;
  }
  return false;
}

TLegend* Plotter::GetLegend(const vector<TH1*>& hists,TString option){
  PAll("[Plotter::GetLegend(const vector<TH1*>& hists,TString option)]");
  vector<TH1*> hists_th1;
  for(const auto& hist:hists){
    if(hist){
      if(strstr(hist->ClassName(),"THStack")){
	TIter next(((THStack*)hist)->GetHists(),kIterBackward);
	while(TObject *obj = next()){
	  hists_th1.push_back((TH1*)obj);
	}
      }else hists_th1.push_back(hist);
    }
  }
  int maxlen=0;
  for(const auto& hist:hists_th1){
    TString name=hist->GetName();
    name.ReplaceAll("#rightarrow",">");
    name.ReplaceAll("#gamma","a");
    name.ReplaceAll("#tau","t");
    if(name.Length()>maxlen) maxlen=name.Length();
  }

  if(maxlen<7) maxlen=7;
  if(maxlen>10) maxlen=10;
  int entry_size=hists_th1.size();
  if(entry_size>3) entry_size=3;
  
  vector<double> coordinates;
  double char_unit=0.02,entry_unit=0.13;
  if(option.Contains("topleftleg")||option.Contains("TLleg")) coordinates={0.16,0.89,0.16+maxlen*char_unit,0.89-entry_size*entry_unit};
  else if(option.Contains("topmiddleleg")||option.Contains("TMleg")) coordinates={0.5-0.5*maxlen*char_unit,0.89,0.5+0.5*maxlen*char_unit,0.89-entry_size*entry_unit};
  else if(option.Contains("bottomleftleg")||option.Contains("BLleg")) coordinates={0.16,0.05,0.16+maxlen*char_unit,0.05+entry_size*entry_unit};
  else if(option.Contains("bottommiddleleg")||option.Contains("BMleg")) coordinates={0.5-0.5*maxlen*char_unit,0.05,0.5+0.5*maxlen*char_unit,0.05+entry_size*entry_unit};
  else if(option.Contains("bottomrightleg")||option.Contains("BRleg")) coordinates={0.89,0.05,0.89-maxlen*char_unit,0.05+entry_size*entry_unit};
  else coordinates={0.89,0.89,0.89-maxlen*char_unit,0.89-entry_size*entry_unit};
  TLegend* legend=new TLegend(coordinates[0],coordinates[1],coordinates[2],coordinates[3]);

  for(const auto& hist:hists_th1) legend->AddEntry(hist,hist->GetName());
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

/////////////////////////////////////////////////////////////////////////////
////////////////////////////// Plot /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void Plotter::SetupPlots(TString filename){
  TPRegexp("[ *]").Substitute(filename,"_","g");
  plotdir=Dirname(filename);
  plotfile=Basename(filename);
  gSystem->Exec("mkdir -p "+plotdir);
  ifstream f(filename.Data());
  if(f.fail()){
    //UpdatePlots();
  }else{
    TString str;
    while(str.ReadLine(f)){
      Plot plot(str);
      plots[plot.name]=plot;
    }
  }
}
set<TString> Plotter::GetHistKeys(TDirectory* dir,TRegexp regexp=".*"){
  PAll("[Plotter::GetHistKeys(TDirectory* dir,TRegexp regexp=\".*\")]");
  set<TString> histkeys;
  for(const auto& obj:*(dir->GetListOfKeys())){
    TKey* key=(TKey*)obj;
    if(strcmp(key->GetClassName(),"TDirectoryFile")==0){
      set<TString> this_histkeys=GetHistKeys((TDirectory*)key->ReadObj(),regexp);
      histkeys.insert(this_histkeys.begin(),this_histkeys.end());
    }else{
      TString path=key->GetMotherDir()->GetPath();
      path=path(path.Index(":")+2,path.Length());
      if(path!="") path+="/";
      path+=key->GetName();
      if(path.Contains(regexp)) histkeys.insert(path);
    }
  }
  return histkeys;
}
set<TString> Plotter::GetHistKeys(TString filename,TString regexp=".*"){
  PAll("[Plotter::GetHistKeys(TString filename=\""+filename+",TString regexp=\""+regexp+"\")]");
  if(histkeys_cache.find(filename)==histkeys_cache.end()){
    TFile *f=new TFile(filename);
    histkeys_cache[filename]=GetHistKeys((TDirectory*)f,".*");
    f->Close();delete f;
  }
  if(regexp==".*") return histkeys_cache[filename];
  else{
    set<TString> histkeys;
    for(const TString& key:histkeys_cache[filename])
      if(key.Contains(TRegexp(regexp))) histkeys.insert(key);
    return histkeys;
  }
}
set<TString> Plotter::GetHistKeys(const Sample& sample,TString regexp=".*"){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetHistKeys(const Sample& sample,TString regexp=\".*\")]"<<endl;
  set<TString> histkeys;
  if(sample.type==Sample::Type::FILE){
    set<TString> this_histkeys=GetHistKeys(sample.title,regexp);
    histkeys.insert(this_histkeys.begin(),this_histkeys.end());
  }else if(sample.subs.size()){
    for(const auto& sub:sample.subs){
      set<TString> this_histkeys=GetHistKeys(sub,regexp);
      histkeys.insert(this_histkeys.begin(),this_histkeys.end());
    }
  }
  return histkeys;
}
/*
set<TString> Plotter::ParseHistKeys(set<TString> histkeys,set<TString> fixes,set<TString> excludes={}){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::ParseHistKeys(set<TString> histkeys_raw,set<TString> fixes,set<TString> excludes]"<<endl;
  set<TString> histkeys_out;
  for(auto key:histkeys){
    bool next_flag=false;

    for(const auto& exclude:excludes){
      if(key.Contains(TRegexp(exclude))){
	next_flag=true;
	break;
      }
    }
    if(next_flag) continue;

    for(const auto& fix:fixes){
      TString dirname=Dirname(key);
      TString basename=Basename(key);
      if(basename.Contains(TRegexp(fix))){
	key=dirname+"/"+Replace(basename,fix,"");
	if(key.Contains("^[.]/")) key=Replace(key,"^[.]/","");
      }
    }
    key.ReplaceAll("forward","AFB");
    key.ReplaceAll("backward","AFB");
    key.ReplaceAll("AFB_num","weightedAFB");
    key.ReplaceAll("AFB_den","weightedAFB");
    histkeys_out.insert(key);
  }
  return histkeys_out;
}
set<TString> Plotter::GetParsedHistKeys(set<TString> excludes={}){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetParsedHistKeys(set<TString> excludes)]"<<endl;
  set<TString> histkeys;
  histkeys=GetHistKeys(entries[0]);
  set<TString> sys_fixes,sample_fixes;

  for(const auto& [sysname,sys]:systematics)
    for(const auto& suf:sys.suffixes)
      sys_fixes.insert(suf+"$");
  histkeys=ParseHistKeys(histkeys,sys_fixes,excludes);

  for(const auto& [samplename,sample]:samples)
    for(const auto& [filename,weight,prefix,suffix]:sample.files){
      if(prefix!="") sample_fixes.insert("^"+prefix);
      if(suffix!="") sample_fixes.insert(suffix+"$");
    }
  histkeys=ParseHistKeys(histkeys,sample_fixes);

  return histkeys;
}
*/
void Plotter::AddPlot(TString plotkey,TString plotoption){
  Plot plot;
  plot.name=plotkey;plot.histname=plotkey;
  if(entries.size()==1) plot.type=Plot::Type::Compare;
  plot.SetOption(plotoption);
  if(plot.type==Plot::Type::DoubleRatio){
    for(auto const& [k,p]:plots){
      if(k.Contains(TRegexp(plotkey+"$"))) plot.subplots.push_back(p);
    }
    if(plot.subplots.size()!=2){
      if(DEBUG>0) std::cout<<"###ERROR### DoubleRatio should have 2 subplots ("<<plot.subplots.size()<<"!=2"<<endl;
      plot.Print();
      exit(1);
    }
  }
  plots[plotkey]=plot;
  //SetPlotRebinXminXmaxAuto(key);
} 
/*
void Plotter::AddPlots(TRegexp regexp,TString plotoption){
  set<TString> histkeys=GetParsedHistKeys();
  for(const auto& key:histkeys){
    if(key.Contains(regexp)) AddPlot(key,plotoption);
  }
}
*/
void Plotter::RemovePlot(TString plotkey){
  plots.erase(plotkey);
}
void Plotter::RemovePlots(TRegexp regexp){
  set<TString> to_remove;
  for(const auto& [key,plot]:plots)
    if(key.Contains(regexp)) to_remove.insert(key);
  for(const auto& key:to_remove)
    RemovePlot(key);
}   
  
void Plotter::SetPlotRebinXminXmaxAuto(TString plotkey){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::SetPlotRebinXminXmaxAuto()]"<<endl;
  Plot& plot=plots[plotkey];
  TH1 *hh=GetHist(entries[0],plot.histname);
  double maxrelerr=0;
  int ifirst=-1;
  int ilast=hh->GetNbinsX();
  for(int i=hh->GetXaxis()->GetFirst();i<hh->GetXaxis()->GetLast()+1;i++){
    double val=hh->GetBinContent(i);
    if(val==0) continue;
    if(ifirst==-1) ifirst=i;
    ilast=i;
    double err=hh->GetBinError(i);
    double relerr=fabs(err/val);
    if(maxrelerr<relerr) maxrelerr=relerr;
  }
  plot.rebin=pow(maxrelerr/0.05,2)>10?10:pow(maxrelerr/0.05,2);
  if(ifirst!=1||ilast!=hh->GetNbinsX()){
    plot.xmin=hh->GetBinLowEdge(ifirst);
    plot.xmax=hh->GetBinLowEdge(ilast)+hh->GetBinWidth(ilast)*(1-0.0000001);
    plot.Print();
  }
}
void Plotter::SetPlotsRebinXminXmaxAuto(){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::SetPlotsRebinXminXmaxAuto()]"<<endl;
  for(auto& [plotname,plot]:plots) SetPlotRebinXminXmaxAuto(plotname);
}
void Plotter::SetPlotOption(TString plotname,TString option,TString suffix){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::SetPlotOption(TString plotname,TString option)]"<<endl;
  Plot temp(plots[plotname]);
  temp.SetOption(option);
  temp.SetOption(Form("name:%s",(temp.name+suffix).Data()));
  plots[plotname+suffix]=temp;
}
void Plotter::SetPlotsOption(TRegexp plotname_,TString option,TString suffix){
  vector<TString> list;
  for(auto& [plotname,plot]:plots){
    if(plotname.Contains(plotname_)) list.push_back(plotname); 
  }
  for(const auto& plotname:list) SetPlotOption(plotname,option,suffix);
}
void Plotter::RemovePlotOption(TString plotname,TString option){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::RemovePlotOption(TString plotname,TString option)]"<<endl;
  plots[plotname].RemoveOption(option);
}
void Plotter::RemovePlotsOption(TRegexp plotname_,TString option){
  for(auto& [plotname,plot]:plots){
    if(plotname.Contains(plotname_)) RemovePlotOption(plotname,option);
  }
}
void Plotter::SavePlotFile(){
  if(plotfile!=""&&plots.size()>0){
    ofstream f(plotdir+"/"+plotfile);
    for(const auto& [plotname,plot]:plots) plot.Print(f);
  }
  plotfile="";
}  

///////////////////////////////////////////////////////////////////////
///////////////////// MakeSysteamtic //////////////////////////////////
///////////////////////////////////////////////////////////////////////
Systematic Plotter::MakeSystematic(TString name,Systematic::Type type,int varibit,vector<TString> includes){
  Systematic this_sys=Systematic(name,type,varibit);
  if(this_sys.type==Systematic::Type::MULTI){
    this_sys.sysbit=0;
    for(const auto& [syskey,sys]:systematics)
      for(int j=0;j<includes.size();j++)
        if(syskey==includes[j]){
          this_sys.sysbit|=sys.sysbit;
          break;
        }
  }else{
    this_sys.sysbit=1<<systematics.size();
    this_sys.suffixes=includes;
  }
  //if(Verbosity) this_sys.Print();
  return this_sys;
}
Systematic Plotter::MakeSystematic(TString name,Systematic::Type type,int varibit,TString includes){
  return MakeSystematic(name,type,varibit,Split(includes," "));
}


///////////////////////////////////////////////////////////////////////
///////////////////// etc /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int Plotter::GetAutoColor() const {
  set<int> colors={1,2,3,4,5,6,7,8,9};
  for(const auto& color:colors){
    bool exist=false;
    for(const auto& entry:entries){
      if(entry.style.linecolor==color){
	exist=true;
	break;
      }
    }
    if(!exist){
      return color;
    }
  }
  return 0;
}
double Plotter::GetChi2(TH1* h1,TH1* h2) const {
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

//////////////////////working/////////////////////////
  

#endif













/*
TCanvas* Plotter::GetCompareAFBAll(vector<TString> histnames,int sysbit=0,TString option=""){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetCompareAFBAll(vector<TString> histnames,int sysbit=0,TString option="")]"<<endl;
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
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::GetCompareAFBAll(TRegexp regexp,int sysbit=0,TString option="")]"<<endl;
  vector<TString> histnames;
  for(auto it=plots.begin();it!=plots.end();it++){
    if(it->first.Contains(regexp)){
      histnames.push_back(it->second.name);
    }
  }
  return GetCompareAFBAll(histnames,sysbit,option);
}


void Plotter::SavePlots(TString outputdir="plot",int njob=1,int ijob=0){
  if(DEBUG>3) std::cout<<"###DEBUG### [Plotter::SavePlots(TString outputdir="plot",int njob=1,int ijob=0)]"<<endl;
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
    if(DEBUG>1) std::cout<<"###INFO#### [SavePlots] save "<<outputdir+"/"+histname+".png"<<endl;
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
      if(DEBUG>1) std::cout<<"###INFO#### [SavePlots] save "<<this_dir+"/"+histname(histname.Last('/'),histname.Length())+".png"<<endl;
      c=GetCompare(histname,systematics[is].sysbit,plot->rebin,plot->xmin,plot->xmax,plot->option);
      if(c){
	c->SaveAs(this_dir+"/"+histname(histname.Last('/'),histname.Length())+".png");
	THStack* hstack=(THStack*)c->GetPad(1)->GetPrimitive("hstack");
	if(hstack) hstack->GetHists()->Clear();
	delete c;
      }
      if(systematics[is].suffixes.size()==1){
	if(DEBUG>1) std::cout<<"###INFO#### [SavePlots] save "<<this_dir+"/"+histname(histname.Last('/'),histname.Length())+"_raw.png"<<endl;
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
