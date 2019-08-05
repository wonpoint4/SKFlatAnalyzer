#ifndef __PLOT_CC__
#define __PLOT_CC__
#include"Utils.h"
class Plot{
public:
  enum Type{UNDEFINE,Compare,Ratio,Diff,CompareAndRatio,CompareAndDiff,DoubleRatio,Collection};
  TString name;
  TString histname;
  TString sysname;
  Type type;
  int rebin=0;
  double xmin=0,xmax=0;
  TString option;
  vector<Plot> subplots;
  Plot(vector<TString> words);
  Plot(TString line=""):Plot(Split(line," ")){};
  ~Plot();
  void Print(std::ostream& out=cout) const;
  void SetOption(TString option_);
  void RemoveOption(TString option_);
};
void Plot::SetOption(TString option_){
  for(const auto& opt:Split(option_," ")){
    if(opt.Contains(TRegexp("^name:"))) name=opt(5,999);
    else if(opt.Contains(TRegexp("^histname:"))) histname=opt(9,999);
    else if(opt.Contains(TRegexp("^sysname:"))) sysname=opt(8,999);
    else if(opt.Contains(TRegexp("^type:"))) type=(Type)TString(opt(5,999)).Atoi();
    else if(opt.Contains(TRegexp("^rebin:"))) rebin=TString(opt(6,999)).Atoi();
    else if(opt.Contains(TRegexp("^xmin:"))) xmin=TString(opt(5,999)).Atof();
    else if(opt.Contains(TRegexp("^xmax:"))) xmax=TString(opt(5,999)).Atof();
    else option+=" "+opt;
  }
}   
void Plot::RemoveOption(TString option_){
  vector<TString> options=Split(option," ");
  vector<TString> options_remove=Split(option_," ");
  for(const auto& remove:options_remove){
    if(remove=="rebin") rebin=0;
    else if(remove=="xmin") xmin=0;
    else if(remove=="xmax") xmax=0;
    else{
      for(int i=0;i<options.size();i++){
        if(options[i].Contains(TRegexp("^"+remove))){
          options.erase(options.begin()+i);
          i--;
        }
      }
    }
  }
  TString newoption;
  for(const auto& opt:options) newoption+=opt+" ";
  option=newoption;
}
void Plot::Print(std::ostream& out) const{
  if(DEBUG>3) cout<<"###DEBUG### [void Plot::Print(std::ostream& out) const]"<<endl;
  out<<"<Plot> ";
  out<<"name:"<<name<<" histname:"<<histname;
  if(sysname!="") out<<" sysname:"<<sysname;
  if(type!=Type::UNDEFINE) out<<" type:"<<type;
  out<<" rebin:"<<rebin<<" xmin:"<<xmin<<" xmax:"<<xmax<<" "<<option<<" </Plot>"<<endl;
}
Plot::~Plot(){
  if(DEBUG>3) cout<<"###DEBUG### [Plot::~Plot()]"<<endl;
}
Plot::Plot(vector<TString> words){
  if(DEBUG>3) cout<<"###DEBUG### [Plot::Plot(vector<TString> words)]"<<endl;
  int imax=words.size();
  if(imax==0) return;
  if(words[0]!="<Plot>") cout<<"###ERROR### [Plot::Plot] Wrong format"<<endl;
  int depth=1;
  for(int i=1;i<imax;i++){
    if(words[i]=="<Plot>") subplots.push_back(Plot(vector<TString>(words.begin()+i,words.end())));
    else if(words[i]=="</Plot>") depth--;
    else SetOption(words[i]);
    if(depth==0) break;
  }
  this->Print();
}  
Plot::Plot(TString line="") : Plot(Split(line," ")){
  if(DEBUG>3) cout<<"###DEBUG### [Plot::Plot(TString line="")]"<<endl;
}

#endif
