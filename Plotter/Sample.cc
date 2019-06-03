#ifndef __SAMPLE_CC__
#define __SAMPLE_CC__
#include"Global.h"
class SampleFrag{
public:
  enum Type{DATA,SIGNAL,BG,GEN,SUM,STACK,SYS,UNDEFINE};
  TString GetTypeString() const{
    switch(type){
    case DATA: return "DATA";
    case SIGNAL: return "SIGNAL";
    case BG: return "BG";
    case GEN: return "GEN";
    case SUM: return "SUM";
    case STACK: return "STACK";
    case SYS: return "SYS";
    case UNDEFINE: return "UNDEFINE";
    default: return "###ERROR### Bad SampleFrag::Type";
    }
  }
  TString title;
  Type type;
  int fillcolor;
  int fillstyle;
  int linecolor;
  int linestyle;
  int linewidth;
  int markercolor;
  int markerstyle;
  double markersize;
  vector<tuple<TString,double,TString,TString>> files; //filename,weight,prefix,suffix

  SampleFrag();
  void SetColor(int color);
  virtual void Add(TString samplefragkey,double weight=1.,TString prefix="",TString suffix="");
  virtual void Add(TRegexp samplefragkeyexp,double weight=1.,TString prefix="",TString suffix="");
  void ApplyHistStyle(TH1* hist) const ;
  virtual void Print() const;
};
map<TString,SampleFrag> samplefrags;

SampleFrag::SampleFrag(){
  fillcolor=-1;
  fillstyle=-1;
  linecolor=-1;
  linestyle=-1;
  linewidth=-1.;
  markercolor=-1;
  markerstyle=-1;
  markersize=-1.;
  type=(SampleFrag::Type)UNDEFINE;
}
SampleFrag MakeSampleFrag(TString title,SampleFrag::Type type,int color,vector<tuple<TString,double,TString,TString>> files){
  SampleFrag samplefrag;
  samplefrag.title=title;samplefrag.type=type;
  samplefrag.SetColor(color);
  samplefrag.files=files;
  return samplefrag;
}
SampleFrag MakeSampleFrag(TString title,SampleFrag::Type type,int color,tuple<TString,double,TString,TString> file1=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file2=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file3=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file4=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file5=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file6=make_tuple("",0.,"",""),tuple<TString,double,TString,TString> file7=make_tuple("",0.,"","")){
  vector<tuple<TString,double,TString,TString>> files;
  if(get<0>(file1)!="") files.push_back(file1);
  if(get<0>(file2)!="") files.push_back(file2);
  if(get<0>(file3)!="") files.push_back(file3);
  if(get<0>(file4)!="") files.push_back(file4);
  if(get<0>(file5)!="") files.push_back(file5);
  if(get<0>(file6)!="") files.push_back(file6);
  if(get<0>(file7)!="") files.push_back(file7);
  return MakeSampleFrag(title,type,color,files);
}
SampleFrag MakeSampleFrag(TString title,SampleFrag::Type type,int color,tuple<TString,double> frag1,tuple<TString,double> frag2=make_tuple("",0.),tuple<TString,double> frag3=make_tuple("",0.),tuple<TString,double> frag4=make_tuple("",0.),tuple<TString,double> frag5=make_tuple("",0.),tuple<TString,double> frag6=make_tuple("",0.),tuple<TString,double> frag7=make_tuple("",0.)){
  SampleFrag frag=MakeSampleFrag(title,type,color);
  if(get<0>(frag1)!="") frag.Add(get<0>(frag1),get<1>(frag1),"","");
  if(get<0>(frag2)!="") frag.Add(get<0>(frag2),get<1>(frag2),"","");
  if(get<0>(frag3)!="") frag.Add(get<0>(frag3),get<1>(frag3),"","");
  if(get<0>(frag4)!="") frag.Add(get<0>(frag4),get<1>(frag4),"","");
  if(get<0>(frag5)!="") frag.Add(get<0>(frag5),get<1>(frag5),"","");
  if(get<0>(frag6)!="") frag.Add(get<0>(frag6),get<1>(frag6),"","");
  if(get<0>(frag7)!="") frag.Add(get<0>(frag7),get<1>(frag7),"","");
  return frag;
}
SampleFrag MakeSampleFrag(TString title,SampleFrag::Type type,int color,TString subsamplefragkey,double weight=1.,TString prefix="",TString suffix=""){
  SampleFrag samplefrag=MakeSampleFrag(title,type,color);
  samplefrag.Add(subsamplefragkey,weight,prefix,suffix);
  return samplefrag;
}
SampleFrag MakeSampleFrag(TString title,SampleFrag::Type type,int color,TRegexp subsamplefragregexp,double weight=1.,TString prefix="",TString suffix=""){
  SampleFrag samplefrag=MakeSampleFrag(title,type,color);
  samplefrag.Add(subsamplefragregexp,weight,prefix,suffix);
  return samplefrag;
}
void SampleFrag::SetColor(int color){
  linecolor=color;
  markercolor=color;
}
void SampleFrag::Add(TString samplefragkey,double weight,TString prefix,TString suffix){
  auto it=samplefrags.find(samplefragkey);
  if(it!=samplefrags.end()){
    auto tempfiles=it->second.files;
    for(unsigned int i=0;i<tempfiles.size();i++){
      get<1>(tempfiles[i])=get<1>(tempfiles[i])*weight;
      get<2>(tempfiles[i])=prefix+get<2>(tempfiles[i]);
      get<3>(tempfiles[i])=get<3>(tempfiles[i])+suffix;
    }
    files.insert(files.end(),tempfiles.begin(),tempfiles.end());
  }else{
    if(DEBUG>0) cout<<"###WARNING### [SampleFrag::Add] no samplefrag "<<samplefragkey<<endl;
  }
  return;
}
void SampleFrag::Add(TRegexp samplefragkeyexp,double weight,TString prefix,TString suffix){
  for(auto it=samplefrags.begin();it!=samplefrags.end();it++){
    if(it->first.Contains(samplefragkeyexp)) Add(it->first,weight,prefix,suffix);
  }
}
void SampleFrag::ApplyHistStyle(TH1* hist) const {
  if(hist){
    if(type!=Type::STACK){
      if(linecolor>=0) hist->SetLineColor(linecolor);
      if(linestyle>=0) hist->SetLineStyle(linestyle);
      if(linewidth>=0) hist->SetLineWidth(linewidth);
      if(fillcolor>=0) hist->SetFillColor(fillcolor);
      if(fillstyle>=0) hist->SetFillStyle(fillstyle);
      if(markercolor>=0) hist->SetMarkerColor(markercolor);
      if(markerstyle>=0) hist->SetMarkerStyle(markerstyle);
      if(markersize>=0) hist->SetMarkerSize(markersize);
    }
    hist->SetNameTitle(title,title);
  }
}
void SampleFrag::Print() const{
  cout<<"Title: "<<title<<" "<<"Type: "<<GetTypeString()<<endl;
  cout<<"Style: "<<fillcolor<<" "<<fillstyle<<" "<<linecolor<<" "<<linestyle<<" "<<linewidth<<" "<<markercolor<<" "<<markerstyle<<" "<<markersize<<endl;
  for(auto it=files.begin();it!=files.end();it++){
    cout<<get<0>(*it)<<" "<<get<1>(*it)<<" "<<get<2>(*it)<<" "<<get<3>(*it)<<endl;
  }
}

class Sample:public SampleFrag{
public:
  vector<tuple<SampleFrag,double>> frags;
  Sample();
  void Add(TString samplefragkey,double weight);
  void Print(bool detail=false) const;
};
Sample::Sample(){
}
Sample MakeSample(TString title,Sample::Type type,int color){
  Sample sample;
  sample.title=title;sample.type=type;
  sample.SetColor(color);
  return sample;
}
Sample MakeSample(TString title,Sample::Type type,int color,tuple<TString,double> frag1,tuple<TString,double> frag2=make_tuple("",0.),tuple<TString,double> frag3=make_tuple("",0.),tuple<TString,double> frag4=make_tuple("",0.),tuple<TString,double> frag5=make_tuple("",0.),tuple<TString,double> frag6=make_tuple("",0.),tuple<TString,double> frag7=make_tuple("",0.)){
  Sample sample=MakeSample(title,type,color);
  if(get<0>(frag1)!="") sample.Add(get<0>(frag1),get<1>(frag1));
  if(get<0>(frag2)!="") sample.Add(get<0>(frag2),get<1>(frag2));
  if(get<0>(frag3)!="") sample.Add(get<0>(frag3),get<1>(frag3));
  if(get<0>(frag4)!="") sample.Add(get<0>(frag4),get<1>(frag4));
  if(get<0>(frag5)!="") sample.Add(get<0>(frag5),get<1>(frag5));
  if(get<0>(frag6)!="") sample.Add(get<0>(frag6),get<1>(frag6));
  if(get<0>(frag7)!="") sample.Add(get<0>(frag7),get<1>(frag7));
  if(type==Sample::Type::SYS){
    for(auto& it:sample.frags){
      get<0>(it).type=SampleFrag::Type::SYS;
    }
  }
  return sample;
}
void Sample::Add(TString samplefragkey,double weight){
  auto it=samplefrags.find(samplefragkey);
  if(it!=samplefrags.end()){
    frags.push_back(make_tuple(it->second,weight));
    if(type==Type::STACK) get<0>(frags.back()).fillcolor=get<0>(frags.back()).linecolor;
  }else{
    if(DEBUG>0) cout<<"###WARNING### [Sample::Add] no samplefrag "<<samplefragkey<<endl;
  }
  return;
}
void Sample::Print(bool detail) const{
  cout<<"Title: "<<title<<" "<<"Type: "<<GetTypeString()<<endl;
  cout<<"Style: "<<fillcolor<<" "<<fillstyle<<" "<<linecolor<<" "<<linestyle<<" "<<linewidth<<" "<<markercolor<<" "<<markerstyle<<" "<<markersize<<endl;
  for(auto it=frags.begin();it!=frags.end();it++){
    cout<<get<0>(*it).title<<" "<<get<1>(*it)<<endl;
    if(detail) get<0>(*it).Print();
  }
}

#endif
