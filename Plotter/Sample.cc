#include"Global.h"
class SampleFrag{
public:
  enum Type{DATA,SIGNAL,BG,COLL,STACK};
  TString GetTypeString(){
    switch(type){
    case DATA: return "DATA";
    case SIGNAL: return "SIGNAL";
    case BG: return "BG";
    case COLL: return "COLL";
    case STACK: return "STACK";
    default: return "###ERROR### Bad Sample::Type";
    }
  }
  TString title;
  Type type;
  int fillcolor;
  int linecolor;
  int markercolor;
  int markerstyle;
  double markersize;
  vector<tuple<TString,TString,double>> files; //filename,prefix,weight
  void DefaultSetup(int color);
  virtual void Add(TString samplefragkey,double weight);
  virtual void Add(TRegexp samplefragkeyexp,double weight);
};
map<TString,SampleFrag> samplefrags;

void SampleFrag::DefaultSetup(int color){
  fillcolor=color;
  linecolor=color;
  markercolor=color;
  if(type==1){
    markerstyle=20;
    markersize=0.7;
  }
}
void SampleFrag::Add(TString samplefragkey,double weight){
  auto it=samplefrags.find(samplefragkey);
  if(it!=samplefrags.end()){
    auto tempfiles=it->second.files;
    for(unsigned int i=tempfiles.size();i<tempfiles.size();i++){
      get<2>(files[i])=get<2>(files[i])*weight;
    }
    files.insert(files.end(),tempfiles.begin(),tempfiles.end());
  }else{
    if(DEBUG>0) cout<<"###WARNING### [SampleFrag::Add] no samplefrag "<<samplefragkey<<endl;
  }
  return;
}
void SampleFrag::Add(TRegexp samplefragkeyexp,double weight){
  for(auto it=samplefrags.begin();it!=samplefrags.end();it++){
    if(it->first.Contains(samplefragkeyexp)) Add(it->first,weight);
  }
}

class SampleFrag SampleFrag(TString title,SampleFrag::Type type,int color,vector<tuple<TString,TString,double>> files){
  class SampleFrag samplefrag;
  samplefrag.title=title;samplefrag.type=type;
  samplefrag.DefaultSetup(color);
  samplefrag.files=files;
  return samplefrag;
}
class SampleFrag SampleFrag(TString title,SampleFrag::Type type,int color,tuple<TString,TString,double> file1=make_tuple("","",0.),tuple<TString,TString,double> file2=make_tuple("","",0.),tuple<TString,TString,double> file3=make_tuple("","",0.),tuple<TString,TString,double> file4=make_tuple("","",0.),tuple<TString,TString,double> file5=make_tuple("","",0.),tuple<TString,TString,double> file6=make_tuple("","",0.),tuple<TString,TString,double> file7=make_tuple("","",0.)){
  vector<tuple<TString,TString,double>> files;
  if(get<0>(file1)!="") files.push_back(file1);
  if(get<0>(file2)!="") files.push_back(file2);
  if(get<0>(file3)!="") files.push_back(file3);
  if(get<0>(file4)!="") files.push_back(file4);
  if(get<0>(file5)!="") files.push_back(file5);
  if(get<0>(file6)!="") files.push_back(file6);
  if(get<0>(file7)!="") files.push_back(file7);
  return SampleFrag(title,type,color,files);
}
class SampleFrag SampleFrag(TString title,SampleFrag::Type type,int color,tuple<TString,double> frag1,tuple<TString,double> frag2=make_tuple("",0.),tuple<TString,double> frag3=make_tuple("",0.),tuple<TString,double> frag4=make_tuple("",0.),tuple<TString,double> frag5=make_tuple("",0.),tuple<TString,double> frag6=make_tuple("",0.),tuple<TString,double> frag7=make_tuple("",0.)){
  class SampleFrag frag=SampleFrag(title,type,color);
  if(get<0>(frag1)!="") frag.Add(get<0>(frag1),get<1>(frag1));
  if(get<0>(frag2)!="") frag.Add(get<0>(frag2),get<1>(frag2));
  if(get<0>(frag3)!="") frag.Add(get<0>(frag3),get<1>(frag3));
  if(get<0>(frag4)!="") frag.Add(get<0>(frag4),get<1>(frag4));
  if(get<0>(frag5)!="") frag.Add(get<0>(frag5),get<1>(frag5));
  if(get<0>(frag6)!="") frag.Add(get<0>(frag6),get<1>(frag6));
  if(get<0>(frag7)!="") frag.Add(get<0>(frag7),get<1>(frag7));
  return frag;
}




class Sample:public SampleFrag{
public:
  vector<tuple<SampleFrag,double>> frags;
  virtual void Add(TString samplefragkey,double weight);
};
void Sample::Add(TString samplefragkey,double weight){
  auto it=samplefrags.find(samplefragkey);
  if(it!=samplefrags.end()){
    frags.push_back(make_tuple(it->second,weight));
  }else{
    if(DEBUG>0) cout<<"###WARNING### [Sample::Add] no samplefrag "<<samplefragkey<<endl;
  }
  return;
}
class Sample Sample(TString title,Sample::Type type,int color){
  class Sample sample;
  sample.title=title;sample.type=type;
  sample.DefaultSetup(color);
  return sample;
}
class Sample Sample(TString title,Sample::Type type,int color,tuple<TString,double> frag1,tuple<TString,double> frag2=make_tuple("",0.),tuple<TString,double> frag3=make_tuple("",0.),tuple<TString,double> frag4=make_tuple("",0.),tuple<TString,double> frag5=make_tuple("",0.),tuple<TString,double> frag6=make_tuple("",0.),tuple<TString,double> frag7=make_tuple("",0.)){
  class Sample sample=Sample(title,type,color);
  if(get<0>(frag1)!="") sample.Add(get<0>(frag1),get<1>(frag1));
  if(get<0>(frag2)!="") sample.Add(get<0>(frag2),get<1>(frag2));
  if(get<0>(frag3)!="") sample.Add(get<0>(frag3),get<1>(frag3));
  if(get<0>(frag4)!="") sample.Add(get<0>(frag4),get<1>(frag4));
  if(get<0>(frag5)!="") sample.Add(get<0>(frag5),get<1>(frag5));
  if(get<0>(frag6)!="") sample.Add(get<0>(frag6),get<1>(frag6));
  if(get<0>(frag7)!="") sample.Add(get<0>(frag7),get<1>(frag7));
  return sample;
}
