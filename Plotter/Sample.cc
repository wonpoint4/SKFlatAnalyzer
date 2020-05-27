#ifndef __SAMPLE_CC__
#define __SAMPLE_CC__
#include"Global.h"
#include"Style.cc"
class Sample{
public:
  enum Type{UNDEF,DATA,SIGNAL,BG,GEN,STACK,SUM,A,B,C,D,FILE};
  TString GetTypeString() const{
    switch(type){
    case UNDEF: return "UNDEF";
    case DATA: return "DATA";
    case SIGNAL: return "SIGNAL";
    case BG: return "BG";
    case GEN: return "GEN";
    case STACK: return "STACK";
    case SUM: return "SUM";
    case A: return "A";
    case B: return "B";
    case C: return "C";
    case D: return "D";
    case FILE: return "FILE";
    default: return "###ERROR### Bad Sample::Type";
    }
  }
  TString title;
  Type type;
  Style style;
  Style style_alt;
  map<TString,TString> replace;
  double weight=1.;
  vector<Sample> subs;

  Sample(TString title_="",Sample::Type type_=Sample::Type::UNDEF,Style styleA=Style(),Style styleB=Style());
  Sample(TString title_,Sample::Type type_,int color,int marker=-1,int fill=-1,TString drawoption="");
  Sample operator+(const Sample& sam);
  Sample operator+(const char* key);
  Sample operator+(const TString key);
  Sample operator+(const TRegexp reg);
  Sample operator-(const Sample& sam);
  Sample operator-(const char* key);
  Sample operator-(const TString key);
  Sample operator-(const TRegexp reg);
  friend Sample operator%(const char* prefix,const Sample& sam);
  friend Sample operator%(const TString prefix,const Sample& sam);
  friend Sample operator%(const Sample& sam,const char* suffix);
  friend Sample operator%(const Sample& sam,const TString suffix);

  Sample operator*(double w) const;
  friend Sample operator*(double w,const Sample& sam);
  void SetStyle(int color,int marker=-1,int fill=-1,TString drawoption="");
  void SetType(int type_);
  void Add(TRegexp sampleregexp,double weight=1.,TString prefix="",TString suffix="");
  void ApplyStyle(TH1* hist,bool alt=false) const;
  void Print(bool detail=false,TString prefix="") const;
  TString ReplaceToString() const;
  bool IsCollection() const;
  bool IsSample() const;
  bool IsFile() const;
  static bool CheckAttributes(const Sample& sam1, const Sample& sam2);
};
map<TString,Sample> samples;

Sample::Sample(TString title_,Sample::Type type_,Style styleA,Style styleB){
  title=title_;
  type=type_;
  style=styleA;
  style_alt=styleB;
}
Sample::Sample(TString title_,Sample::Type type_,int color,int marker,int fill,TString drawoption){
  title=title_;
  type=type_;
  SetStyle(color,marker,fill,drawoption);
}
TString Sample::ReplaceToString() const {
  TString out;
  for(const auto& [reg,newstr]:replace){
    out+="{'"+reg+"'->'"+newstr+"'} ";
  }
  return out;
}
bool Sample::IsCollection() const {
  if(type==Sample::Type::STACK||type==Sample::Type::SUM) return true;
  else return false;
}
bool Sample::IsFile() const {
  if(type==Sample::Type::FILE) return true;
  return false;
}
bool Sample::IsSample() const{
  if(!IsCollection()&&!IsFile()) return true;
  return false;
}
bool Sample::CheckAttributes(const Sample& sam1,const Sample& sam2){
  bool ret=true;
  if(sam1.weight!=sam2.weight){
    PWarning("[Sample::CheckAttributes] "+sam1.title+Form("%f",sam1.weight)+" != "+sam2.title+Form("%f",sam2.weight));
    ret=false;
  }
  if(sam1.replace!=sam2.replace){
    PWarning("[Sample::CheckAttributes] "+sam1.title+" "+sam1.ReplaceToString()+" != "+sam2.title+" "+sam2.ReplaceToString());
    ret=false;
  }
  return ret;
}
Sample Sample::operator+(const Sample& sam){
  Sample temp(*this);
  if(temp.IsCollection()&&sam.IsCollection()){
    if(!CheckAttributes(temp,sam))
      PWarning("Some attributes of "+sam.title+" will be ignored");
    temp.subs.insert(temp.subs.end(),sam.subs.begin(),sam.subs.end());
  }else if(temp.IsCollection()&&sam.IsSample()){
    temp.subs.push_back(sam);
    if(temp.type==Sample::Type::STACK) temp.subs.back().style.fillcolor=temp.subs.back().style.linecolor;
  }else if(temp.IsCollection()&&sam.IsFile()){
    Sample sample(sam.title,Type::UNDEF);
    sample=sample+sam;
    temp=temp+sample;
  }else if(temp.IsSample()&&sam.IsSample()){
    if(!CheckAttributes(temp,sam))
      PWarning("Some attributes of "+sam.title+" will be ignored");
    temp.subs.insert(temp.subs.end(),sam.subs.begin(),sam.subs.end());
  }else if(temp.IsSample()&&sam.IsFile()){ 
    temp.subs.push_back(sam);
  }else{
    PError("[Sample::operator+] Cannot add samples; "+temp.title+"("+temp.GetTypeString()+") + "+sam.title+"("+sam.GetTypeString()+")");
  }
  return temp;
}
Sample Sample::operator+(const char* key){
  for(const auto& [k,sample]:samples){
    if(k==key) return (*this)+sample;
  }
  PError((TString)"[Sample::operator+(const char* key,const Sample& sam)] no sample with key "+key);
  exit(1);
}  
Sample Sample::operator+(const TString key){
  return (*this)+key.Data();
}
Sample Sample::operator+(const TRegexp reg){
  Sample temp(*this);
  for(auto const& [key,sample]:samples)
    if(key.Contains(reg)) temp=temp+sample;
  return temp;
}

Sample Sample::operator-(const Sample& sam){
  return (*this)+(-1.)*sam;
}
Sample Sample::operator-(const char* key){
  for(const auto& [k,sample]:samples){
    if(k==key) return (*this)-sample;
  }
  PError((TString)"[Sample operator-(const char* key,const Sample& sam)] no sample with key "+key);
  exit(1);
}  
Sample Sample::operator-(const TString key){
  return (*this)-key.Data();
}
Sample Sample::operator-(const TRegexp reg){
  Sample temp(*this);
  for(auto const& [key,sample]:samples)
    if(key.Contains(reg)) temp=temp-sample;
  return temp;
}
Sample operator%(const char* prefix,const Sample& sam){
  Sample temp(sam);
  if(temp.replace.find("/([^/]*)$")==temp.replace.end())
    temp.replace["/([^/]*)$"]="/"+TString(prefix)+"$1";
  else if(TPRegexp("/(.*)\\$1$").MatchB(temp.replace["/([^/]*)$"])){
    TObjArray* array=TPRegexp("/(.*)\\$1$").MatchS(temp.replace["/([^/]*)$"]);
    TString old_prefix=((TObjString*)array->At(1))->GetString();
    temp.replace["/([^/]*)$"]="/"+TString(prefix)+old_prefix+"$1";
    array->Delete();
  }else{
    PWarning("[Sample::operator%] Cannot add prefix... Overwrite...");
    temp.replace["/([^/]*)$"]="/"+TString(prefix)+"$1";    
  }
  return temp;
}  
Sample operator%(const TString prefix,const Sample& sam){
  return prefix.Data()%sam;
}  
Sample operator%(const Sample& sam,const char* suffix){
  Sample temp(sam);
  temp.replace["$"]=temp.replace["$"]+TString(suffix);
  return temp;
}  
Sample operator%(const Sample& sam,const TString suffix){
  return sam%suffix.Data();
}  
Sample Sample::operator*(double f) const {
  Sample temp(*this);
  temp.weight*=f;
  return temp;
}
Sample operator*(double f,const Sample& sam){
  return sam*f;
}
void Sample::SetStyle(int color,int marker,int fill,TString drawoption){
  style=Style(color,marker,fill,drawoption);
  
  if(type==Type::STACK){
    style.drawoption="e";
  }
}
void Sample::SetType(int type_){
  if(type==Type::STACK){
    for(auto& sub:subs){
      sub.SetType(type_);
    }
  }else type=(Type)type_;
}
void Sample::Add(TRegexp regexp,double weight,TString prefix,TString suffix){
  for(auto it=samples.begin();it!=samples.end();it++){
    if(it->first.Contains(regexp)) (*this)=(*this)+weight*(prefix%(it->second)%suffix);
  }
}
void Sample::ApplyStyle(TH1* hist,bool alt) const {
  PAll("[Sample::ApplyStyle(TH1* hist)]");
  if(alt) style_alt.Apply(hist);
  else style.Apply(hist);

  if(hist){
    hist->SetName(title);
  }
}
void Sample::Print(bool detail,TString pre) const{
  if(type==Type::FILE){
    cout<<pre<<title<<" "<<gSystem->GetFromPipe("stat -c %.19y "+title)<<endl;
    return;
  }
  if(detail){
    TString weightstring="";
    if(weight==1) weightstring="+";
    else if(weight==-1) weightstring="-";
    else weightstring=Form("%+.1f ",weight);
    cout<<pre<<weightstring<<"Title:"<<title<<" Type:"<<GetTypeString()<<" "<<ReplaceToString()<<" ";
    style.Print();
    for(const auto& sub:subs){
      sub.Print(detail,pre+"  ");
    }
  }else cout<<pre<<"Title: "<<title<<" "<<"Type: "<<GetTypeString()<<endl;
}
#endif
