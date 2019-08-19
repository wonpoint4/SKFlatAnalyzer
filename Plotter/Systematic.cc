#ifndef __SYSTEMATIC_CC__
#define __SYSTEMATIC_CC__
#include "Utils.h"
class Systematic{
public:
  enum Type{UNDEFINE,ENVELOPE,GAUSSIAN,HESSIAN,MULTI};
  TString name;
  Type type=Type::UNDEFINE;
  vector<TString> suffixes;
  int sysbit=0;
  int varibit=0;

  Systematic(TString name_="",Type type_=Type::UNDEFINE,int varibit_=0);
  ~Systematic();
  TString GetTypeString() const;
  void Print() const;
};
Systematic::Systematic(TString name_,Type type_,int varibit_){
  name=name_;
  type=type_;
  varibit=varibit_;
}
Systematic::~Systematic(){}

TString Systematic::GetTypeString() const{
  switch(type){
  case UNDEFINE: return "UNDEFINE";
  case ENVELOPE: return "ENVELOPE";
  case GAUSSIAN: return "GAUSSIAN";
  case HESSIAN: return "HESSIAN";
  case MULTI: return "MULTI";
  default: return "###ERROR### Bad Systematic::Type";
  }
}
void Systematic::Print() const{
  cout<<"#######################"<<endl;
  cout<<name<<" "<<GetTypeString()<<endl;
  cout<<"sysbit:"<<sysbit<<" varibit:"<<varibit<<endl;
  for(const auto& suffix:suffixes) cout<<suffix<<" ";
  cout<<"\n"<<endl;
}
#endif
