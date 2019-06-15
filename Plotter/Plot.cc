#ifndef __PLOT_CC__
#define __PLOT_CC__
class Plot{
public:
  enum Type={Compare,Ratio,Diff,CompareAndRatio,CompareAndDiff};
  TString name;
  TString histname;
  TString sysname;
  Type type;
  int rebin;
  double xmin,xmax;
  TString option;
  void Print(std::ostream& out=cout);
};
void Print(std::ostream& out){
  out<<"name:"<<name<<endl;
  out<<"histname:"<<histname<<" sysname:"<<sysname<<" type:"<<type<<endl;
  out<<"rebin:"<<rebin<<" xmin:"<<xmin<<" xmax:"<<xmax<<endl;
  out<<"option:"<<endl;
}
#endif
