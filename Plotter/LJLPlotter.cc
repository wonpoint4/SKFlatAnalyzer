#include"Plotter.cc"
#ifndef __LJLPLOTTER_CC__
#define __LJLPLOTTER_CC__
class LJLPlotter:public Plotter{
public:
  void SetupEntries();
  void SetupSystematics();
  int Setup(TString mode_="");
  TString mode;
  TString analyzer;
  LJLPlotter(TString mode_="");
  double GetChi2(TH1* h1,TH1* h2=NULL);
};
LJLPlotter::LJLPlotter(TString mode_=""){
  analyzer="LJLAnalyzer";
  vector<TString> files=Split(gSystem->GetFromPipe("find $SKFlatOutputDir$SKFlatV/"+analyzer+"/ -type f"),"\n");
  for(auto file:files){
    TString key=Replace(file,TString()+getenv("SKFlatOutputDir")+getenv("SKFlatV")+"/"+analyzer+"/","");
    key.ReplaceAll(".root","");
    Sample sample;
    sample.files.push_back(make_tuple(file,1.,"",""));
    samples[key]=sample;
  } 
  
  vector<TString> syears={"2016","2017","2018"};
  for(const auto& syear:syears){
    TString skimprefix="SkimTree_LJL_";
    samples["muon"+syear]=Sample("muon"+syear,Sample::Type::DATA,kBlack)+TRegexp(syear+"/DATA/"+analyzer+"_"+skimprefix+"SingleMuon_[A-Z].*");
    samples["electron"+syear]=Sample("electron"+syear,Sample::Type::DATA,kBlack)+TRegexp(syear+"/DATA/"+analyzer+"_"+skimprefix+".*EG.*_[A-Z].*");

    samples["amc"+syear]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp(syear+"/"+analyzer+"_"+skimprefix+"DYJets$");
    samples["amctt"+syear]="tau_"%(Sample("#gamma*/Z#rightarrow#tau#tau",Sample::Type::SIGNAL,kGreen)+TRegexp(syear+"/"+analyzer+"_"+skimprefix+"DYJets$"));
    samples["vv"+syear]=Sample("Diboson",Sample::Type::BG,kBlue)+TRegexp(syear+"/"+analyzer+"_"+skimprefix+"[W-Z][W-Z]_pythia$");
    samples["wjets"+syear]=Sample("W",Sample::Type::BG,kYellow)+TRegexp(syear+"/"+analyzer+"_"+skimprefix+"WJets_MG$");
    samples["tt"+syear]=Sample("t#bar{t}",Sample::Type::BG,kMagenta)+TRegexp(syear+"/"+analyzer+"_"+skimprefix+"TTLL_powheg$");

  }
  samples["data"]=Sample("data",Sample::Type::DATA,kBlack)+TRegexp("^muon201[6-8]$")+TRegexp("^electron201[6-8]$");
  samples["sim"]=Sample("sim",Sample::Type::STACK,kRed)+TRegexp("^amc201[6-8]$")+TRegexp("^amctt201[6-8]$")+TRegexp("^vv201[6-8]$")+TRegexp("^wjets201[6-8]$")+TRegexp("^tt201[6-8]$");
  samples["qcd"]="pp_"%(Sample("qcd",Sample::Type::BG,kCyan)+"data"+(-1)*(Sample()+TRegexp("^amc201[6-8]$")+TRegexp("^amctt201[6-8]$")+TRegexp("^vv201[6-8]$")+TRegexp("^wjets201[6-8]$")+TRegexp("^tt201[6-8]$")))+"mm_"%(Sample("qcd",Sample::Type::BG,kCyan)+"data"+(-1)*(Sample()+TRegexp("^amc201[6-8]$")+TRegexp("^amctt201[6-8]$")+TRegexp("^vv201[6-8]$")+TRegexp("^wjets201[6-8]$")+TRegexp("^tt201[6-8]$")));
  Setup(mode_);
}

int LJLPlotter::Setup(TString mode_=""){
  entries.clear();
  systematics.clear();
  plots.clear();
  mode=mode_;

  SetupEntries();
  SetupSystematics();
  SetupPlots("plots_LJLAnalyzer/"+mode+"/"+mode+".dat");

  if(DEBUG) cout<<"[Setup] nentries: "<<entries.size()<<endl;
  if(DEBUG) cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}

void LJLPlotter::SetupEntries(){
  if(DEBUG)  cout<<"[LJLPlotter::SetupEntries]"<<endl;
  if(mode==""){
    entries.push_back(samples["data"]);
    entries.push_back(samples["sim"]);      
  }else if(mode=="qcd"){
    entries.push_back(samples["data"]);
    entries.push_back(samples["sim"]+"qcd");      
  }
  return;
}
void LJLPlotter::SetupSystematics(){
  if(DEBUG)  cout<<"[SetupSystematics]"<<endl;
  systematics["RECOSF"]=MakeSystematic("RECOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_RECOSF_up _RECOSF_down");
  systematics["IDSF"]=MakeSystematic("IDSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_IDSF_up _IDSF_down");
  systematics["ISOSF"]=MakeSystematic("ISOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_ISOSF_up _ISOSF_down");
  systematics["triggerSF"]=MakeSystematic("triggerSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_triggerSF_up _triggerSF_down");

  systematics["noRECOSF"]=MakeSystematic("noRECOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noRECOSF");
  systematics["noIDSF"]=MakeSystematic("noIDSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noIDSF");
  systematics["noISOSF"]=MakeSystematic("noISOSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noISOSF");
  systematics["notriggerSF"]=MakeSystematic("notriggerSF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_notriggerSF");

  systematics["nozptcor"]=MakeSystematic("nozptcor",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_nozptcor");
  systematics["noefficiencySF"]=MakeSystematic("noefficiencySF",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL)+(1<<Sample::Type::BG),"_noefficiencySF");
  systematics["efficiencySF"]=MakeSystematic("efficiencySF",Systematic::Type::MULTI,0,"RECOSF IDSF ISOSF triggerSF");

}

double LJLPlotter::GetChi2(TH1* h1,TH1* h2){
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
  
#endif
