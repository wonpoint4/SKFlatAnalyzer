#include"Plotter.cc"
#ifndef __LQ766PLOTTER_CC__
#define __LQ766PLOTTER_CC__
class LQ766Plotter:public Plotter{
public:
  void SetupEntries();
  void SetupSystematics();
  int Setup(TString mode_);
  TString mode;
  LQ766Plotter(TString mode_="data sim_stack");
  double GetChi2(TH1* h1,TH1* h2=NULL);
};
LQ766Plotter::LQ766Plotter(TString mode_){
  vector<TString> files=Split(gSystem->GetFromPipe("find ${LQANALYZER_OUTPUT_PATH}CAT/ExampleAnalyzerDiElectron/periodCtoD/ -type f"),"\n");
  for(int i=0;i<files.size();i++){
    TString file=files[i];
    TString key=Replace(file,TString()+getenv("LQANALYZER_OUTPUT_PATH")+"CAT/ExampleAnalyzerDiElectron/periodCtoD/","");
    Sample sample(key,Sample::Type::UNDEF,i%8+2);
    sample.files.push_back(make_tuple(file,1.,"",""));
    samples[key]=sample;
  } 

  samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("data_DoubleEG");
  samples["amc"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("SKDY50plus");
  samples["amctt"]="tau_"%(Sample("#gamma*/Z#rightarrow#tau#tau",Sample::Type::BG,kGreen)+TRegexp("SKDY50plus"));
  samples["vv"]=Sample("Diboson",Sample::Type::BG,kBlue)+TRegexp("SK[WZ][WZ]_pythia8");
  samples["wjets"]=Sample("W",Sample::Type::BG,kYellow)+TRegexp("SKWJets");
  samples["tt"]=Sample("t#bar{t}",Sample::Type::BG,kMagenta)+TRegexp("SKTTToLLNuNu");

  samples["sim_stack"]=Sample("sim",Sample::Type::STACK,kRed)+"amc"+"amctt"+"vv"+"wjets"+"tt";
  samples["sim"]=Sample("sim",Sample::Type::SUM,kRed)+"amc"+"amctt"+"vv"+"wjets"+"tt";
  samples["sim_noSF"]=(Sample("sim",Sample::Type::SUM,kBlue)+"amc"+"amctt"+"vv"+"wjets"+"tt")%"_noefficiencySF";

  Setup(mode_);
}

int LQ766Plotter::Setup(TString mode_){
  Reset();

  mode=mode_;

  SetupEntries();
  SetupSystematics();
  SetupPlots("plot_EfficiencyValidation/"+mode+"/plot.dat");

  if(DEBUG) std::cout<<"[Setup] nentries: "<<entries.size()<<endl;
  if(DEBUG) std::cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) std::cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}

void LQ766Plotter::SetupEntries(){
  if(DEBUG)  cout<<"[LQ766Plotter::SetupEntries] mode="<<mode<<endl;
  vector<TString> entry_keys=Split(mode," ");
  for(auto entry_key:entry_keys){
    if(samples.find(entry_key)!=samples.end())
      entries.push_back(samples[entry_key]);
    else{
      if(DEBUG>0) std::cout<<"###ERROR### [LQ766Plotter::SetupEntries] No "<<entry_key<<" in samples"<<endl;
    }
  }
  if(DEBUG>1) PrintEntries();
  return;
}
void LQ766Plotter::SetupSystematics(){
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

double LQ766Plotter::GetChi2(TH1* h1,TH1* h2){
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
