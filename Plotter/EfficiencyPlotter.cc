#include"Plotter.cc"
#ifndef __EFFICIENCYPLOTTER_CC__
#define __EFFICIENCYPLOTTER_CC__
class EfficiencyPlotter:public Plotter{
public:
  void SetupSystematics();
  int Setup(TString mode_);
  TString mode;
  EfficiencyPlotter(TString mode_="data sim_stack");
};
EfficiencyPlotter::EfficiencyPlotter(TString mode_){
  ScanFiles((TString)getenv("SKFlatOutputDir")+getenv("SKFlatV")+"/EfficiencyValidation/");

  samples["muon"]=Sample("data (ee)",Sample::Type::DATA,kBlack,20)+TRegexp("/DATA/EfficiencyValidation_DoubleMuon_[A-Z]")+TRegexp("/DATA/EfficiencyValidation_SingleMuon_[A-Z]");
  samples["electron"]=Sample("data (#mu#mu)",Sample::Type::DATA,kBlack,20)+TRegexp("/EfficiencyValidation_.*EG.*_[A-Z]")+TRegexp("/EfficiencyValidation_SingleElectron_[A-Z]");
  samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+"muon"+"electron";
  samples["amc"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/EfficiencyValidation_DYJets$");
  samples["amctt"]="tau_"%(Sample("#gamma*/Z#rightarrow#tau#tau",Sample::Type::BG,kGreen)+TRegexp("/EfficiencyValidation_DYJets$"));
  samples["vv"]=Sample("Diboson",Sample::Type::BG,kBlue)+TRegexp("/EfficiencyValidation_[W-Z][W-Z]_pythia$");
  samples["wjets"]=Sample("W",Sample::Type::BG,kYellow)+TRegexp("/EfficiencyValidation_WJets_MG$");
  samples["tt"]=Sample("t#bar{t}",Sample::Type::BG,kMagenta)+TRegexp("/EfficiencyValidation_TTLL_powheg$");

  samples["sim_stack"]=Sample("sim",Sample::Type::STACK,Style(kRed,-1,3001,"e2"),Style(kCyan,-1,3001,"e2"))+"amc"+"amctt"+"vv"+"wjets"+"tt";
  samples["sim"]=Sample("simulation",Sample::Type::SUM,Style(kRed,22,3001,"e2"),Style(kCyan,-1,3001,"e2"))+"amc"+"amctt"+"vv"+"wjets"+"tt";
  samples["sim_noSF"]=(Sample("w/o efficiency SF",Sample::Type::SUM,kBlue)+"amc"+"amctt"+"vv"+"wjets"+"tt")%"_noefficiencySF";
  for(auto& sub:samples["sim_noSF"].subs) sub.type=Sample::Type::A;
  samples["sim_noSF"].style.linewidth=1;

  Setup(mode_);
}

int EfficiencyPlotter::Setup(TString mode_){
  Reset();

  mode=mode_;

  SetupEntries(mode_);
  SetupSystematics();
  SetupPlots("plot_EfficiencyValidation/"+mode+"/plot.dat");

  if(DEBUG) std::cout<<"[Setup] nentries: "<<entries.size()<<endl;
  if(DEBUG) std::cout<<"[Setup] nsys: "<<systematics.size()<<endl;
  if(DEBUG) std::cout<<"[Setup] nplot: "<<plots.size()<<endl;

  return 1;
}
void EfficiencyPlotter::SetupSystematics(){
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
#endif
