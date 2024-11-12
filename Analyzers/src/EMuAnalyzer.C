#include "EMuAnalyzer.h"

EMuAnalyzer::EMuAnalyzer(){
}
EMuAnalyzer::~EMuAnalyzer(){
}
void EMuAnalyzer::executeEvent(){
  //////// Gen level ////////
  if(GetSkimName()==""){
    if(MCSample.Contains("TTLL")){
      FillHist("normcheck",0,reductionweight*MCweight()*_event.GetTriggerLumi("Full"),5,0,5);
      FillHist("normcheck",1,reductionweight*MCweight()*_event.GetTriggerLumi("Full")*mcCorr->GetTopPtReweight(gens),5,0,5);
      FillHist("normcheck",2,reductionweight*MCweight()*_event.GetTriggerLumi("Full")*GetTopPtReweight2(gens),5,0,5);
    }
    return;
  }
  
  //////// nominal channels //////////
  if(!IsDATA||DataStream.Contains("SingleMuon")){
    executeEventWithParameter(MakeParameter("me"));
    //executeEventWithParameter(MakeParameter("me","DeepJet::Medium"));
    //executeEventWithParameter(MakeParameter("me","DeepJet::Tight::mujets"));
    //executeEventWithParameter(MakeParameter("me","DeepCSV::Medium"));
    //executeEventWithParameter(MakeParameter("me","DeepCSV"));
  }
  if(!IsDATA||DataStream.Contains("SingleElectron")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("em"));
    //executeEventWithParameter(MakeParameter("em","DeepJet::Medium"));
    //executeEventWithParameter(MakeParameter("em","DeepJet::Tight::mujets"));
    //executeEventWithParameter(MakeParameter("em","DeepCSV::Medium"));
    //executeEventWithParameter(MakeParameter("em","DeepCSV"));
  }
}
SMPAnalyzerCore::Parameter EMuAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  if(p.suffix==""){
    p.variationbits|=NominalWeight|SystematicWeight|EfficiencyWeight;
    if(IsDYSample||IsTTLLSample){
      p.variationbits|=PDFWeight;
    }
  }
  if(option.Contains("DeepJet::Medium")){
    p.variationbits=NominalWeight;
    p.prefix="medium/"+p.prefix;
  }else if(option.Contains("DeepJet::Tight::mujets")){
    p.variationbits=NominalWeight;
    p.prefix="mujets/"+p.prefix;
  }else if(option.Contains("DeepCSV::Medium")){
    p.variationbits=NominalWeight;
    p.prefix="csvmedium/"+p.prefix;
  }else if(option.Contains("DeepCSV")){
    p.variationbits=NominalWeight;
    p.prefix="csv/"+p.prefix;
  }
  return p;
}
void EMuAnalyzer::FillHists(Parameter& p){
  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass<52) return;
    
  TString region=p.bjets.size()&&p.bjets.at(0).Pt()>p.c.jetpt ? "nbjet/" : "0bjet/";
  TString pre=p.prefix+region+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double w=p.weight;
    
  //for leptons
  for(int i=0;i<(int)p.leptons.size();i++){
    double pt=p.leptons.at(i)->Pt();
    double eta=p.leptons.at(i)->Eta();
    FillHist(Form("%sl%dpt%s",pre.Data(),i,suf.Data()),pt,w,nptbin,ptbins);
    FillHist(Form("%sl%deta%s",pre.Data(),i,suf.Data()),eta,w,50,-2.5,2.5);
    
    FillHist(Form("%slpt%s",pre.Data(),suf.Data()),pt,w,nptbin,ptbins);
    FillHist(Form("%sleta%s",pre.Data(),suf.Data()),eta,w,50,-2.5,2.5);
    
    if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::MUON){
      FillHist(Form("%smpt%s",pre.Data(),suf.Data()),pt,w,nptbin,ptbins);
      FillHist(Form("%smeta%s",pre.Data(),suf.Data()),eta,w,50,-2.5,2.5);
    }else if(p.leptons.at(i)->LeptonFlavour()==Lepton::Flavour::ELECTRON){
      Electron* el=(Electron*)p.leptons.at(i);
      FillHist(Form("%slsceta%s",pre.Data(),suf.Data()),el->scEta(),w,50,-2.5,2.5);
      FillHist(Form("%sept%s",pre.Data(),suf.Data()),pt,w,nptbin,ptbins);
      FillHist(Form("%seeta%s",pre.Data(),suf.Data()),eta,w,50,-2.5,2.5);
    }
  }
  
  double dipt=dilepton.Pt();
  double dirap=dilepton.Rapidity();
  FillHist(pre+"dimass"+suf,dimass,w,nmbin,mbins);
  FillHist(pre+"dipt"+suf,dipt,w,nptbin,ptbins);
  FillHist(pre+"dirap"+suf,dirap,w,50,-2.5,2.5);

}
