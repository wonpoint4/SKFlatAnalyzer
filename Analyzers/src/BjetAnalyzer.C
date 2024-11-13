#include "BjetAnalyzer.h"

BjetAnalyzer::BjetAnalyzer(){
}
BjetAnalyzer::~BjetAnalyzer(){
}
void BjetAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void BjetAnalyzer::executeEvent(){
  for(TString option:{"","AFBCharge_old","AFBCharge_no_pperpcut","AFBCharge_no_isocut","AFBCharge_no_ipcut","AFBCharge_no_drcut","AFBCharge_no_3chargecut"}){
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
      executeEventWithParameter(MakeParameter("mm",option));
    }
    if(!IsDATA||DataStream.Contains("SingleMuon")){
      executeEventWithParameter(MakeParameter("me",option));
    }
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
      executeEventWithParameter(MakeParameter("ee",option));
    }
  }
}
SMPAnalyzerCore::Parameter BjetAnalyzer::MakeParameter(TString key,TString option){
  Parameter p=SMPAnalyzerCore::MakeParameter(key,option);
  for(TString this_option:{"AFBCharge_old","AFBCharge_no_pperpcut","AFBCharge_no_isocut","AFBCharge_no_ipcut","AFBCharge_no_drcut","AFBCharge_no_3chargecut"}){
    if(option.Contains(this_option)){
      p.prefix=this_option(10,999)+"/"+p.prefix;
    }
  }
  return p;
}
void BjetAnalyzer::EvalDefaultWeight(Parameter& p){
  SMPAnalyzerCore::EvalDefaultWeight(p);
  p.default_weight/=p.w.bchargeSF;
  p.weight=p.default_weight;
}
SMPAnalyzerCore::Variations BjetAnalyzer::MakeVariations(const Parameter& p){
  Variations v;
  AddVariationWeight(v,"",p.default_weight);
  return v;
}
void BjetAnalyzer::FillHists(Parameter& p){
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double weight=p.weight;

  int nbjet=count_if(p.bjets.begin(),p.bjets.end(),[&p](auto& jet){return jet.Pt()>p.c.jetpt;});
  if(nbjet<1) return;

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass<52) return;
  FillHist(pre+"dimass"+suf,dimass,weight,nmassbin,massbins);

  double bjet_eta=p.bjets.at(0).Eta();
  double bjet_pt=p.bjets.at(0).Pt();
  FillHist(pre+"bjetapt"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
  int flavour=p.bjets.at(0).GenHFHadronMatcherFlavour();
  int origin=p.bjets.at(0).GenHFHadronMatcherOrigin();
  double AFBCharge=p.bjets.at(0).userFloat["AFBCharge"];
  if(flavour==5 && abs(origin)==6){
    int correct=AFBCharge*origin<0;
    TString scharge=origin<0 ? "p" : "m";
    int type=0;
    if(fabs(AFBCharge)<1) type=0;
    else if(fabs(AFBCharge)<3) type=1;
    else type=2;
    
    FillHist(pre+Form("b%setapt_correct%d_type%d",scharge.Data(),correct,type)+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
    FillHist(pre+Form("b%scharge_correct%d_type%d",scharge.Data(),correct,type)+suf,AFBCharge,weight,100,-5,5);
    if(type>0){
      FillHist(pre+Form("b%slpt_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_pt"],weight,100,0,100);
      FillHist(pre+Form("b%slpperp_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_pperp"],weight,100,0,10);
      FillHist(pre+Form("b%slriso_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_riso"],weight,100,0,0.2);
      FillHist(pre+Form("b%slsip_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_sip"],weight,100,0,10);
      FillHist(pre+Form("b%sldr_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_dr"],weight,80,0,0.8);
      FillHist(pre+Form("b%sl3charge_correct%d_type%d",scharge.Data(),correct,type)+suf,p.bjets.at(0).userFloat["AFBCharge_softlepton_3charge"],weight,2,0,2);      
    }
  }
}
