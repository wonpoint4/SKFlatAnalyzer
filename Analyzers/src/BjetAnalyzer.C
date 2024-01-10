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
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
  }
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

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  if(dimass<52) return;
  FillHist(pre+"dimass"+suf,dimass,weight,nmassbin,massbins);

  if(p.bjets.size()&&p.bjets.at(0).Pt()>p.c.jetpt){
    double bjet_eta=p.bjets.at(0).Eta();
    double bjet_pt=p.bjets.at(0).Pt();
    FillHist(pre+"bjetapt"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
    int flavour=p.bjets.at(0).GenHFHadronMatcherFlavour();
    int origin=p.bjets.at(0).GenHFHadronMatcherOrigin();
    if(flavour==-999){
      FillHist(pre+"bjetapt_p"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);      
    }else if(flavour==5){
      if(origin==-999){
	FillHist(pre+"bjetapt_n"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);	
      }else if(origin*p.bjets.at(0).userFloat["AFBCharge"]<0){
	FillHist(pre+"bjetapt_o"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
      }else{
	FillHist(pre+"bjetapt_x"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);
      }
    }else{
      FillHist(pre+"bjetapt_l"+suf,bjet_eta,bjet_pt,weight,netabin,etabins,nptbin,ptbins);      
    }
  }
}
