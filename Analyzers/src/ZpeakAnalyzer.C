#include "ZpeakAnalyzer.h"

ZpeakAnalyzer::ZpeakAnalyzer(){
}
ZpeakAnalyzer::~ZpeakAnalyzer(){
}
void ZpeakAnalyzer::initializeAnalyzer(){
  SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
  fChain->SetBranchStatus("pfMET_*",false);
  fChain->SetBranchStatus("pfMET_Type1_pt",true);
  fChain->SetBranchStatus("jet_*",false);
  fChain->SetBranchStatus("fatjet_*",false);
  fChain->SetBranchStatus("photon_*",false);
}
void ZpeakAnalyzer::executeEvent(){
  if(!IsDATA||DataStream.Contains("DoubleMuon")){
    executeEventWithParameter(MakeParameter("mm"));
    executeEventWithParameter(MakeParameter("MM"));
  }
  if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
    executeEventWithParameter(MakeParameter("ee"));
    executeEventWithParameter(MakeParameter("EE"));
  }
}

SMPAnalyzerCore::Variations ZpeakAnalyzer::MakeVariations(const Parameter& p){
  Variations v;
  AddVariationWeight(v,"",p.default_weight);
  if(!IsDATA){
    EvalVariationsCF(p,v);
  }
  return v;
}
void ZpeakAnalyzer::FillHists(Parameter& p){
  TString pre=p.prefix+p.hprefix;
  TString suf=p.suffix+p.vsuffix;
  double weight=p.weight;

  TLorentzVector dilepton=(*p.lepton0)+(*p.lepton1);
  double dimass=dilepton.M();
  FillHist(pre+"dimass"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
  FillHist(pre+"dimass"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
  TH1* h=GetHist3D(pre+"dimass"+suf);
  int l0ptbin=h->GetXaxis()->FindBin(p.lepton0->Pt());
  int l0etabin=h->GetYaxis()->FindBin(p.lepton0->Eta());
  int l1ptbin=h->GetXaxis()->FindBin(p.lepton1->Pt());
  int l1etabin=h->GetYaxis()->FindBin(p.lepton1->Eta());
  if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0){
    FillHist(pre+"dimass_even"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    FillHist(pre+"dimass_even"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
  }else{
    FillHist(pre+"dimass_odd"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    FillHist(pre+"dimass_odd"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
  }
  if(!IsDATA){
    if(p.lepton0->Charge()*p.truth_lepton0.Charge()<0){
      FillHist(pre+"dimass_cf"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0)
	FillHist(pre+"dimass_cf_even"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      else
	FillHist(pre+"dimass_cf_odd"+suf,p.lepton0->Pt(),p.lepton0->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    }
    if(p.lepton1->Charge()*p.truth_lepton1.Charge()<0){
      FillHist(pre+"dimass_cf"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      if((l0ptbin-l1ptbin+l0etabin-l1etabin)%2==0)
	FillHist(pre+"dimass_cf_even"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
      else
	FillHist(pre+"dimass_cf_odd"+suf,p.lepton1->Pt(),p.lepton1->Eta(),dimass,weight,nptbin,ptbins,netabin,etabins,nmassbin,massbins);
    }
  }
}
